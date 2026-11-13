"""Докинг «вслепую»: PDB + SMILES → карманы P2Rank → Vina в каждый карман.

Сценарий 2 целиком. Ключевое отличие от прежней версии — центр бокса больше не
берётся по кристаллическому лиганду. Раньше система получала координаты ответа
на вход и считала один бокс вокруг них, поэтому «карман найден» ничего не
значило: его не искали. Теперь карманы предсказывает P2Rank по очищенной
структуре, в которой лиганда нет, а кристаллический лиганд используется только
после расчёта — как эталон для контрольных метрик.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from models.docking import (
    PocketRow,
    PoseRow,
    best_pose_overall,
    write_pockets_csv,
    write_poses_csv,
)
from tools.p2rank_runner import predict_pockets
from tools.pdbqt_prep import (
    clean_receptor_pdb,
    crystal_ligand_center,
    ligand_pdbqt_from_smiles,
    protein_center,
    protonated_receptor_pdbqt,
    receptor_pdbqt,
    write_text,
)
from tools.rmsd import RmsdError, min_distance_to_atoms, rmsd_to_crystal, select_crystal_reference
from tools.vina_runner import run_vina_modes, split_pdbqt_models, vina_binary

_ROOT = Path(__file__).resolve().parent.parent

# Размер бокса вокруг центра кармана. 22 Å покрывает типичный карман вместе с
# подступами к нему и при этом не превращает докинг в перебор по всему белку.
DEFAULT_BOX = (22.0, 22.0, 22.0)

# Общепринятый порог успешного редокинга: поза №1 не дальше 2 Å от кристаллической.
RMSD_SUCCESS_A = 2.0
# Порог контрольной метрики защиты: карман считается верным, если его центр
# лежит не дальше этого расстояния от ближайшего тяжёлого атома кристалла.
POCKET_HIT_RADIUS_A = 4.0


@dataclass
class DockingRun:
    """Результат одного прогона сценария 2."""

    pdb_id: str
    ligand_id: str
    pdb_path: Path
    receptor_pdb: Path
    receptor_pdbqt: Path
    ligand_pdbqt: Path
    pockets: list[PocketRow]
    poses: list[PoseRow]
    best: PoseRow
    poses_csv: Path
    pockets_csv: Path
    run_meta: Path
    box_source: str
    box_size: tuple[float, float, float]
    ligand_prep: str
    exhaustiveness: int
    num_modes: int
    seed: int | None
    control: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def vina_score(self) -> float:
        return self.best.vina_score

    @property
    def pose_file(self) -> Path:
        return Path(self.best.pose_file)

    @property
    def center(self) -> tuple[float, float, float]:
        return self.best.pocket_center or (0.0, 0.0, 0.0)


def locate_pdb(pdb_id: str) -> Path:
    clean = pdb_id.strip().upper()
    if len(clean) != 4:
        path = Path(pdb_id)
        if path.is_file():
            return path
        raise FileNotFoundError(f"Нет PDB-файла: {pdb_id}")
    for folder in (_ROOT / "data" / "external" / "rcsb", _ROOT / "data" / "cache" / "pdb"):
        candidate = folder / f"{clean}.pdb"
        if candidate.is_file():
            return candidate
    dest = _ROOT / "data" / "cache" / "pdb" / f"{clean}.pdb"
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{clean}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "DEARMO/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        dest.write_bytes(resp.read())
    return dest


def _fallback_pockets(pdb_text: str, reason: str) -> tuple[list[PocketRow], str, list[str]]:
    """Один бокс по центру белка, когда карманы предсказать не удалось.

    Это не поиск карманов, а честное «ищем по всему белку». Помечается
    отдельным box_source, чтобы такой прогон нельзя было выдать за результат
    P2Rank на защите.
    """
    center = protein_center(pdb_text)
    pocket = PocketRow(pocket_id="whole_protein", rank=1, p2rank_score=None, center=center)
    return [pocket], "protein_center", [f"Карманы не предсказаны ({reason}). Бокс поставлен по центру белка."]


def dock_pockets(
    pdb_id: str,
    smiles: str,
    ligand_id: str,
    output_dir: Path,
    *,
    chain: str | None = None,
    top_n: int = 3,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int | None = 42,
    box_size: tuple[float, float, float] = DEFAULT_BOX,
    crystal_comp_id: str | None = None,
) -> DockingRun:
    """Найти карманы, посчитать докинг в каждом и собрать таблицу поз.

    crystal_comp_id не влияет на расчёт: он нужен только для контрольных
    метрик после него. Докинг запускается одинаково независимо от того,
    известен нам ответ или нет.
    """
    if vina_binary() is None:
        raise RuntimeError(
            "AutoDock Vina не найден. Положите бинарник в tools/bin/vina или задайте VINA_BIN."
        )
    pdb_path = locate_pdb(pdb_id)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    notes: list[str] = []

    receptor_pdb_path = write_text(output_dir / "receptor.pdb", clean_receptor_pdb(pdb_text, chain))
    try:
        receptor_text, prep_notes = protonated_receptor_pdbqt(receptor_pdb_path, output_dir / "receptor_prep")
        receptor_prep = "pdb2pqr+meeko"
        notes.extend(prep_notes)
    except Exception as exc:  # noqa: BLE001 — граница внешних программ
        # Откат на грубую подготовку возможен, но молча делать его нельзя: без
        # водородов и корректных типов атомов оценка Vina и позы заметно хуже.
        receptor_text = receptor_pdbqt(pdb_text, chain)
        receptor_prep = "coordinates-only"
        notes.append(
            "Рецептор подготовлен упрощённо, без водородов и типизации Meeko "
            f"({str(exc)[:160]}). Это занижает качество скоринга: у белка нет доноров "
            "водородных связей, а RMSD до кристалла может оказаться завышенным."
        )
    receptor_path = write_text(output_dir / "receptor.pdbqt", receptor_text)
    ligand_text, ligand_prep = ligand_pdbqt_from_smiles(smiles)
    ligand_path = write_text(output_dir / "ligand.pdbqt", ligand_text)
    if ligand_prep.startswith("rdkit-rigid"):
        # Жёсткий лиганд не меняет торсионы, поэтому гибкая молекула докируется
        # в одной конформации. Тихо это делать нельзя: низкий RMSD или его
        # отсутствие нельзя будет честно объяснить.
        notes.append(
            "Meeko не подготовил лиганд, он докирован как жёсткое тело (TORSDOF 0): "
            "гибкие торсионы не варьировались, оценка и RMSD занижают возможности Vina."
        )

    try:
        pockets = predict_pockets(receptor_pdb_path, output_dir / "p2rank", top_n=top_n)
        box_source = "p2rank"
    except Exception as exc:  # noqa: BLE001 — граница внешнего бинарника
        pockets, box_source, extra = _fallback_pockets(pdb_text, str(exc)[:200])
        notes.extend(extra)

    poses_dir = output_dir / "poses"
    poses: list[PoseRow] = []
    for pocket in pockets:
        multi = output_dir / f"{pocket.pocket_id}_out.pdbqt"
        try:
            modes, multi = run_vina_modes(
                receptor_path,
                ligand_path,
                pocket.center,
                size=box_size,
                out_pose=multi,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                seed=seed,
            )
        except Exception as exc:  # noqa: BLE001
            notes.append(f"Карман {pocket.pocket_id}: Vina не дала поз ({str(exc)[:160]}).")
            continue
        files = split_pdbqt_models(multi, poses_dir, pocket.pocket_id)
        for mode in modes:
            # Если файлов меньше, чем строк таблицы, позу без координат не
            # записываем: строка с оценкой, но без файла, не проверяема.
            if mode.mode > len(files):
                notes.append(
                    f"Карман {pocket.pocket_id}: поза {mode.mode} есть в таблице Vina, но не в файле координат."
                )
                continue
            poses.append(
                PoseRow(
                    ligand_id=ligand_id,
                    pocket_id=pocket.pocket_id,
                    pocket_rank=pocket.rank,
                    p2rank_score=pocket.p2rank_score,
                    pose_id=mode.mode,
                    vina_score=mode.affinity,
                    # Абсолютный путь: таблицу читают другие сценарии из своих каталогов.
                    pose_file=str(files[mode.mode - 1].resolve()),
                    rmsd_lb=mode.rmsd_lb,
                    rmsd_ub=mode.rmsd_ub,
                    rmsd_to_first_pose=mode.rmsd_lb,
                    pocket_center=pocket.center,
                )
            )

    if not poses:
        raise RuntimeError(
            f"Докинг не дал ни одной позы ни в одном из {len(pockets)} карманов. " + " ".join(notes)
        )

    best = best_pose_overall(poses)
    assert best is not None  # poses непуст, проверено выше

    control = _control_metrics(poses, pockets, best, pdb_text, crystal_comp_id, chain, notes)

    poses_csv = write_poses_csv(output_dir / "poses.csv", poses)
    pockets_csv = write_pockets_csv(output_dir / "pockets.csv", pockets)

    meta = {
        "pdb_id": str(pdb_id).strip().upper(),
        "pdb_chain": chain,
        "ligand_id": ligand_id,
        "smiles": smiles,
        "engine": "AutoDock Vina",
        "scoring": "vina",
        "exhaustiveness": exhaustiveness,
        "num_modes": num_modes,
        "seed": seed,
        "box_size": list(box_size),
        "box_source": box_source,
        "pocket_source": box_source,
        "ligand_prep": ligand_prep,
        "receptor_prep": receptor_prep,
        # Абсолютные пути: паспорт прогона читают сценарии 3 и 4 из своего
        # рабочего каталога, и относительный путь у них бы не разрешился.
        "receptor_pdb_path": str(receptor_pdb_path.resolve()),
        "receptor_pdbqt": str(receptor_path.resolve()),
        "pdb_path": str(pdb_path.resolve()),
        "pocket_count": len(pockets),
        "pose_count": len(poses),
        "best_pose": {
            "pocket_id": best.pocket_id,
            "pose_id": best.pose_id,
            "vina_score": best.vina_score,
            "pose_file": best.pose_file,
        },
        "control": control,
        "notes": notes,
    }
    meta_path = output_dir / "run_meta.json"
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    return DockingRun(
        pdb_id=str(pdb_id).strip().upper(),
        ligand_id=ligand_id,
        pdb_path=pdb_path,
        receptor_pdb=receptor_pdb_path,
        receptor_pdbqt=receptor_path,
        ligand_pdbqt=ligand_path,
        pockets=pockets,
        poses=poses,
        best=best,
        poses_csv=poses_csv,
        pockets_csv=pockets_csv,
        run_meta=meta_path,
        box_source=box_source,
        box_size=box_size,
        ligand_prep=ligand_prep,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        control=control,
        notes=notes,
    )


def _control_metrics(
    poses: list[PoseRow],
    pockets: list[PocketRow],
    best: PoseRow,
    pdb_text: str,
    crystal_comp_id: str | None,
    chain: str | None,
    notes: list[str],
) -> dict:
    """Контрольные метрики защиты. Считаются только при наличии эталона.

    Если кристаллического лиганда нет, метрики отсутствуют. Подставлять вместо
    них ноль или «успех по умолчанию» нельзя: отсутствие проверки — это не
    пройденная проверка.
    """
    if not crystal_comp_id:
        return {"available": False, "reason": "Эталонный лиганд не задан"}
    try:
        reference, reference_note = select_crystal_reference(pdb_text, crystal_comp_id, chain)
    except RmsdError as exc:
        notes.append(f"Контрольные метрики не посчитаны: {exc}")
        return {"available": False, "reason": str(exc)}

    # Ранг первого кармана, попавшего в окрестность кристаллического лиганда.
    hit_rank: int | None = None
    hit_distance: float | None = None
    for pocket in sorted(pockets, key=lambda p: p.rank):
        distance = min_distance_to_atoms(pocket.center, reference)
        if distance <= POCKET_HIT_RADIUS_A:
            hit_rank, hit_distance = pocket.rank, distance
            break

    best_rmsd: float | None = None
    try:
        best_rmsd = rmsd_to_crystal(best.pose_file, pdb_text, crystal_comp_id, chain)
    except RmsdError as exc:
        notes.append(f"RMSD позы №1 до кристалла не посчитан: {exc}")

    # RMSD до кристалла заполняется для всех поз, а не только для позы №1:
    # иначе нельзя показать, что лучшая по Vina поза — действительно лучшая
    # и геометрически, а не просто по оценке.
    for pose in poses:
        try:
            pose.rmsd_to_crystal = rmsd_to_crystal(pose.pose_file, pdb_text, crystal_comp_id, chain)
        except RmsdError:
            pose.rmsd_to_crystal = None

    passed = hit_rank is not None and hit_rank <= 3 and best_rmsd is not None and best_rmsd <= RMSD_SUCCESS_A
    return {
        "available": True,
        "crystal_comp_id": crystal_comp_id,
        "reference": reference_note,
        "pocket_hit_rank": hit_rank,
        "pocket_hit_distance_a": hit_distance,
        "pocket_hit_radius_a": POCKET_HIT_RADIUS_A,
        "best_pose_rmsd_to_crystal_a": best_rmsd,
        "rmsd_threshold_a": RMSD_SUCCESS_A,
        "demo_check_passed": passed,
    }


def dock_pair(
    pdb_id: str,
    smiles: str,
    ligand_id: str,
    output_dir: Path,
    *,
    exhaustiveness: int = 8,
    **kwargs,
) -> DockingRun:
    """Прежнее имя точки входа. Теперь считает карманы, а не один бокс."""
    return dock_pockets(
        pdb_id,
        smiles,
        ligand_id,
        output_dir,
        exhaustiveness=exhaustiveness,
        **kwargs,
    )


# Прежнее имя типа результата: часть кода и тестов ссылается на него.
LiveDockResult = DockingRun

__all__ = [
    "DEFAULT_BOX",
    "POCKET_HIT_RADIUS_A",
    "DockingRun",
    "LiveDockResult",
    "crystal_ligand_center",
    "dock_pair",
    "dock_pockets",
    "locate_pdb",
]
