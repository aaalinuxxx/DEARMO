"""Агент докинга: считает AutoDock Vina для любой пары «структура + лиганд».

Сохранённый прогон используется только как запасной путь и только при точном
совпадении пары (PDB, лиганд). Чужой скоринг под видом ответа не подставляется.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

from models.docking import best_pose_overall
from models.result import Evidence
from orchestration.state import WorkflowState
from tools.docking_io import load_run_directory, poses_from_dir, read_report_summary
from tools.live_docking import dock_pair

_ROOT = Path(__file__).resolve().parent.parent
_EXAMPLES = _ROOT / "data" / "examples"


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_") or "ligand"


def run_key(pdb_id: str | None, ligand_id: str | None) -> str:
    return f"{str(pdb_id or 'nopdb').upper()}_{_slug(ligand_id or 'ligand')}"


def _saved_run_for(pdb_id: str | None, ligand_id: str | None) -> Path | None:
    """Каталог сохранённого прогона ровно этой пары, если он есть.

    Пара подтверждается паспортом run_meta.json внутри каталога: совпасть должны
    и структура, и лиганд. Каталог без паспорта не используется — иначе нельзя
    доказать, что лежащие в нём позы относятся именно к этому запросу.
    """
    if not pdb_id or not ligand_id:
        return None
    wanted_pdb = str(pdb_id).strip().upper()
    wanted_ligand = _slug(ligand_id)
    for base in (_EXAMPLES, _ROOT / "outputs" / "docking"):
        if not base.is_dir():
            continue
        for candidate in sorted(base.iterdir()):
            meta_path = candidate / "run_meta.json"
            if not candidate.is_dir() or not (candidate / "poses.csv").is_file():
                continue
            if not meta_path.is_file():
                continue
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if str(meta.get("pdb_id", "")).strip().upper() != wanted_pdb:
                continue
            if _slug(meta.get("ligand_id", "")) != wanted_ligand:
                continue
            return candidate
    return None


def _live_disabled() -> bool:
    return os.environ.get("DOCKING_LIVE", "1").lower() in ("0", "false", "no")


def _record_score(
    state: WorkflowState,
    score: float,
    ligand_id: str,
    pocket_id: str,
    pose_file: str,
    note: str,
) -> None:
    state.observations["best_pose"] = pose_file
    state.observations["pose_file"] = pose_file
    state.add_evidence(
        Evidence(
            kind="calculation",
            value=score,
            unit="kcal/mol",
            method="vina_score",
            source="AutoDock Vina",
            context={
                "ligand_id": ligand_id,
                "pocket_id": pocket_id,
                "pose_file": pose_file,
                "note": note,
            },
        )
    )


def _first_comp_id(state: WorkflowState) -> str | None:
    """Код кристаллического лиганда, если сценарий 1 его нашёл.

    Нужен только для контрольных метрик после расчёта. На сам докинг не влияет:
    карманы ищутся по структуре, из которой лиганд вырезан.
    """
    # Код, подтверждённый составом именно этой структуры, надёжнее общего списка
    # компонентов, которые RCSB сопоставила молекуле.
    confirmed = state.observations.get("crystal_comp_id")
    if confirmed:
        return str(confirmed)
    comp_ids = state.observations.get("ligand_comp_ids")
    if isinstance(comp_ids, (list, tuple)) and comp_ids:
        return str(comp_ids[0])
    if isinstance(comp_ids, str) and comp_ids.strip():
        return comp_ids.strip()
    return None


def _record_pockets(state: WorkflowState, result) -> None:
    """Карманы-кандидаты: ранг, балл и центр каждого."""
    for pocket in result.pockets:
        state.add_evidence(
            Evidence(
                # Балл P2Rank — выход обученной модели машинного обучения, то есть
                # предсказание, а не расчёт по формуле. Кейс требует их различать.
                kind="prediction",
                value=pocket.p2rank_score,
                method="p2rank_pocket",
                source="P2Rank",
                context={
                    "note": "Предсказание модели машинного обучения P2Rank; не измерение.",
                    "pocket_id": pocket.pocket_id,
                    "rank": pocket.rank,
                    "center": [round(c, 3) for c in pocket.center],
                    "surf_atom_count": pocket.surf_atom_count,
                    "residue_ids": list(pocket.residue_ids[:12]),
                },
            )
        )


def _record_control(state: WorkflowState, result) -> None:
    """Контрольные метрики защиты, только если эталон действительно был."""
    control = result.control or {}
    if not control.get("available"):
        return
    if control.get("pocket_hit_rank") is not None:
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=control["pocket_hit_rank"],
                method="pocket_hit_rank",
                source="P2Rank",
                context={
                    "distance_a": round(float(control["pocket_hit_distance_a"]), 3),
                    "radius_a": control["pocket_hit_radius_a"],
                    "crystal_comp_id": control.get("crystal_comp_id"),
                    "note": "Ранг первого кармана, центр которого попал в окрестность кристаллического лиганда.",
                },
            )
        )
    rmsd = control.get("best_pose_rmsd_to_crystal_a")
    if rmsd is not None:
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=round(float(rmsd), 3),
                unit="Å",
                method="rmsd_to_crystal",
                source="DEARMO rmsd",
                context={
                    "pose_file": result.best.pose_file,
                    "crystal_comp_id": control.get("crystal_comp_id"),
                    "threshold_a": control.get("rmsd_threshold_a"),
                    "note": "RMSD позы №1 до кристаллической по тяжёлым атомам, с учётом симметрии, без выравнивания.",
                },
            )
        )
        # Голое число RMSD вердикта не даёт. Провал контроля должен быть сказан
        # словами: иначе ответ «success» с позой №1 выглядит как надёжный результат.
        threshold = float(control.get("rmsd_threshold_a") or 2.0)
        code = control.get("crystal_comp_id")
        if float(rmsd) > threshold:
            state.limitations.append(
                f"Контроль редокинга не пройден: поза №1 лежит в {float(rmsd):.1f} Å от кристаллической "
                f"({code}), порог {threshold:g} Å. Оценке Vina и положению позы №1 следует доверять "
                "с осторожностью: найденный карман верен, но скоринг мог поставить лиганд в иную ориентацию."
            )
            state.observations["redock_check"] = "failed"
        else:
            state.limitations.append(
                f"Контроль редокинга пройден: поза №1 в {float(rmsd):.1f} Å от кристаллической ({code}), "
                f"порог {threshold:g} Å. Это проверка метода на этой структуре, а не гарантия для других лигандов."
            )
            state.observations["redock_check"] = "passed"


def _run_live(state: WorkflowState) -> bool:
    """Запустить Vina для текущей пары. True, если скоринг получен в этом запуске."""
    pdb_id = state.observations.get("pdb_id") or state.request.pdb_id
    smiles = state.observations.get("smiles")
    ligand_id = str(
        state.observations.get("ligand_id") or state.observations.get("ligand_query") or "ligand"
    )
    if not pdb_id or not smiles:
        state.limitations.append(
            "Живой Vina не запущен: нужны PDB ID и SMILES лиганда. "
            "Чужой сохранённый скоринг не подставляется."
        )
        return False
    output_dir = _ROOT / "outputs" / "docking" / run_key(pdb_id, ligand_id)
    chain = state.observations.get("pdb_chain")
    try:
        result = dock_pair(
            str(pdb_id),
            str(smiles),
            ligand_id,
            output_dir,
            chain=str(chain) if chain else None,
            crystal_comp_id=_first_comp_id(state),
        )
    except Exception as exc:  # noqa: BLE001 — граница внешнего бинарника
        state.limitations.append(f"Живой AutoDock Vina не завершился: {exc}")
        return False

    state.artifacts["poses_csv"] = str(result.poses_csv)
    state.artifacts["pockets_csv"] = str(result.pockets_csv)
    state.artifacts["pose_pdbqt"] = str(result.pose_file)
    state.artifacts["receptor_pdbqt"] = str(result.receptor_pdbqt)
    state.artifacts["receptor_pdb"] = str(result.receptor_pdb)
    state.artifacts["run_meta"] = str(result.run_meta)
    state.observations["poses_csv_path"] = str(result.poses_csv)
    state.observations["pockets_csv_path"] = str(result.pockets_csv)
    state.observations["receptor_pdb_path"] = str(result.receptor_pdb)
    state.observations["docking_mode"] = "live"
    state.observations["receptor_prep"] = (
        json.loads(Path(result.run_meta).read_text(encoding="utf-8")).get("receptor_prep")
        if Path(result.run_meta).is_file()
        else None
    )
    state.observations["pocket_count"] = len(result.pockets)
    state.observations["pose_count"] = len(result.poses)

    _record_pockets(state, result)
    _record_score(
        state,
        result.best.vina_score,
        ligand_id,
        result.best.pocket_id,
        str(result.pose_file),
        (
            "Поза №1: самая отрицательная оценка Vina среди всех карманов этого прогона "
            f"({len(result.poses)} поз в {len(result.pockets)} карманах) при одинаковых настройках скоринга. "
            "Это расчёт, а не экспериментальная аффинность. "
            f"Источник карманов: {result.box_source}. Подготовка лиганда: {result.ligand_prep}."
        ),
    )
    _record_control(state, result)

    # Оговорка зависит от того, как рецептор был подготовлен на самом деле: при
    # протонировании pdb2pqr+Meeko она неверна, а при откате на координаты — нужна.
    # Об откате живой докинг уже сообщает отдельной заметкой.
    if state.observations.get("receptor_prep") == "coordinates-only":
        state.limitations.append(
            "Рецептор переведён в PDBQT из координат ATOM без водородов и типизации Meeko."
        )
    requested = int(state.observations.get("pocket_top_n") or 3)
    found = len(getattr(result, "pockets", []) or [])
    if 0 < found < requested:
        state.limitations.append(
            f"P2Rank нашёл {found} {'карман' if found == 1 else 'кармана'} из запрошенных {requested}: "
            "список не дополнялся выдуманными центрами, докинг выполнен только в найденных."
        )
    if result.box_source == "p2rank":
        state.limitations.append(
            "Докинг выполнен вслепую: центры боксов взяты у P2Rank по структуре без лиганда, "
            "координаты кристаллического лиганда в расчёт не подавались."
        )
    else:
        state.limitations.append(
            "Карманы предсказать не удалось, бокс поставлен по центру белка. "
            "Это поиск по всему белку, а не результат P2Rank."
        )
    state.limitations.extend(result.notes)
    return True


def _load_saved(state: WorkflowState, run_dir: Path) -> bool:
    """Прочитать сохранённый прогон ровно этой пары."""
    artifacts = load_run_directory(run_dir)
    state.artifacts.update(artifacts)
    if "report_md" in artifacts:
        state.observations["docking_report_excerpt"] = read_report_summary(
            Path(artifacts["report_md"])
        )[:800]
    try:
        rows = poses_from_dir(run_dir)
    except FileNotFoundError:
        state.add_error(f"poses.csv не найден в {run_dir}")
        return False
    if not rows:
        state.add_error(f"Пустая таблица поз в {run_dir}")
        return False

    state.observations["docking_mode"] = "saved"
    if "poses_csv" in artifacts:
        state.observations["poses_csv_path"] = artifacts["poses_csv"]
    if "receptor_pdb" in artifacts:
        state.observations["receptor_pdb_path"] = artifacts["receptor_pdb"]
    # Поза №1 выбирается тем же правилом, что и в живом расчёте: иначе
    # сохранённый и свежий прогоны называли бы «лучшей» разные позы.
    best = best_pose_overall(rows)
    if best is None:
        state.add_error(f"Пустая таблица поз в {run_dir}")
        return False
    state.observations["best_pose"] = best.pose_file
    state.observations["pose_file"] = best.pose_file
    state.add_evidence(
        Evidence(
            kind="calculation",
            value=best.vina_score,
            unit="kcal/mol",
            method="vina_score",
            source="AutoDock Vina",
            context={
                "ligand_id": best.ligand_id,
                "pocket_id": best.pocket_id,
                "pose_id": best.pose_id,
                "run_dir": str(run_dir),
                "note": "Оценка Vina из сохранённого прогона этой же пары, не экспериментальная аффинность.",
            },
        )
    )
    if best.rmsd_to_crystal is not None:
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=best.rmsd_to_crystal,
                unit="Å",
                method="rmsd_to_crystal",
                source="docking",
                context={"pose_file": best.pose_file},
            )
        )
    state.limitations.append(
        f"Докинг не считался в этом запуске: взят сохранённый прогон той же пары из {run_dir}."
    )
    return True


def run_or_load(state: WorkflowState, example_dir: Path | None = None) -> WorkflowState:
    """Сначала живой расчёт, при неудаче — сохранённый прогон ровно этой пары."""
    with state.agent_step("docking_agent", "vina"):
        pdb_id = state.observations.get("pdb_id") or state.request.pdb_id
        ligand_id = state.observations.get("ligand_id") or state.observations.get("ligand_query")

        if not _live_disabled() and _run_live(state):
            return state

        saved = _saved_run_for(pdb_id, ligand_id)
        if saved is not None and _load_saved(state, saved):
            return state

        state.add_error(
            f"Докинг не выполнен для пары {run_key(pdb_id, ligand_id)}: живой Vina не дал результата, "
            "а сохранённого прогона именно этой пары нет."
        )
        return state
