"""Агент докинга: считает AutoDock Vina для любой пары «структура + лиганд».

Сохранённый прогон используется только как запасной путь и только при точном
совпадении пары (PDB, лиганд). Чужой скоринг под видом ответа не подставляется.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

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
    try:
        result = dock_pair(str(pdb_id), str(smiles), ligand_id, output_dir)
    except Exception as exc:  # noqa: BLE001 — граница внешнего бинарника
        state.limitations.append(f"Живой AutoDock Vina не завершился: {exc}")
        return False

    state.artifacts["poses_csv"] = str(result.poses_csv)
    state.artifacts["pose_pdbqt"] = str(result.pose_file)
    state.artifacts["receptor_pdbqt"] = str(result.receptor_pdbqt)
    state.observations["poses_csv_path"] = str(result.poses_csv)
    state.observations["docking_mode"] = "live"
    _record_score(
        state,
        result.vina_score,
        ligand_id,
        "site1",
        str(result.pose_file),
        (
            "Оценка Vina посчитана в этом запуске. Это не экспериментальная аффинность. "
            f"Бокс: {result.box_source}. Подготовка лиганда: {result.ligand_prep}."
        ),
    )
    state.limitations.append(
        "Рецептор переведён в PDBQT из координат ATOM без отдельной оптимизации зарядов. "
        "Карман взят по лиганду кристалла, если он есть, иначе по центру белка."
    )
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
    best = min(rows, key=lambda r: r.vina_score)
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
