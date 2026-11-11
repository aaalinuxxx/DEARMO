"""Агент докинга: загрузка кеша или live (DOCKING_LIVE=1)."""
from __future__ import annotations

import os
from pathlib import Path

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.docking_io import load_run_directory, poses_from_dir, read_report_summary


def run_or_load(state: WorkflowState, example_dir: Path) -> WorkflowState:
    with state.agent_step("docking_agent", "load_cached_or_run"):
        example_dir = Path(example_dir)
        live = os.environ.get("DOCKING_LIVE", "").lower() in ("1", "true", "yes")
        if live:
            state.limitations.append("DOCKING_LIVE запрошен, но полный live-пайплайн требует PDBQT; используйте docking_pipeline CLI.")
        if not example_dir.is_dir():
            state.add_error(f"Нет каталога примера: {example_dir}")
            return state

        artifacts = load_run_directory(example_dir)
        state.artifacts.update(artifacts)
        if "report_md" in artifacts:
            summary = read_report_summary(Path(artifacts["report_md"]))
            state.observations["docking_report_excerpt"] = summary[:800]

        try:
            rows = poses_from_dir(example_dir)
        except FileNotFoundError:
            state.add_error("poses.csv не найден в примере.")
            return state

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
                    "note": "Оценка Vina — расчёт, не экспериментальная аффинность.",
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
            "Докинг загружен из сохранённого прогона (data/examples/). RMSD к первой позе Vina и к кристаллу подписываются отдельно в poses.csv."
        )
    return state
