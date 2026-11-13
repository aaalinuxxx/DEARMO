"""Обработчики сценариев — детерминированные цепочки агентов."""
from __future__ import annotations

from pathlib import Path

from orchestration.state import WorkflowState
from orchestration.tool_router import run_selected_tools

_REPO_ROOT = Path(__file__).resolve().parent.parent
# Каталог сохранённых прогонов. Конкретная молекула по умолчанию не подставляется.
_EXAMPLES_DIR = _REPO_ROOT / "data" / "examples"


def handle_target_for_ligand(state: WorkflowState) -> WorkflowState:
    return run_selected_tools(state, "target_for_ligand")


def handle_dock_pockets(state: WorkflowState) -> WorkflowState:
    if not state.observations.get("pdb_id") and state.request.pdb_id:
        state.observations["pdb_id"] = state.request.pdb_id
    if not state.observations.get("ligand_id") and state.request.ligand_ids:
        state.observations["ligand_id"] = state.request.ligand_ids[0]
    return run_selected_tools(state, "dock_pockets")


def handle_compare_docking(state: WorkflowState) -> WorkflowState:
    return run_selected_tools(state, "compare_docking")


def handle_validate_pose(state: WorkflowState) -> WorkflowState:
    """Поза берётся из результата докинга этого же запроса, а не из демо-каталога."""
    if not state.observations.get("pose_file"):
        best = state.observations.get("best_pose")
        if best:
            state.observations["pose_file"] = best
            state.request.pose_file = best
    return run_selected_tools(state, "validate_pose")
