"""Обработчики сценариев — детерминированные цепочки агентов."""
from __future__ import annotations

from pathlib import Path

from agents import analysis_agent, chemistry_agent, data_agent, docking_agent, structure_agent, validation_agent
from orchestration.state import WorkflowState

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_POSES = _REPO_ROOT / "data" / "examples" / "acetazolamide" / "poses.csv"
_DEFAULT_EXAMPLE_DIR = _REPO_ROOT / "data" / "examples" / "acetazolamide"


def handle_target_for_ligand(state: WorkflowState) -> WorkflowState:
    state = chemistry_agent.resolve_ligand(state)
    state = data_agent.collect_target_data(state)
    state = structure_agent.select_structure(state)
    state = data_agent.attach_literature_rag(state)
    return state


def handle_dock_pockets(state: WorkflowState) -> WorkflowState:
    if not state.observations.get("pdb_id") and state.request.pdb_id:
        state.observations["pdb_id"] = state.request.pdb_id
    if not state.observations.get("ligand_id") and state.request.ligand_ids:
        state.observations["ligand_id"] = state.request.ligand_ids[0]
    state = chemistry_agent.resolve_ligand(state)
    state = structure_agent.select_structure(state)
    state = docking_agent.run_or_load(state, example_dir=_DEFAULT_EXAMPLE_DIR)
    return state


def handle_compare_docking(state: WorkflowState) -> WorkflowState:
    csv_path = state.request.poses_csv_path or state.observations.get("poses_csv_path")
    if not csv_path:
        csv_path = str(_DEFAULT_POSES)
        state.observations["poses_csv_path"] = csv_path
    state.artifacts["poses_csv"] = csv_path
    state = analysis_agent.compare_pockets(state, Path(csv_path))
    state = chemistry_agent.append_ligand_descriptors_for_poses(state, Path(csv_path))
    return state


def handle_validate_pose(state: WorkflowState) -> WorkflowState:
    if not state.observations.get("pose_file"):
        best = _DEFAULT_EXAMPLE_DIR / "poses" / "acetazolamide_pocket1_pose1.pdbqt"
        if best.is_file():
            state.observations["pose_file"] = str(best)
            state.request.pose_file = str(best)
    if not state.observations.get("pdb_id"):
        state.observations["pdb_id"] = state.request.pdb_id or "3HS4"
    if not state.observations.get("ligand_id"):
        state.observations["ligand_id"] = "acetazolamide"
    state = structure_agent.select_structure(state)
    state = validation_agent.validate_pose(state)
    return state
