"""Валидация позы."""
from __future__ import annotations

from pathlib import Path

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.pose_validation import pass_rate, validate_pose_file

_REPO_ROOT = Path(__file__).resolve().parent.parent


def validate_pose(state: WorkflowState) -> WorkflowState:
    with state.agent_step("validation_agent", "pose_checks"):
        pose_file = state.observations.get("pose_file") or state.request.pose_file
        if not pose_file:
            state.add_error("Укажите pose_file для валидации.")
            return state
        pose_path = Path(pose_file)
        if not pose_path.is_absolute():
            pose_path = _REPO_ROOT / pose_path
        pdb_id = state.observations.get("pdb_id") or state.request.pdb_id
        protein_path = None
        if pdb_id:
            for folder in ("data/cache/pdb", "data/external/rcsb"):
                protein = _REPO_ROOT / folder / f"{str(pdb_id).upper()}.pdb"
                if protein.is_file():
                    protein_path = protein
                    break

        results = validate_pose_file(pose_path, protein_path)
        pct = pass_rate(results)
        state.observations["validation_checks"] = [r.__dict__ for r in results]
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=pct,
                unit="%",
                method="checks_passed_percent",
                source="pose_validation",
                context={"checks": state.observations["validation_checks"]},
            )
        )
        for r in results:
            if not r.passed:
                state.limitations.append(f"Не пройдена проверка {r.name}: {r.detail}")
        meta = state.observations.get("pdb_meta") or {}
        if meta.get("ligand_rmsd_reference_angstrom"):
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=meta["ligand_rmsd_reference_angstrom"],
                    unit="Å",
                    method="rmsd_to_crystal",
                    source="reference",
                    context={"note": "Эталон из кеша PDB, если есть кокристалл."},
                )
            )
    return state
