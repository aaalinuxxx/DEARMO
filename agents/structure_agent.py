"""Агент структур: PDB метаданные и рекомендации."""
from __future__ import annotations

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.pdb_client import get_structure_metadata, search_by_ligand_target


def select_structure(state: WorkflowState) -> WorkflowState:
    with state.agent_step("structure_agent", "pdb_select"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else "acetazolamide"
        )
        pdb_id = state.observations.get("pdb_id") or state.request.pdb_id
        rec = state.observations.get("recommended_target") or {}
        if not pdb_id and rec.get("pdb_id"):
            pdb_id = rec["pdb_id"]
        if not pdb_id:
            candidates = search_by_ligand_target(str(ligand), state.observations.get("target_hint"))
            if candidates:
                pdb_id = candidates[0].get("pdb_id")
                state.observations["pdb_candidates"] = candidates
        if not pdb_id:
            state.add_error("Не выбран PDB ID. Укажите структуру явно.")
            return state

        state.observations["pdb_id"] = pdb_id
        state.request.pdb_id = pdb_id
        meta = get_structure_metadata(pdb_id)
        state.observations["pdb_meta"] = meta
        state.add_evidence(
            Evidence(
                kind="database",
                value=pdb_id,
                method="pdb_structure",
                source="RCSB PDB",
                context=meta,
            )
        )
        if meta.get("resolution_angstrom"):
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=meta["resolution_angstrom"],
                    unit="Å",
                    method="resolution",
                    source="RCSB PDB",
                    context={"pdb_id": pdb_id},
                )
            )
    return state
