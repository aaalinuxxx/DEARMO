"""Агент химии: RDKit дескрипторы."""
from __future__ import annotations

from pathlib import Path

from models.docking import load_poses_csv
from models.result import Evidence
from orchestration.state import WorkflowState
from tools.pubchem_client import resolve_compound
from tools.rdkit_tools import calculate_descriptors


def resolve_ligand(state: WorkflowState) -> WorkflowState:
    with state.agent_step("chemistry_agent", "resolve_ligand"):
        ligand = state.observations.get("ligand_query") or state.observations.get("ligand_id")
        if not ligand and state.request.molecules:
            ligand = state.request.molecules[0]
        if not ligand:
            return state
        state.observations["ligand_query"] = ligand
        info = resolve_compound(ligand)
        smiles = info.get("connectivity_smiles") or state.observations.get("smiles")
        if smiles:
            state.observations["smiles"] = smiles
            desc = calculate_descriptors(smiles)
            if desc.status == "success":
                for item in desc.evidence:
                    state.add_evidence(item)
                state.limitations.extend(desc.limitations)
    return state


def append_ligand_descriptors_for_poses(state: WorkflowState, csv_path: Path) -> WorkflowState:
    with state.agent_step("chemistry_agent", "ligand_descriptors"):
        rows = load_poses_csv(csv_path)
        ligands = sorted({r.ligand_id for r in rows})
        smiles_map = {
            "acetazolamide": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
            "benzamidine": "NC(=N)C1=CC=CC=C1",
            "methylacetate": "CC(=O)OC",
        }
        for lig in ligands:
            smi = smiles_map.get(lig)
            if not smi:
                continue
            desc = calculate_descriptors(smi)
            for item in desc.evidence:
                if item.method in ("molecular_weight", "logp"):
                    tagged = Evidence(
                        kind=item.kind,
                        value=item.value,
                        unit=item.unit,
                        source=item.source,
                        method=f"{item.method}_{lig}",
                        context={**item.context, "ligand_id": lig},
                    )
                    state.add_evidence(tagged)
        state.limitations.append("Молекулярная масса и logP относятся к лиганду и не меняются при смене кармана.")
    return state


def analyze(state: WorkflowState) -> WorkflowState:
    return resolve_ligand(state)
