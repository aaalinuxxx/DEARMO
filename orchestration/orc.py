"""ORC: кластеры-карманы, параллельные воркеры, writer и summary."""
from __future__ import annotations

import asyncio
import os
from pathlib import Path

from agents import analysis_agent, chemistry_agent, docking_agent, structure_agent
from agents.pocket_worker import run_pocket
from agents.writer_agent import write_report
from models.docking import PoseRow, load_poses_csv
from orchestration.pipeline import CLUSTER, INGEST, PARALLEL_SEARCH, PROPERTY_SUMMARY
from orchestration.state import WorkflowState

_ROOT = Path(__file__).resolve().parent.parent


def _max_pockets() -> int:
    return int(os.environ.get("ORC_MAX_POCKETS", "2"))


def _cluster_pockets(rows: list[PoseRow]) -> list[str]:
    ordered: list[str] = []
    for row in sorted(rows, key=lambda item: item.pocket_rank):
        if row.pocket_id not in ordered:
            ordered.append(row.pocket_id)
        if len(ordered) >= _max_pockets():
            break
    return ordered


async def _search(rows: list[PoseRow], pockets: list[str], ligand: str | None, terms: list[str]):
    async def one(pocket_id: str):
        return await asyncio.to_thread(
            run_pocket,
            pocket_id=pocket_id,
            rows=rows,
            ligand_id=ligand,
            query_terms=terms,
        )

    return await asyncio.gather(*(one(pocket_id) for pocket_id in pockets))


def run_dock_pipeline(state: WorkflowState, example_dir: Path | None = None) -> WorkflowState:
    state.observations["pipeline_phase"] = INGEST
    if not state.observations.get("pdb_id") and state.request.pdb_id:
        state.observations["pdb_id"] = state.request.pdb_id
    if not state.observations.get("ligand_id"):
        ligand = state.observations.get("ligand_query") or (
            state.request.ligand_ids[0] if state.request.ligand_ids else None
        )
        if ligand:
            state.observations["ligand_id"] = ligand
    state = chemistry_agent.resolve_ligand(state)
    state = structure_agent.select_structure(state)
    state = docking_agent.run_or_load(state, example_dir)
    if not state.artifacts.get("poses_csv"):
        return state

    csv_path = Path(state.artifacts["poses_csv"])
    rows = load_poses_csv(csv_path)
    state.observations["pipeline_phase"] = CLUSTER
    pockets = _cluster_pockets(rows)
    state.observations["pockets"] = pockets
    with state.agent_step("orc", f"cluster:{','.join(pockets)}"):
        pass

    ligand = state.observations.get("ligand_id")
    terms = [
        str(item)
        for item in (
            ligand,
            state.observations.get("target_hint"),
            state.observations.get("pdb_id"),
        )
        if item
    ]
    state.observations["pipeline_phase"] = PARALLEL_SEARCH
    items = asyncio.run(_search(rows, pockets, ligand, terms))
    state.work_queue.extend(items)
    for item in items:
        state.trace("orc", f"queue:{item.pocket_id}:{item.status}", detail=item.agent)
        state.evidence.extend(item.evidence)
        for error in item.errors:
            state.add_error(error)

    state.observations["pipeline_phase"] = PROPERTY_SUMMARY
    state.observations["poses_csv_path"] = str(csv_path)
    state = analysis_agent.compare_pockets(state, csv_path)
    state = write_report(state)
    if state.observations.get("docking_mode") == "live":
        state.limitations.append(
            "Карманы обработаны параллельными воркерами ORC по позам, которые в этом запуске посчитал AutoDock Vina."
        )
    else:
        state.limitations.append(
            "Карманы обработаны параллельными воркерами ORC по сохранённому прогону этой же пары."
        )
    return state
