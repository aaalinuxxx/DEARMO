"""Планировщик: сценарий, сущности, ограничения (GigaChat + keyword fallback)."""
from __future__ import annotations

import re
from pathlib import Path

from orchestration.gigachat_client import PlannerExtraction, plan_with_gigachat
from orchestration.scenarios import SCENARIO_IDS, classify_all_keywords
from orchestration.task_graph import order_scenarios
from orchestration.state import WorkflowState

_PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "planner.txt"

_PDB_RE = re.compile(r"\b([0-9][A-Za-z0-9]{3})\b")


def _load_system_prompt() -> str:
    if _PROMPT_PATH.is_file():
        return _PROMPT_PATH.read_text(encoding="utf-8")
    return "Ты планировщик drug discovery. Извлеки сценарий и сущности из запроса."


def _keyword_extraction(text: str) -> PlannerExtraction:
    ids = classify_all_keywords(text)
    scenario_id = ids[0] if ids else None
    ligands: list[str] = []
    lower = text.lower()
    if "ацетазоламид" in lower or "acetazolamide" in lower:
        ligands.append("acetazolamide")
    pdb_match = _PDB_RE.search(text)
    pdb_id = pdb_match.group(1).upper() if pdb_match else None
    ligand_id = ligands[0] if ligands else None
    return PlannerExtraction(
        scenario_id=scenario_id,
        ligands=ligands,
        target_hint="carbonic anhydrase" if "ацетазоламид" in lower or "acetazolamide" in lower else None,
        pdb_id=pdb_id,
        poses_csv_path=None,
        pose_file=None,
        ligand_id=ligand_id,
        explanation_hint=text[:200],
        source="keyword",
    )


def _apply_extraction(state: WorkflowState, ext: PlannerExtraction) -> WorkflowState:
    keyword_ids = classify_all_keywords(state.request.text)
    chosen = []
    if ext.scenario_id and ext.scenario_id in SCENARIO_IDS:
        chosen.append(ext.scenario_id)
    for sid in keyword_ids:
        if sid not in chosen:
            chosen.append(sid)
    chosen = order_scenarios(chosen)
    state.observations["scenario_ids"] = chosen
    state.scenario_id = chosen[0] if chosen else None
    state.request.scenario = state.scenario_id
    state.planner_source = ext.source
    if ext.ligands:
        state.request.molecules = ext.ligands
        state.request.ligand_ids = ext.ligands
        state.observations["ligand_query"] = ext.ligands[0]
    if ext.ligand_id:
        state.observations["ligand_id"] = ext.ligand_id
    if ext.target_hint:
        state.observations["target_hint"] = ext.target_hint
    if ext.pdb_id:
        state.request.pdb_id = ext.pdb_id
        state.observations["pdb_id"] = ext.pdb_id
    if ext.poses_csv_path:
        state.request.poses_csv_path = ext.poses_csv_path
        state.observations["poses_csv_path"] = ext.poses_csv_path
    if ext.pose_file:
        state.request.pose_file = ext.pose_file
        state.observations["pose_file"] = ext.pose_file
    state.plan = [f"{i + 1}:{sid}" for i, sid in enumerate(chosen)] + [f"planner={ext.source}"]
    return state


def plan(state: WorkflowState) -> WorkflowState:
    with state.agent_step("planner", "classify"):
        text = state.request.text
        ext = plan_with_gigachat(text, _load_system_prompt())
        if ext is None:
            ext = _keyword_extraction(text)
        return _apply_extraction(state, ext)
