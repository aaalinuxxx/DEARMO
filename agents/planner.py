"""Планировщик: сценарий, сущности, проверка неоднозначности и неподтверждённых выводов (GigaChat + keyword fallback)."""
from __future__ import annotations

import re
from pathlib import Path

from orchestration.gigachat_client import PlannerExtraction, plan_with_gigachat
from orchestration.scenarios import SCENARIO_IDS, SCENARIO_TARGET_FOR_LIGAND, classify_all_keywords
from orchestration.state import WorkflowState
from orchestration.task_graph import order_scenarios
from tools.pubchem_client import resolve_compound

_PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "planner.txt"

_PDB_RE = re.compile(r"\b([0-9][A-Za-z0-9]{3})\b")

_KNOWN_LIGANDS: tuple[tuple[tuple[str, ...], str], ...] = (
    (("ацетазоламид", "acetazolamide"), "acetazolamide"),
    (("метазоламид", "methazolamide"), "methazolamide"),
    (("этоксзоламид", "ethoxzolamide"), "ethoxzolamide"),
    (("дорзоламид", "dorzolamide"), "dorzolamide"),
    (("бринзоламид", "brinzolamide"), "brinzolamide"),
    (("бензамидин", "benzamidine"), "benzamidine"),
    (("аспирин", "aspirin"), "aspirin"),
)

_AMBIGUOUS_LIGAND_TERMS = ("сульфонамид", "sulfonamide", "любой ингибитор", "неизвестный лиганд")
_AMBIGUOUS_TARGET_TERMS = ("любая карбоангидраза", "какой-нибудь белок", "неизвестная мишень")

_UNSUPPORTED_CLAIM_KEYWORDS = (
    "клиническ",
    "токсичност",
    "in vivo",
    "побочн",
    "гарантированно вылеч",
    "безопасен для человека",
    "пройдет клинические испытания",
    "clinical efficacy",
)


def _load_system_prompt() -> str:
    if _PROMPT_PATH.is_file():
        return _PROMPT_PATH.read_text(encoding="utf-8")
    return "Ты планировщик drug discovery. Извлеки сценарий и сущности из запроса."


def _keyword_extraction(text: str) -> PlannerExtraction:
    ids = classify_all_keywords(text)
    lower = text.lower()
    if not ids and any(k in lower for k in ("аналог", "похож", "активност", "свойств")):
        ids = [SCENARIO_TARGET_FOR_LIGAND]
    scenario_id = ids[0] if ids else None
    ligands: list[str] = []

    for aliases, canonical in _KNOWN_LIGANDS:
        if any(a in lower for a in aliases):
            ligands.append(canonical)

    if not ligands:
        for amb in _AMBIGUOUS_LIGAND_TERMS:
            if amb in lower:
                ligands.append(amb)
                break

    pdb_match = _PDB_RE.search(text)
    pdb_id = pdb_match.group(1).upper() if pdb_match else None
    ligand_id = ligands[0] if ligands else None
    return PlannerExtraction(
        scenario_id=scenario_id,
        ligands=ligands,
        target_hint="carbonic anhydrase" if any(x in lower for x in ("ацетазоламид", "acetazolamide", "карбоангидраз")) else None,
        pdb_id=pdb_id,
        poses_csv_path=None,
        pose_file=None,
        ligand_id=ligand_id,
        explanation_hint=text[:200],
        source="keyword",
    )


def _check_ambiguity_and_unsupported_claims(state: WorkflowState) -> None:
    """Проверить запрос на неоднозначность мишени/соединения и неподтверждённые данными выводы."""
    lower = state.request.text.lower()

    if not state.request.pdb_id and not state.request.target_id:
        for term in _AMBIGUOUS_TARGET_TERMS:
            if term in lower:
                state.add_error(
                    f"Уточните мишень: формулировка '{term}' неоднозначна (укажите конкретную изоформу, например CA II / CHEMBL205, или PDB ID)."
                )

    # Проверяем, не задал ли пользователь уже конкретную (не неоднозначную) молекулу
    active_mol = (state.observations.get("ligand_query") or "").strip().lower()
    has_concrete_mol = bool(active_mol) and not any(amb in active_mol for amb in _AMBIGUOUS_LIGAND_TERMS)

    if not has_concrete_mol:
        for term in _AMBIGUOUS_LIGAND_TERMS:
            if term in lower or term in active_mol:
                pub_info = resolve_compound("сульфонамид" if "сульфонамид" in term or "sulfonamide" in term else "ингибитор")
                if pub_info.get("candidates"):
                    state.observations["clarification_candidates"] = pub_info["candidates"]
                state.add_error(
                    f"Уточните соединение: '{term}' — это класс веществ, а не конкретная молекула. Укажите название или SMILES."
                )
                break

    if any(kw in lower for kw in _UNSUPPORTED_CLAIM_KEYWORDS):
        note = (
            "Ограничение интерпретации: имеющиеся in vitro данные (Ki/IC50), расчётные дескрипторы RDKit "
            "и оценки докинга Vina НЕ позволяют делать выводы о клинической эффективности, токсичности in vivo "
            "или безопасности для человека."
        )
        if note not in state.limitations:
            state.limitations.append(note)
        state.observations["unsupported_conclusion_requested"] = True


def _apply_extraction(state: WorkflowState, ext: PlannerExtraction) -> WorkflowState:
    keyword_ids = classify_all_keywords(state.request.text)
    chosen = []
    if ext.scenario_id and ext.scenario_id in SCENARIO_IDS:
        chosen.append(ext.scenario_id)
    for sid in keyword_ids:
        if sid not in chosen:
            chosen.append(sid)
    chosen = order_scenarios(chosen)
    for sid in keyword_ids:
        if sid not in chosen:
            chosen.append(sid)
    # --- ВСТАВИТЬ ЭТОТ БЛОК ---
    lower_text = state.request.text.lower()
    if any(kw in lower_text for kw in ("подбери белок", "обоснуй выбор", "для докинга")):
        if "target_for_ligand" not in chosen:
            chosen.insert(0, "target_for_ligand")
    # --------------------------
    chosen = order_scenarios(chosen)
    state.observations["scenario_ids"] = chosen
    state.scenario_id = chosen[0] if chosen else None
    state.request.scenario = state.scenario_id
    state.planner_source = ext.source

    if not state.request.molecules and ext.ligands:
        state.request.molecules = ext.ligands
    if not state.request.ligand_ids and ext.ligands:
        state.request.ligand_ids = ext.ligands

    if state.request.molecules:
        state.observations["ligand_query"] = state.request.molecules[0]
    elif ext.ligands:
        state.observations["ligand_query"] = ext.ligands[0]

    if state.request.ligand_ids:
        state.observations["ligand_id"] = state.request.ligand_ids[0]
    elif ext.ligand_id:
        state.observations["ligand_id"] = ext.ligand_id

    if ext.target_hint:
        state.observations["target_hint"] = ext.target_hint
    if ext.pdb_id and not state.request.pdb_id:
        state.request.pdb_id = ext.pdb_id
    if state.request.pdb_id:
        state.observations["pdb_id"] = state.request.pdb_id
    if ext.poses_csv_path:
        state.request.poses_csv_path = ext.poses_csv_path
        state.observations["poses_csv_path"] = ext.poses_csv_path
    if ext.pose_file:
        state.request.pose_file = ext.pose_file
        state.observations["pose_file"] = ext.pose_file

    _check_ambiguity_and_unsupported_claims(state)
    state.plan = [f"{i + 1}:{sid}" for i, sid in enumerate(chosen)] + [f"planner={ext.source}"]
    return state


def plan(state: WorkflowState) -> WorkflowState:
    with state.agent_step("planner", "classify"):
        text = state.request.text
        ext = plan_with_gigachat(text, _load_system_prompt())
        if ext is None:
            ext = _keyword_extraction(text)
        return _apply_extraction(state, ext)
