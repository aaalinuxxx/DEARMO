"""Планировщик: сценарий и сущности. Qwen через Ollama, иначе ключевые слова."""
from __future__ import annotations

import re
from pathlib import Path

from orchestration.llm_client import PlannerExtraction, plan_with_llm
from orchestration.scenarios import SCENARIO_IDS, SCENARIO_TARGET_FOR_LIGAND, classify_all_keywords
from orchestration.state import WorkflowState
from orchestration.task_graph import order_scenarios
from tools.pubchem_client import resolve_compound
from tools.translit import looks_russian, resolve_russian_name

_PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "planner.txt"

_PDB_RE = re.compile(r"\b([0-9][A-Za-z0-9]{3})\b")

_AMBIGUOUS_LIGAND_TERMS = ("сульфонамид", "sulfonamide", "любой ингибитор", "неизвестный лиганд")
_AMBIGUOUS_TARGET_TERMS = ("какой-нибудь белок", "неизвестная мишень", "любой белок")

# Латинское слово или SMILES в кавычках — всё, что резервный разбор может узнать о молекуле
# без словаря. Нормализацию русских и торговых названий делает LLM-планировщик.
_QUOTED_RE = re.compile(r"[«\"\']([^«»\"\']{2,60})[»\"\']")
_LATIN_TOKEN_RE = re.compile(r"\b([A-Za-z][A-Za-z0-9\-]{3,40})\b")
_SMILES_RE = re.compile(r"\b([A-Za-z0-9@+\-\[\]()=#\\/%]{6,120})\b")
_STOPWORDS = {
    "pdb", "smiles", "vina", "docking", "protein", "ligand", "chembl", "bindingdb",
    "pubchem", "rdkit", "csv", "json", "posebusters", "qed", "logp", "tpsa",
}

_RUSSIAN_TOKEN_RE = re.compile(r"\b([а-яё]{6,40})\b", re.I)
# Частые слова запроса, которые заведомо не являются названием вещества.
# Список лингвистический, а не химический: он не содержит ни одной молекулы.
_RU_STOPWORDS = {
    "подбери", "выполни", "проверь", "сравни", "структуру", "структура", "структуре",
    "правдоподобность", "обоснуй", "предложи", "покажи", "рассчитай", "посчитай",
    "докинга", "докинг", "докинге", "белок", "белка", "мишень", "мишени", "лиганд",
    "лиганда", "лиганды", "лигандов", "карман", "кармане", "карманах", "карманов",
    "молекулы", "молекула", "молекулу", "соединения", "соединение", "соединений",
    "свойства", "свойств", "активность", "активности", "экспериментальными",
    "данными", "литературу", "аналоги", "аналогов", "порядок", "оценке", "оценка",
    "позы", "позу", "поза", "лучшей", "таблице", "таблицу", "результаты",
}

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


def _ligands_without_dictionary(text: str) -> list[str]:
    """Достать молекулу из текста без списка известных соединений.

    Берём то, что в кавычках, иначе — латинские токены, не являющиеся служебными
    словами. Русские названия здесь намеренно не распознаются: их переводит LLM,
    а резервный режим честно попросит уточнить.
    """
    found: list[str] = []
    for match in _QUOTED_RE.finditer(text):
        token = match.group(1).strip()
        if token and token.lower() not in _STOPWORDS:
            found.append(token)
    if found:
        return found[:3]
    for match in _LATIN_TOKEN_RE.finditer(text):
        token = match.group(1)
        if token.lower() in _STOPWORDS or _PDB_RE.fullmatch(token):
            continue
        found.append(token)
    if found:
        return found[:3]

    # Латиницы нет: пробуем русские слова через транслитерацию. Название считается
    # распознанным, только если PubChem точно нашёл по нему соединение, поэтому
    # случайное слово из запроса молекулой не станет.
    for token in _RUSSIAN_TOKEN_RE.findall(text):
        if len(token) < 6 or token.lower() in _RU_STOPWORDS:
            continue
        info = resolve_russian_name(token, resolve_compound)
        if info:
            return [info.get("title") or info["translit_variant"]]
    return []


def _keyword_extraction(text: str) -> PlannerExtraction:
    ids = classify_all_keywords(text)
    lower = text.lower()
    if not ids and any(k in lower for k in ("аналог", "похож", "активност", "свойств")):
        ids = [SCENARIO_TARGET_FOR_LIGAND]
    scenario_id = ids[0] if ids else None
    ligands = _ligands_without_dictionary(text)

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
        smiles=None,
        target_hint=None,
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
                    f"Уточните мишень: формулировка '{term}' неоднозначна. Назовите белок, его ChEMBL ID или PDB ID."
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

    if ext.smiles and not state.observations.get("smiles"):
        state.observations["smiles"] = ext.smiles
        state.observations["smiles_source"] = ext.source
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
        ext = plan_with_llm(text, _load_system_prompt())
        if ext is None:
            ext = _keyword_extraction(text)
        return _apply_extraction(state, ext)
