"""Выбор инструментов агентов: обязательный маршрут сценария и дополнение от локальной Qwen."""
from __future__ import annotations

from pathlib import Path

from agents import analysis_agent, chemistry_agent, data_agent, structure_agent, validation_agent
from orchestration.llm_client import ask_json, ollama_model
from orchestration.orc import run_dock_pipeline
from orchestration.state import WorkflowState

_PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "tool_router.txt"

_ORDER = (
    "rdkit_ligand",
    "activity_sources",
    "pdb_select",
    "literature",
    "dock_pipeline",
    "pocket_compare",
    "pose_descriptors",
    "pose_checks",
)

_PURPOSE = {
    "rdkit_ligand": "PubChem и расчётные дескрипторы RDKit лиганда",
    "activity_sources": "BindingDB по соединению, затем ChEMBL по мишени или по молекуле; аналоги, если их просят",
    "pdb_select": "Выбор PDB: указанный код, рекомендация по измерениям или поиск RCSB",
    "literature": "Локальные выдержки литературы по лиганду и мишени",
    "dock_pipeline": "Подготовка рецептора и лиганда, AutoDock Vina, отчёт по карманам",
    "pocket_compare": "Сравнение оценок Vina по карманам и тепловая карта",
    "pose_descriptors": "Дескрипторы RDKit для лигандов из таблицы поз",
    "pose_checks": "Геометрические проверки позы",
}

_PROPERTY_WORDS = ("logp", "qed", "дескриптор", "молекулярн", "липинск", "pains", "tpsa", "свойств")
_BIOLOGY_WORDS = ("белок", "мишен", "докинг", "pdb", "карман", "активност", "структур", "позу", "позы", "pose", "аффин")
_SKIP_LITERATURE = ("без литератур", "не ищи стать", "без статей")


def _load_prompt() -> str:
    if _PROMPT_PATH.is_file():
        return _PROMPT_PATH.read_text(encoding="utf-8")
    return "Выбери инструменты из списка и верни JSON с полем tools."


def property_only(text: str) -> bool:
    lower = text.lower()
    asks_property = any(word in lower for word in _PROPERTY_WORDS)
    asks_biology = any(word in lower for word in _BIOLOGY_WORDS)
    return asks_property and not asks_biology


def heuristic_tools(state: WorkflowState, scenario: str | None = None) -> list[str]:
    """Маршрут без модели: только те программы, которые отвечают на этот запрос."""
    scenario = scenario if scenario is not None else (state.scenario_id or "")
    text = state.request.text.lower()
    if scenario == "target_for_ligand" and property_only(text):
        tools = ["rdkit_ligand"]
    elif scenario == "target_for_ligand":
        tools = ["rdkit_ligand", "activity_sources", "pdb_select", "literature"]
    elif scenario == "dock_pockets":
        tools = ["dock_pipeline"]
    elif scenario == "compare_docking":
        tools = ["pocket_compare", "pose_descriptors"]
    elif scenario == "validate_pose":
        tools = ["pdb_select", "pose_checks"]
    else:
        tools = ["rdkit_ligand", "activity_sources"]
    if any(word in text for word in _SKIP_LITERATURE):
        tools = [name for name in tools if name != "literature"]
    return tools


def required_tools(state: WorkflowState, scenario: str | None = None) -> list[str]:
    """Шаги, без которых сценарий не отвечает на вопрос. Qwen не может их выкинуть."""
    scenario = scenario if scenario is not None else (state.scenario_id or "")
    if scenario == "target_for_ligand" and property_only(state.request.text):
        return ["rdkit_ligand"]
    if scenario == "target_for_ligand":
        required = ["rdkit_ligand", "activity_sources", "pdb_select"]
        if any(word in state.request.text.lower() for word in _SKIP_LITERATURE):
            return required
        return required
    return heuristic_tools(state, scenario)


def optional_tools(state: WorkflowState, scenario: str | None = None) -> set[str]:
    scenario = scenario if scenario is not None else (state.scenario_id or "")
    if scenario == "target_for_ligand":
        return {"literature"}
    return set()


def _order(names: list[str]) -> list[str]:
    rank = {name: index for index, name in enumerate(_ORDER)}
    seen: list[str] = []
    for name in sorted(names, key=lambda item: rank.get(item, 100)):
        if name not in seen and name in _PURPOSE:
            seen.append(name)
    return seen


def choose_tool_names(state: WorkflowState, scenario: str | None = None) -> tuple[list[str], str]:
    scenario = scenario if scenario is not None else state.scenario_id
    fallback = heuristic_tools(state, scenario)
    allowed = set(fallback) | set(required_tools(state, scenario)) | optional_tools(state, scenario)
    asked = _ask_qwen(state, [name for name in _ORDER if name in allowed], scenario)
    if not asked:
        return _order(fallback), "heuristic"
    chosen = [name for name in asked if name in allowed]
    for name in required_tools(state, scenario):
        if name not in chosen:
            chosen.append(name)
    if not chosen:
        chosen = list(fallback)
    return _order(chosen), f"ollama:{ollama_model()}"


def _ask_qwen(state: WorkflowState, allowed: list[str], scenario: str | None) -> list[str] | None:
    catalog = "\n".join(f"- {name}: {_PURPOSE[name]}" for name in allowed)
    known = ", ".join(
        f"{key}={state.observations.get(key)}"
        for key in ("ligand_query", "target_hint", "pdb_id", "smiles")
        if state.observations.get(key)
    )
    user = (
        f"Сценарий: {scenario}\n"
        f"Запрос: {state.request.text}\n"
        f"Уже известно: {known or 'ничего'}\n"
        f"Доступные инструменты:\n{catalog}"
    )
    data = ask_json(_load_prompt(), user)
    if not data:
        return None
    tools = data.get("tools")
    if isinstance(tools, str):
        tools = [tools]
    if not isinstance(tools, list):
        return None
    names = [str(item).strip() for item in tools if str(item).strip()]
    return names or None


def _poses_csv(state: WorkflowState) -> Path | None:
    """Таблица поз этого запроса: явный путь или результат докинга. Демо-файл не подставляется."""
    for candidate in (
        state.request.poses_csv_path,
        state.observations.get("poses_csv_path"),
        state.artifacts.get("poses_csv"),
    ):
        if candidate and Path(candidate).is_file():
            return Path(candidate)
    return None


def _run_tool(state: WorkflowState, name: str) -> None:
    if name == "rdkit_ligand":
        chemistry_agent.resolve_ligand(state)
    elif name == "activity_sources":
        data_agent.collect_target_data(state)
    elif name == "pdb_select":
        structure_agent.select_structure(state)
    elif name == "literature":
        data_agent.run_rag(state)
    elif name == "dock_pipeline":
        run_dock_pipeline(state)
    elif name == "pocket_compare":
        csv_path = _poses_csv(state)
        if csv_path is None:
            state.add_error(
                "Нет таблицы поз для сравнения. Укажите путь к poses.csv или сначала выполните докинг."
            )
            return
        state.observations["poses_csv_path"] = str(csv_path)
        state.artifacts["poses_csv"] = str(csv_path)
        analysis_agent.compare_pockets(state, csv_path)
    elif name == "pose_descriptors":
        csv_path = _poses_csv(state)
        if csv_path is not None:
            chemistry_agent.append_ligand_descriptors_for_poses(state, csv_path)
    elif name == "pose_checks":
        validation_agent.validate_pose(state)


def run_selected_tools(state: WorkflowState, scenario_id: str | None = None) -> WorkflowState:
    scenario = scenario_id or state.scenario_id
    names, source = choose_tool_names(state, scenario)
    plans = list(state.observations.get("tool_plans") or [])
    plans.append({"scenario": scenario, "tools": names, "source": source})
    state.observations["tool_plans"] = plans
    state.observations["tool_plan"] = names
    state.observations["tool_plan_source"] = source
    with state.agent_step("tool_router", f"{source}:{','.join(names)}"):
        pass
    for name in names:
        _run_tool(state, name)
        if any("уточн" in error.lower() for error in state.errors):
            break
    return state
