"""Идентификаторы сценариев, обязательные поля и демо-запросы."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

SCENARIO_TARGET_FOR_LIGAND = "target_for_ligand"
SCENARIO_DOCK_POCKETS = "dock_pockets"
SCENARIO_COMPARE_DOCKING = "compare_docking"
SCENARIO_VALIDATE_POSE = "validate_pose"

SCENARIO_IDS = (
    SCENARIO_TARGET_FOR_LIGAND,
    SCENARIO_DOCK_POCKETS,
    SCENARIO_COMPARE_DOCKING,
    SCENARIO_VALIDATE_POSE,
)

ScenarioId = str


@dataclass(frozen=True)
class ScenarioSpec:
    id: ScenarioId
    title: str
    description: str
    required_observation_keys: tuple[str, ...]
    demo_request: str | None = None


SCENARIOS: dict[ScenarioId, ScenarioSpec] = {
    SCENARIO_TARGET_FOR_LIGAND: ScenarioSpec(
        id=SCENARIO_TARGET_FOR_LIGAND,
        title="Поиск и обоснование мишени",
        description="Подбор белка и PDB по лиганду с опорой на BindingDB и литературу.",
        required_observation_keys=("ligand_query",),
        demo_request=(
            "Подбери белок для докинга ацетазоламида. "
            "Обоснуй выбор экспериментальными данными и предложи структуру из PDB."
        ),
    ),
    SCENARIO_DOCK_POCKETS: ScenarioSpec(
        id=SCENARIO_DOCK_POCKETS,
        title="Поиск карманов и докинг",
        description="P2Rank + Vina; на демо — загрузка сохранённого прогона.",
        required_observation_keys=("pdb_id", "ligand_id"),
        demo_request="Выполни докинг ацетазоламида в структуру 3HS4 (два верхних кармана P2Rank).",
    ),
    SCENARIO_COMPARE_DOCKING: ScenarioSpec(
        id=SCENARIO_COMPARE_DOCKING,
        title="Анализ результатов докинга",
        description="Сравнение лигандов по карманам на таблице poses.csv.",
        required_observation_keys=(),
        demo_request=(
            "Сравни три лиганда в двух карманах. Покажи, меняется ли их порядок "
            "по оценке Vina при смене кармана."
        ),
    ),
    SCENARIO_VALIDATE_POSE: ScenarioSpec(
        id=SCENARIO_VALIDATE_POSE,
        title="Проверка правдоподобности позы",
        description="Геометрические проверки и RMSD к кристаллу при наличии.",
        required_observation_keys=("pdb_id", "ligand_id", "pose_file"),
        demo_request="Проверь правдоподобность лучшей позы ацетазоламида в 3HS4.",
    ),
}

# Ключевые слова для резервной маршрутизации без LLM (порядок важен при пересечениях).
_KEYWORD_RULES: list[tuple[tuple[str, ...], ScenarioId]] = [
    (("правдоподоб", "posebuster", "валидац", "проверь поз"), SCENARIO_VALIDATE_POSE),
    (("сравни", "теплов", "два карман", "три лиганд", "poses.csv"), SCENARIO_COMPARE_DOCKING),
    (("подбери белок", "обоснуй выбор", "bindingdb", "для докинга"), SCENARIO_TARGET_FOR_LIGAND),
    (("выполни докинг", "p2rank", "запусти vina"), SCENARIO_DOCK_POCKETS),
]


def classify_by_keywords(text: str) -> ScenarioId | None:
    found = classify_all_keywords(text)
    return found[0] if found else None


def classify_all_keywords(text: str) -> list[ScenarioId]:
    """Все подошедшие сценарии в научном порядке, не в порядке ключевых слов."""
    lower = text.lower()
    matched = {
        scenario_id
        for keywords, scenario_id in _KEYWORD_RULES
        if any(k in lower for k in keywords)
    }
    if SCENARIO_COMPARE_DOCKING in matched and "выполни докинг" not in lower and "p2rank" not in lower:
        matched.discard(SCENARIO_DOCK_POCKETS)
    if SCENARIO_VALIDATE_POSE in matched and "подбери белок" not in lower:
        matched.discard(SCENARIO_TARGET_FOR_LIGAND)
    return [sid for sid in SCENARIO_IDS if sid in matched]


def missing_required_keys(observations: dict, scenario_id: ScenarioId) -> list[str]:
    spec = SCENARIOS.get(scenario_id)
    if spec is None:
        return ["unknown_scenario"]
    missing = []
    for key in spec.required_observation_keys:
        val = observations.get(key)
        if val is None or (isinstance(val, str) and not val.strip()):
            missing.append(key)
    return missing
