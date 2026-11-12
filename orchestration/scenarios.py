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
    # Что система должна показать. Это формулировка «ожидаемого ответа» из кейса;
    # фактический результат пишет `cli.py demo` в docs/scenarios.md.
    expected: str | None = None
    # Чем сценарий интересен на защите.
    shows: str | None = None


SCENARIOS: dict[ScenarioId, ScenarioSpec] = {
    SCENARIO_TARGET_FOR_LIGAND: ScenarioSpec(
        id=SCENARIO_TARGET_FOR_LIGAND,
        title="Поиск и обоснование мишени",
        description="Подбор белка и PDB по лиганду с опорой на BindingDB и литературу.",
        required_observation_keys=("ligand_query",),
        demo_request=(
            "Подбери белок для докинга иматиниба. "
            "Обоснуй выбор экспериментальными данными и предложи структуру из PDB."
        ),
        expected=(
            "CID и структура из PubChem, дескрипторы RDKit (calculation), измерения "
            "BindingDB по этому соединению (database), выбранная мишень и PDB с "
            "обоснованием модели (hypothesis)."
        ),
        shows="Живые источники и решение LLM, подкреплённое данными, а не памятью модели.",
    ),
    SCENARIO_DOCK_POCKETS: ScenarioSpec(
        id=SCENARIO_DOCK_POCKETS,
        title="Поиск карманов и докинг",
        description="Подготовка рецептора и лиганда, AutoDock Vina, отчёт по карманам.",
        required_observation_keys=("pdb_id", "ligand_id"),
        demo_request="Выполни докинг эрлотиниба в структуру 1M17.",
        expected=(
            "Оценка Vina в ккал/моль, посчитанная в этом запуске (calculation), файл позы, "
            "poses.csv и report.md в outputs/."
        ),
        shows="Реальный расчёт на сцене за секунды для пары, которой нет в репозитории.",
    ),
    SCENARIO_COMPARE_DOCKING: ScenarioSpec(
        id=SCENARIO_COMPARE_DOCKING,
        title="Анализ результатов докинга",
        description="Сравнение лигандов и карманов на таблице poses.csv.",
        required_observation_keys=("poses_csv_path",),
        demo_request=(
            "Сравни лиганды по карманам для эрлотиниба в 1M17 и покажи, "
            "меняется ли их порядок по оценке Vina."
        ),
        expected=(
            "Таблица poses.csv, тепловая карта PNG и признак смены порядка. "
            "Масса и logP от кармана не зависят."
        ),
        shows="Разрешение зависимостей: таблицы поз нет, система сама добавляет докинг.",
    ),
    SCENARIO_VALIDATE_POSE: ScenarioSpec(
        id=SCENARIO_VALIDATE_POSE,
        title="Проверка правдоподобности позы",
        description="Геометрические проверки и RMSD к кристаллу при наличии.",
        required_observation_keys=("pdb_id", "ligand_id", "pose_file"),
        demo_request="Проверь правдоподобность позы эрлотиниба в 1M17.",
        expected=(
            "Доля пройденных геометрических проверок и их список. Поза берётся только "
            "из расчёта этого же запроса."
        ),
        shows="Цепочка из трёх агентов по запросу об одном шаге; чужая поза не подставляется.",
    ),
}

# Какой сценарий производит недостающий вход. Благодаря этой карте запрос
# «проверь позу X» не падает с «укажите pose_file», а сам достраивается докингом.
# None означает, что значение может прийти только от пользователя.
OBSERVATION_PRODUCERS: dict[str, ScenarioId | None] = {
    "pdb_id": SCENARIO_TARGET_FOR_LIGAND,
    "pose_file": SCENARIO_DOCK_POCKETS,
    "poses_csv_path": SCENARIO_DOCK_POCKETS,
    "ligand_query": None,
    "ligand_id": None,
}


def producer_for(key: str) -> ScenarioId | None:
    return OBSERVATION_PRODUCERS.get(key)


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
