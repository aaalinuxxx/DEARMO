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
            "Подбери белок-мишень для ацетазоламида и обоснуй выбор экспериментальными данными."
        ),
        expected=(
            "Идентичность лиганда: CID и канонический SMILES из PubChem. "
            "Экспериментальное подтверждение: тип измерения (Ki, Kd или IC50), значение, "
            "единицы и ссылка на запись BindingDB; если для точной пары данных нет — "
            "перечисление того, что нашлось. Мишень: название, организм, accession UniProt "
            "и обоснование выбора (hypothesis). Структура: PDB ID, цепь, разрешение и "
            "значимые кофакторы, например ион Zn."
        ),
        shows=(
            "Решение LLM, ограниченное данными: accession приходит из UniProt, структура "
            "ищется по нему уровнями A–D, а выдуманная мишень и ложное обоснование отбрасываются."
        ),
    ),
    SCENARIO_DOCK_POCKETS: ScenarioSpec(
        id=SCENARIO_DOCK_POCKETS,
        title="Поиск карманов и докинг",
        description=(
            "P2Rank предсказывает карманы по структуре без лиганда, AutoDock Vina считает "
            "позы в каждом из топ-3, поза №1 — глобальный минимум Vina."
        ),
        required_observation_keys=("pdb_id", "ligand_id"),
        demo_request=(
            "Найди три наиболее перспективных кармана, проведи докинг в каждом "
            "и покажи лучшую позу по Vina."
        ),
        expected=(
            "Для каждого кармана: ID, ранг, балл P2Rank и координаты центра. Для каждой позы: "
            "номер, оценка Vina в ккал/моль и файл координат. Одна поза №1 с ID кармана; "
            "RMSD l.b./u.b. внутри кармана относительно позы №1 того же кармана. "
            "poses.csv, pockets.csv и run_meta.json в outputs/docking/."
        ),
        shows=(
            "Докинг вслепую: центры боксов берёт P2Rank, координаты кристаллического лиганда "
            "в расчёт не подаются и служат только контролем (ранг кармана, RMSD до кристалла)."
        ),
    ),
    SCENARIO_COMPARE_DOCKING: ScenarioSpec(
        id=SCENARIO_COMPARE_DOCKING,
        title="Сравнение поз",
        description=(
            "Разница оценок Vina между позами, контакты с белком по каждой позе, "
            "расстояния до ключевых остатков и разбор закономерностей в этих рядах."
        ),
        required_observation_keys=("poses_csv_path",),
        demo_request=(
            "Просмотри таблицу универсальных свойств белка, которую мы построили "
            "в докинге, и покажи, какие взаимодействия с белком меняются; найди "
            "интересные случаи, где свойство меняется неожиданно или есть "
            "закономерность в падениях, подъёмах и всплесках, и составь аналитику."
        ),
        expected=(
            "Разница оценок Vina относительно позы №1 в ккал/моль и график по позам. "
            "Карта «поза × остаток»: где появляются и исчезают водородные связи, "
            "солевые мостики и другие контакты. График расстояний в Å между одними "
            "и теми же парами атомов во всех позах с подписями атомов. Таблицы "
            "protein_state.csv и pose_contacts.csv и вывод, какую позу проверить "
            "подробнее и по какой причине."
        ),
        shows=(
            "Числа считает код по координатам, а модель только выбирает позу из "
            "готового списка: у каждого расстояния на графике есть оба атома и "
            "файл позы, поэтому его можно сверить вручную."
        ),
    ),
    SCENARIO_VALIDATE_POSE: ScenarioSpec(
        id=SCENARIO_VALIDATE_POSE,
        title="Проверка геометрии позы №1",
        description=(
            "Столкновения с белком и внутри лиганда, длины связей и валентные углы, "
            "стереохимия, планарность колец и координация металла, если он есть."
        ),
        required_observation_keys=("pdb_id", "ligand_id", "pose_file"),
        demo_request="Проверь геометрию позы №1 и покажи, какие критерии она прошла.",
        expected=(
            "Вердикт по позе, число пройденных проверок из числа применимых и процент "
            "пройденности. По каждой проверке — измеренное значение и порог: пара атомов, "
            "остаток белка и расстояние для столкновения; значение и допустимый диапазон "
            "для связи или угла; отклонение от кристалла для контакта с металлом. "
            "Критические нарушения перечислены отдельно. Поза берётся только из расчёта "
            "этого же запроса."
        ),
        shows=(
            "Вердикт считает код по измерениям, а не LLM по тексту: критическое "
            "столкновение отклоняет позу даже при четырёх пройденных проверках из пяти. "
            "Проверки, которые не к чему применить, не занижают процент."
        ),
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
    (
        (
            "правдоподоб", "posebuster", "валидац", "проверь поз", "проверь геометри",
            "геометри", "столкновен", "стереохими", "планарн", "какие критерии",
        ),
        SCENARIO_VALIDATE_POSE,
    ),
    (
        (
            "сравни",
            "теплов",
            "два карман",
            "три лиганд",
            "poses.csv",
            "сравнение поз",
            "взаимодействия с белком",
            "какие взаимодействия",
            "контакты с белком",
            "свойств белка",
            "закономерност",
            "всплеск",
            "аналитик",
        ),
        SCENARIO_COMPARE_DOCKING,
    ),
    (("подбери белок", "обоснуй выбор", "bindingdb", "для докинга"), SCENARIO_TARGET_FOR_LIGAND),
    (("выполни докинг", "p2rank", "запусти vina", "перспективных карман", "лучшую позу"), SCENARIO_DOCK_POCKETS),
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
