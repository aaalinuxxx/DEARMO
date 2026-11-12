# Сценарии DEARMO (кейс СИНТЕЗ №1)

Файл сгенерирован командой `python cli.py demo` 2026-11-12 15:26 UTC. Колонка «фактический результат» — вывод этого прогона, а не описание намерений.

| Сценарий | Статус | Шаги | Чем интересен |
| --- | --- | --- | --- |
| Поиск и обоснование мишени | `success` | ok:target_for_ligand | Живые источники и решение LLM, подкреплённое данными, а не памятью модели. |
| Поиск карманов и докинг | `success` | ok:dock_pockets | Реальный расчёт на сцене за секунды для пары, которой нет в репозитории. |
| Анализ результатов докинга | `success` | ok:dock_pockets, ok:compare_docking | Разрешение зависимостей: таблицы поз нет, система сама добавляет докинг. |
| Проверка правдоподобности позы | `success` | ok:dock_pockets, ok:validate_pose | Цепочка из трёх агентов по запросу об одном шаге; чужая поза не подставляется. |

## Поиск и обоснование мишени (`target_for_ligand`)

**Запрос:** Подбери белок для докинга иматиниба. Обоснуй выбор экспериментальными данными и предложи структуру из PDB.

**Ожидаемый ответ:** CID и структура из PubChem, дескрипторы RDKit (calculation), измерения BindingDB по этому соединению (database), выбранная мишень и PDB с обоснованием модели (hypothesis).

**Действия системы:** planner → tool_router → chemistry_agent → data_agent → structure_agent → data_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `molecular_weight` = 493.615 g/mol [calculation, RDKit]
- `logp` = 4.59 [calculation, RDKit]
- `tpsa` = 86.28 Å² [calculation, RDKit]
- `hbd` = 2 [calculation, RDKit]
- `hba` = 7 [calculation, RDKit]
- `rotatable_bonds` = 7 [calculation, RDKit]

Артефакты: `evidence_csv`, `report_json`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Контекст RAG — выдержки из заранее загруженных источников, не полнотекстовый поиск
- Критик отметил неполноту ответа: отсутствует обоснование выбора белка экспериментальными данными. Автоматически это не достраивается.

## Поиск карманов и докинг (`dock_pockets`)

**Запрос:** Выполни докинг эрлотиниба в структуру 1M17.

**Ожидаемый ответ:** Оценка Vina в ккал/моль, посчитанная в этом запуске (calculation), файл позы, poses.csv и report.md в outputs/.

**Действия системы:** planner → tool_router → chemistry_agent → structure_agent → docking_agent → orc → analysis_agent → writer → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `molecular_weight` = 393.443 g/mol [calculation, RDKit]
- `logp` = 3.405 [calculation, RDKit]
- `tpsa` = 74.73 Å² [calculation, RDKit]
- `hbd` = 1 [calculation, RDKit]
- `hba` = 7 [calculation, RDKit]
- `rotatable_bonds` = 10 [calculation, RDKit]

Артефакты: `evidence_csv`, `heatmap_png`, `pose_pdbqt`, `poses_csv`, `receptor_pdbqt`, `report_json`, `report_md`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Рецептор переведён в PDBQT из координат ATOM без отдельной оптимизации зарядов. Карман взят по лиганду кристалла, если он есть, иначе по центру белка.
- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- Карманы обработаны параллельными воркерами ORC по позам, которые в этом запуске посчитал AutoDock Vina.
- Критик отметил неполноту ответа: dock_pockets. Автоматически это не достраивается.

## Анализ результатов докинга (`compare_docking`)

**Запрос:** Сравни лиганды по карманам для эрлотиниба в 1M17 и покажи, меняется ли их порядок по оценке Vina.

**Ожидаемый ответ:** Таблица poses.csv, тепловая карта PNG и признак смены порядка. Масса и logP от кармана не зависят.

**Действия системы:** planner → orchestrator → tool_router → chemistry_agent → structure_agent → docking_agent → orc → analysis_agent → writer → orchestrator → tool_router → analysis_agent → chemistry_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Шаги, добавленные системой самостоятельно:** dock_pockets (не хватало poses_csv_path)

**Фактический результат:** `success`

Ключевые значения:

- `molecular_weight` = 393.443 g/mol [calculation, RDKit]
- `logp` = 3.405 [calculation, RDKit]
- `tpsa` = 74.73 Å² [calculation, RDKit]
- `hbd` = 1 [calculation, RDKit]
- `hba` = 7 [calculation, RDKit]
- `rotatable_bonds` = 10 [calculation, RDKit]

Артефакты: `evidence_csv`, `heatmap_png`, `pose_pdbqt`, `poses_csv`, `receptor_pdbqt`, `report_json`, `report_md`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Рецептор переведён в PDBQT из координат ATOM без отдельной оптимизации зарядов. Карман взят по лиганду кристалла, если он есть, иначе по центру белка.
- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- Карманы обработаны параллельными воркерами ORC по позам, которые в этом запуске посчитал AutoDock Vina.
- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- Молекулярная масса, logP и QED относятся к лиганду и не меняются при смене кармана.
- Критик отметил неполноту ответа: Порядок лигандов по оценке Vina в карманах 1M17 не указан. Автоматически это не достраивается.

## Проверка правдоподобности позы (`validate_pose`)

**Запрос:** Проверь правдоподобность позы эрлотиниба в 1M17.

**Ожидаемый ответ:** Доля пройденных геометрических проверок и их список. Поза берётся только из расчёта этого же запроса.

**Действия системы:** planner → orchestrator → tool_router → chemistry_agent → structure_agent → docking_agent → orc → analysis_agent → writer → orchestrator → tool_router → structure_agent → validation_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Шаги, добавленные системой самостоятельно:** dock_pockets (не хватало pose_file)

**Фактический результат:** `success`

Ключевые значения:

- `molecular_weight` = 393.443 g/mol [calculation, RDKit]
- `logp` = 3.405 [calculation, RDKit]
- `tpsa` = 74.73 Å² [calculation, RDKit]
- `hbd` = 1 [calculation, RDKit]
- `hba` = 7 [calculation, RDKit]
- `rotatable_bonds` = 10 [calculation, RDKit]

Артефакты: `evidence_csv`, `heatmap_png`, `pose_pdbqt`, `poses_csv`, `receptor_pdbqt`, `report_json`, `report_md`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Рецептор переведён в PDBQT из координат ATOM без отдельной оптимизации зарядов. Карман взят по лиганду кристалла, если он есть, иначе по центру белка.
- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- Карманы обработаны параллельными воркерами ORC по позам, которые в этом запуске посчитал AutoDock Vina.
