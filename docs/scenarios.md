# Сценарии DEARMO (кейс СИНТЕЗ №1)

Файл сгенерирован командой `python cli.py demo` 2026-11-12 21:55 UTC. Колонка «фактический результат» — вывод этого прогона, а не описание намерений.

| Сценарий | Статус | Шаги | Чем интересен |
| --- | --- | --- | --- |
| Поиск и обоснование мишени | `success` | ok:target_for_ligand | Решение LLM, ограниченное данными: accession приходит из UniProt, структура ищется по нему уровнями A–D, а выдуманная мишень и ложное обоснование отбрасываются. |
| Поиск карманов и докинг | `success` | ok:dock_pockets | Докинг вслепую: центры боксов берёт P2Rank, координаты кристаллического лиганда в расчёт не подаются и служат только контролем (ранг кармана, RMSD до кристалла). |
| Сравнение поз | `success` | ok:compare_docking | Числа считает код по координатам, а модель только выбирает позу из готового списка: у каждого расстояния на графике есть оба атома и файл позы, поэтому его можно сверить вручную. |
| Проверка геометрии позы №1 | `success` | ok:validate_pose | Вердикт считает код по измерениям, а не LLM по тексту: критическое столкновение отклоняет позу даже при четырёх пройденных проверках из пяти. Проверки, которые не к чему применить, не занижают процент. |

## Поиск и обоснование мишени (`target_for_ligand`)

**Запрос:** Подбери белок-мишень для ацетазоламида и обоснуй выбор экспериментальными данными.

**Ожидаемый ответ:** Идентичность лиганда: CID и канонический SMILES из PubChem. Экспериментальное подтверждение: тип измерения (Ki, Kd или IC50), значение, единицы и ссылка на запись BindingDB; если для точной пары данных нет — перечисление того, что нашлось. Мишень: название, организм, accession UniProt и обоснование выбора (hypothesis). Структура: PDB ID, цепь, разрешение и значимые кофакторы, например ион Zn.

**Действия системы:** planner → tool_router → chemistry_agent → data_agent → target_agent → data_agent → structure_agent → data_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `uniprot_accession` = P00918 [database, UniProt]
- `pdb_structure` = 3HS4 [database, RCSB PDB]
- `resolution` = 1.1 Å [database, RCSB PDB]
- `structure_chain` = A [database, RCSB PDB]
- `metal_ion` = ZN [database, RCSB PDB]
- `Ki` = 25.0 nM [database, BindingDB]
- `target_choice_rationale` = Carbonic anhydrase 2 [hypothesis, ollama:qwen2.5:7b-instruct]
- `molecular_weight` = 222.251 g/mol [calculation, RDKit]

Артефакты: `evidence_csv`, `report_json`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Подсказка планировщика «Carbonic anhydrase II» не принята: в запросе такой мишени нет. Мишень выбрана по данным измерений и структур.
- Измерения BindingDB по «Carbonic anhydrase 2» относятся к мутантным вариантам или фрагментам (A65S,N67Q,E69T,I91L,F131V,K170E,L204A), дикого типа в выдаче нет: значения аффинности не переносятся на природный белок без оговорки.
- Обоснование выбора «Carbonic anhydrase 2» от модели отклонено: в нём утверждалось сильнейшее связывание, но среди кандидатов оно у другой мишени. Выбор сохранён, обоснование заменено на фактическое.
- Контекст RAG — выдержки из заранее загруженных источников, не полнотекстовый поиск

## Поиск карманов и докинг (`dock_pockets`)

**Запрос:** Найди три наиболее перспективных кармана, проведи докинг в каждом и покажи лучшую позу по Vina.

**Ожидаемый ответ:** Для каждого кармана: ID, ранг, балл P2Rank и координаты центра. Для каждой позы: номер, оценка Vina в ккал/моль и файл координат. Одна поза №1 с ID кармана; RMSD l.b./u.b. внутри кармана относительно позы №1 того же кармана. poses.csv, pockets.csv и run_meta.json в outputs/docking/.

**Контекст из предыдущего шага цепочки:** `molecules` = ['acetazolamide'], `pdb_id` = 3HS4

**Действия системы:** planner → tool_router → chemistry_agent → structure_agent → docking_agent → orc → analysis_agent → writer → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `pdb_structure` = 3HS4 [database, RCSB PDB]
- `resolution` = 1.1 Å [database, RCSB PDB]
- `structure_chain` = A [database, RCSB PDB]
- `metal_ion` = ZN [database, RCSB PDB]
- `vina_score` = -6.362 kcal/mol [calculation, AutoDock Vina]
- `vina_score_pocket1` = -6.362 kcal/mol [calculation, AutoDock Vina]
- `vina_score_pocket2` = -6.278 kcal/mol [calculation, AutoDock Vina]
- `p2rank_pocket` = 29.12 [calculation, P2Rank]

Артефакты: `evidence_csv`, `heatmap_png`, `pockets_csv`, `pose_pdbqt`, `poses_csv`, `receptor_pdb`, `receptor_pdbqt`, `report_json`, `report_md`, `run_meta`

Ограничения, которые система назвала сама:

- Нормализация таутомеров и протонирования под pH не выполнялась.
- Дескрипторы не подтверждают активность, эффективность или безопасность.
- Контроль редокинга пройден: поза №1 в 1.3 Å от кристаллической (AZM), порог 2 Å. Это проверка метода на этой структуре, а не гарантия для других лигандов.
- P2Rank нашёл 2 кармана из запрошенных 3: список не дополнялся выдуманными центрами, докинг выполнен только в найденных.
- Докинг выполнен вслепую: центры боксов взяты у P2Rank по структуре без лиганда, координаты кристаллического лиганда в расчёт не подавались.
- Рецептор протонирован pdb2pqr (AMBER, pH 7.4, PROPKA) и типизирован Meeko; состояния протонирования не оптимизировались под конкретный лиганд.
- Ионы металлов (1) возвращены из исходной структуры: pdb2pqr их отбрасывает. Протонирование остатков, координирующих металл, не уточнялось.
- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- Карманы обработаны параллельными воркерами ORC по позам, которые в этом запуске посчитал AutoDock Vina.
- Замечание LLM-критика, не подтверждённое проверкой: докинг в карманах и анализ лучшей позы по Vina. Шаг для устранения критик не назвал, автоматически это не достраивается.

## Сравнение поз (`compare_docking`)

**Запрос:** Просмотри таблицу универсальных свойств белка, которую мы построили в докинге, и покажи, какие взаимодействия с белком меняются; найди интересные случаи, где свойство меняется неожиданно или есть закономерность в падениях, подъёмах и всплесках, и составь аналитику.

**Ожидаемый ответ:** Разница оценок Vina относительно позы №1 в ккал/моль и график по позам. Карта «поза × остаток»: где появляются и исчезают водородные связи, солевые мостики и другие контакты. График расстояний в Å между одними и теми же парами атомов во всех позах с подписями атомов. Таблицы protein_state.csv и pose_contacts.csv и вывод, какую позу проверить подробнее и по какой причине.

**Контекст из предыдущего шага цепочки:** `molecules` = ['acetazolamide'], `pdb_id` = 3HS4, `poses_csv_path` = outputs/docking/3HS4_acetazolamide/poses.csv, `pose_file` = outputs/docking/3HS4_acetazolamide/poses/pocket1_pose1.pdbqt

**Действия системы:** planner → tool_router → analysis_agent → pose_comparison_agent → chemistry_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `vina_delta_to_best_pose` = 0.226 kcal/mol [calculation, poses.csv]
- `contact_gained_in_pose` = Asn62 [calculation, pose_contacts.csv]
- `contact_lost_in_pose` = Gln92 [calculation, pose_contacts.csv]
- `key_distance_by_pose` = 2.231 Å [calculation, protein_state.csv]
- `pose_to_inspect` = pocket1/поза 8 [hypothesis, ollama:qwen2.5:7b-instruct]
- `best_vina_acetazolamide_pocket1` = -6.362 kcal/mol [calculation, poses.csv]
- `best_vina_acetazolamide_pocket2` = -6.278 kcal/mol [calculation, poses.csv]
- `vina_score_by_pose` = -6.362 kcal/mol [calculation, poses.csv]

Артефакты: `contact_map_png`, `evidence_csv`, `heatmap_png`, `key_distances_png`, `pose_contacts_csv`, `poses_csv`, `protein_state_csv`, `report_json`, `vina_by_pose_png`

Ограничения, которые система назвала сама:

- Тепловая карта строится по лучшему Vina score на пару лиганд×карман.
- На графиках показаны 10 поз из 18: карман позы №1 целиком и лучшая поза каждого из остальных карманов. Полные таблицы содержат все позы.
- Контакты определяются по геометрии тяжёлых атомов и свойствам элементов без явного расчёта водородов: пороги 3,5 Å для водородной связи, 4,0 Å для солевого мостика и гидрофобного контакта, 2,8 Å для координации металла.
- Сравниваются позы одного прогона при одинаковых настройках скоринга (vina, exhaustiveness=8); разница оценок Vina не является разницей измеренной аффинности.
- RMSD до кристаллической позы показан как контрольная метрика и в выбор позы не входит.
- Молекулярная масса, logP и QED относятся к лиганду и не меняются при смене кармана.
- Замечание LLM-критика, не подтверждённое проверкой: таблица универсальных свойств белка и аналитика. Шаг для устранения критик не назвал, автоматически это не достраивается.

## Проверка геометрии позы №1 (`validate_pose`)

**Запрос:** Проверь геометрию позы №1 и покажи, какие критерии она прошла.

**Ожидаемый ответ:** Вердикт по позе, число пройденных проверок из числа применимых и процент пройденности. По каждой проверке — измеренное значение и порог: пара атомов, остаток белка и расстояние для столкновения; значение и допустимый диапазон для связи или угла; отклонение от кристалла для контакта с металлом. Критические нарушения перечислены отдельно. Поза берётся только из расчёта этого же запроса.

**Контекст из предыдущего шага цепочки:** `molecules` = ['acetazolamide'], `pdb_id` = 3HS4, `poses_csv_path` = outputs/docking/3HS4_acetazolamide/poses.csv, `pose_file` = outputs/docking/3HS4_acetazolamide/poses/pocket1_pose1.pdbqt

**Действия системы:** planner → tool_router → structure_agent → validation_agent → orchestrator → critic → review_agent

**Планировщик:** `ollama:qwen2.5:7b-instruct`

**Фактический результат:** `success`

Ключевые значения:

- `pdb_structure` = 3HS4 [database, RCSB PDB]
- `resolution` = 1.1 Å [database, RCSB PDB]
- `structure_chain` = A [database, RCSB PDB]
- `metal_ion` = ZN [database, RCSB PDB]
- `geometry_verdict` = поза принята [calculation, pose_geometry]
- `checks_passed_percent` = 100.0 % [calculation, pose_geometry]
- `geometry_checks_passed` = 6 проверок [calculation, pose_geometry]
- `bound_ligand` = AZM [database, RCSB PDB]

Артефакты: `evidence_csv`, `geometry_report_json`, `geometry_report_md`, `report_json`
