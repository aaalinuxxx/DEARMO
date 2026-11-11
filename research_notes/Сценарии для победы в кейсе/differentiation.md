# Сценарии, которые ещё не закрыты открытыми инструментами

Проверено 11 ноября 2026. HTML-страницы GitHub отдали только описание репозитория, без тела README. Утверждения README ниже сняты с raw-файла ветки по умолчанию. Лицензия — с поля `license` GitHub API, а где API расходится с файлом, указаны оба источника. Локальный обзор `repos.pdf` использован как проверяемый черновик, не как источник фактов о репозиториях.

Защита по календарю кейса — 11–13 ноября 2026, то есть у команды двое суток от даты проверки. Кейс 1 не требует докинга: специализация на одном направлении допустима. Он требует не меньше трёх сценариев, которые прототип реально выполняет, запрос на естественном языке, координацию агентов, LLM как компонент управления и явное разделение расчётов, предсказаний и гипотез без вымышленных результатов.

## Что README каждого проекта называет пользовательским сценарием и какая у него лицензия?

### Takeaway

Ни один из семи README не обещает диалог, который уточняет неоднозначный запрос, отказывается смешивать несопоставимые измерения и помечает каждое число как запись базы, расчёт или гипотезу, включая кейс без уверенного ответа. Ближе всего к естественного языку — Prompt-to-Pill и dock_assist, к проверяемому следу — delta и AutoDock-Pipeline, к готовому P2Rank+Vina — PocketDock и AutoDock-Pipeline. Лицензии: MIT у AutoDock-Pipeline, delta, PocketDock, dock_assist и AgentD; у CLADD файл LICENSE заявляет Apache-2.0, а GitHub API ставит NOASSERTION; у Prompt-to-Pill лицензии в репозитории нет.

### Cited Findings

- Prompt-to-Pill описывает пользовательский сценарий так: после запуска MCP-серверов человек вводит задачу в терминале, и система сама проходит генерацию молекул, ADMET, докинг, оптимизацию, генерацию протокола испытания и подбор пациентов. Заявлена «ask-driven» симуляция разработки и осуществимости испытания для терапевтической мишени. Докинг в таблице файлов назван AutoDock Vina с поиском кармана P2Rank; отдельный MCP-сервер считает дескрипторы RDKit. В корне нет файла LICENSE, поле `license` API равно `null`, в README раздела о лицензии нет. — [README](https://raw.githubusercontent.com/ChatMED/Prompt-to-Pill/main/README.md); [GitHub API](https://api.github.com/repos/ChatMED/Prompt-to-Pill); [список корня](https://api.github.com/repos/ChatMED/Prompt-to-Pill/contents/)
- Имя файла в README Prompt-to-Pill, `chemical_properties_mcp_server.py`, не совпадает с файлом в корне `chemical_properties_mcp_sever.py`. — [README](https://raw.githubusercontent.com/ChatMED/Prompt-to-Pill/main/README.md); [список корня](https://api.github.com/repos/ChatMED/Prompt-to-Pill/contents/)
- AutoDock-Pipeline — один запуск `python run_pipeline.py`. Вход: PDB рецептора и список PubChem CID, не фраза на естественном языке. Цепочка: подготовка, P2Rank, бокс, пакетный AutoDock Vina, ранжирование, ProLIF, отчёты. Отдельно заявлены redocking сокристаллизованного лиганда с RMSD, отказ считать прогон успешным, если анализ взаимодействий не удался, и манифест версий и хешей. Скрипт не утверждает, что скор Vina есть экспериментальная аффинность. Лицензия API: MIT. — [README](https://raw.githubusercontent.com/RahulRoktim/AutoDock-Pipeline/main/README.md); [GitHub API](https://api.github.com/repos/RahulRoktim/AutoDock-Pipeline)
- delta описывает многоагентный цикл оптимизации молекул против AKT1: Principal Researcher задаёт цель, агенты базы, генерации, медицинской химии с вызовом докинга, ранжирования и критики обмениваются записями происхождения. Запуск — `src/demo/multi_agent_pipeline.py` с флагами провайдера LLM и числа итераций, не свободный диалог. Инструменты README: UniProt, PubMed, PDB, ChEMBL, Vina, DiffDock, PLIP, QED, SA, logP, молекулярная масса. P2Rank в списке инструментов README нет. Лицензия и в README, и в API: MIT. — [README](https://raw.githubusercontent.com/deltawave-tech/delta/main/README.md); [GitHub API](https://api.github.com/repos/deltawave-tech/delta)
- PocketDock — веб-платформа: белок и лиганд или библиотека до 100 молекул, P2Rank, AutoDock Vina, дескрипторы RDKit на каждой задаче, опционально OpenMM и MM-GBSA, просмотр в браузере. Пример в репозитории: EGFR PDB 4HJO и erlotinib. LLM, агенты и запрос на естественном языке в перечне функций README не названы; Claude упомянут как помощь при написании кода. Лицензия README и API: MIT. — [README](https://raw.githubusercontent.com/gozsari/PocketDock/main/README.md); [GitHub API](https://api.github.com/repos/gozsari/PocketDock)
- dock_assist в README назван Modrag-dock: один tool-using LLM-агент. Пример запроса: «dock dopamine into SULT1A3». Агент ищет белок в PDB, скачивает структуру, находит карман сеткой buriedness без P2Rank, fpocket и SiteMap и докает AutoDock Vina. Есть докинг по известному центроиду и по металлу. Описание репозитория на GitHub: «AI-assisted blind docking for SMILES strings on any protein». Лицензия API: MIT; в прочитанной верхней части README отдельного раздела License нет. — [README](https://raw.githubusercontent.com/MauricioCafiero/dock_assist/main/README.md); [GitHub API](https://api.github.com/repos/MauricioCafiero/dock_assist)
- CLADD — мультиагентный RAG для drug discovery без доменного дообучения. Документированный запуск — не свободная фраза, а скрипты на датасетах: `caption_multi_agent.py --dataset bbbp` и `multi_agent.py --dataset hERG`. Нужны подготовленные PubChem, PrimeKG и веса GraphMVP. P2Rank и Vina в README не названы; в зависимостях есть RDKit и LangGraph. Файл LICENSE начинается с Apache License 2.0. Поле GitHub API: `other` / SPDX `NOASSERTION`. — [README](https://raw.githubusercontent.com/Genentech/CLADD/main/README.md); [LICENSE](https://raw.githubusercontent.com/Genentech/CLADD/main/LICENSE); [GitHub API](https://api.github.com/repos/Genentech/CLADD)
- AgentD README описывает пакет модульных агентов для извлечения литературы, предсказания свойств и генерации молекул, плюс MCP-пайплайн и режим Q&A. Запуск пайплайна читает YAML: белок, болезнь, число итераций, REINVENT, фильтры QED и предсказанного pKd, опционально 3D через Boltz. Ветка по умолчанию — `release`, не `main`. Команда клонирования в README указывает `hoon-ock/llm-dd`, а не `hoon-ock/AgentD`. P2Rank и AutoDock Vina в этом README не названы. Лицензия README и API: MIT. — [README ветки release](https://raw.githubusercontent.com/hoon-ock/AgentD/release/README.md); [GitHub API](https://api.github.com/repos/hoon-ock/AgentD)
- Кейс 1 хакатона «СИНТЕЗ», стр. 2–3 PDF: прототип МАС на LLM; запрос на естественном языке; координация агентов; не меньше трёх сценариев вида «запрос → действия → ожидаемый ответ», отражающих реальные возможности; прозрачное различие реальных расчётов, предсказаний и гипотез; отсутствие вымышленных результатов. Докинг перечислен как возможное направление, список не обязателен, система может специализироваться на одном направлении. — [СИНТЕЗ_кейсы.pdf](file:///Users/forallxyz/DEARMO/%D0%A1%D0%98%D0%9D%D0%A2%D0%95%D0%97_%D0%BA%D0%B5%D0%B9%D1%81%D1%8B.pdf)
- Локальный обзор `repos.pdf` в целом совпадает с README: Prompt-to-Pill как LLM плюс P2Rank и Vina; AutoDock-Pipeline как пайплайн без мультиагентности; delta как мультиагентный каркас с Vina и PLIP; PocketDock как готовый P2Rank+Vina без проверенной LLM-МАС; dock_assist как один агент и не P2Rank; CLADD как RAG, не докинг; лицензию CLADD обзор не раскрывает и просит читать файл отдельно. Обзор прямо пишет, что проекты локально не запускались. — [repos.pdf](file:///Users/forallxyz/DEARMO/repos.pdf)

### Inferences

- Для жюри «уже сделано открытым кодом» закрыты отдельно: пакетный P2Rank+Vina с ADMET (PocketDock, AutoDock-Pipeline), естественный язык → Vina (dock_assist и Prompt-to-Pill), мультиагентная оптимизация со следом вызовов (delta) и RAG-ответы по свойствам и токсичности (CLADD). Не закрыта связка «уточнить или отказаться» плюс подпись типа каждого числа в ответе пользователю.
- Копировать код Prompt-to-Pill в решение команды нельзя обосновать лицензией: объявления нет. MIT-репозитории можно использовать с сохранением уведомления; Apache-2.0 у CLADD следует из текста LICENSE, даже хотя GitHub не проставил SPDX.

### Gaps

- Не проверялось, есть ли LICENSE у Prompt-to-Pill вне корня ветки `main` или в коммите после списка contents. API на момент запроса вернул `license: null`.
- Полный текст Apache-2.0 у CLADD не сверялся побайтно с каноническим текстом: прочитан заголовок файла LICENSE, который называет Apache License, Version 2.0. Почему API ставит NOASSERTION, не выяснялось.
- README dock_assist длинный; лицензионный абзац ниже раздела про blind docking в выгрузке не встретился. Вывод «MIT» опирается на GitHub API, не на цитату из README.
- Репозиторий `hoon-ock/llm-dd` из команды clone AgentD отдельно не открывался. Нельзя сказать, совпадает ли он с веткой `release` репозитория AgentD.
- Работоспособность семи проектов не запускалась. Обзор `repos.pdf` тоже этого не делал.

## Какие проекты — один пайплайн, какие — мультиагентные, и кто уже соединяет вызов инструментов LLM с P2Rank и Vina?

### Takeaway

Один скриптовый или веб-пайплайн без заявленных агентов: AutoDock-Pipeline и PocketDock; оба уже соединяют P2Rank и Vina, но не с LLM. Мультиагентные по README: Prompt-to-Pill, delta, CLADD. Один tool-using агент: dock_assist. AgentD описывает модульных агентов и фиксированный пайплайн, не коллегию ролей в духе delta. Единственный README, который прямо собирает LLM-вызываемые инструменты, P2Rank и Vina в одном сценарии, — Prompt-to-Pill. delta вызывает Vina и DiffDock без P2Rank. dock_assist вызывает Vina без P2Rank.

### Cited Findings

- Prompt-to-Pill: Orchestrator и Planning Agent, специализированные агенты, инструменты как MCP-серверы; докинг — Vina и P2Rank. Планировщик выбирает агентов под задачу. — [README](https://raw.githubusercontent.com/ChatMED/Prompt-to-Pill/main/README.md)
- AutoDock-Pipeline: «single workflow», `run_pipeline.py`, внешние бинарники Vina и P2Rank. LLM в README нет. — [README](https://raw.githubusercontent.com/RahulRoktim/AutoDock-Pipeline/main/README.md)
- delta: несколько ролей, схемы инструментов для LLM, `vina_docker` и `diffdock`. P2Rank в перечне Available Tools нет. — [README](https://raw.githubusercontent.com/deltawave-tech/delta/main/README.md)
- PocketDock: Django, Celery, шаги P2Rank → Meeko → RDKit → Vina. Агентов и tool calling в архитектуре README нет. — [README](https://raw.githubusercontent.com/gozsari/PocketDock/main/README.md)
- dock_assist: «conversational agent loop», девять инструментов у одной LLM, карман — собственная сетка, докинг — Vina. README прямо говорит, что внешний инструмент кармана не нужен. — [README](https://raw.githubusercontent.com/MauricioCafiero/dock_assist/main/README.md)
- CLADD: «collaboration of multiple LLM agents», каталог `multi_agent/`, LangGraph. Докинг и P2Rank не входят в описанный запуск. — [README](https://raw.githubusercontent.com/Genentech/CLADD/main/README.md)
- AgentD: «modular agents», режим Q&A и отдельный конфигурируемый пайплайн extract → pool → predict/refine → Boltz. Это последовательный сценарий с LLM, а не описание нескольких спорящих ролей. P2Rank и Vina не названы. — [README](https://raw.githubusercontent.com/hoon-ock/AgentD/release/README.md)
- Локальный код DEARMO: RDKit считает молекулярную массу, logP, TPSA, HBD, HBA и вращаемые связи и кладёт их в `Evidence` с `kind="calculation"`. Статусы результата уже включают `needs_clarification` и `no_data`, виды свидетельств — `database`, `calculation`, `prediction`, `hypothesis`. Планировщик, агенты химии, данных и проверки, клиенты ChEMBL и PubChem, ранжирование и оркестрация бросают `NotImplementedError`. Черновик сценариев — три TODO: анализ, сравнение, поиск аналогов. — [README проекта](file:///Users/forallxyz/DEARMO/README.md); [rdkit_tools.py](file:///Users/forallxyz/DEARMO/tools/rdkit_tools.py); [result.py](file:///Users/forallxyz/DEARMO/models/result.py); [scenarios.md](file:///Users/forallxyz/DEARMO/docs/scenarios.md)

### Inferences

- Показывать жюри «наш агент вызвал P2Rank и Vina» — повторение заявленного сценария Prompt-to-Pill. Показывать таблицу скоров и ADMET — повторение PocketDock и AutoDock-Pipeline.
- Журнал вызовов сам по себе уже заявлен delta. Отличие должно быть видно в реплике пользователю, а не только в том, что инструменты логируются.

### Gaps

- Не читался исходный код планировщика Prompt-to-Pill и delta, только README. Нельзя подтвердить, задают ли они уточняющий вопрос при дырявом запросе: README этого не обещает.
- Не проверялось, вызывает ли AgentD Vina в ноутбуках `example/test_case/`, которых README только перечисляет. В тексте README Vina и P2Rank отсутствуют.

## Какое взаимодействие останется отличимым, если научный движок заимствован?

### Takeaway

Отличимый слой — не новый скорер, а поведение координации, которого нет в проверенных README: уточняющий вопрос, отказ смешивать несопоставимые измерения, подпись каждого числа как база / расчёт / гипотеза и явный ответ «уверенного результата нет». Ближайшие частичные аналоги — научные оговорки AutoDock-Pipeline и записи происхождения delta — не закрывают этот пользовательский диалог.

### Cited Findings

- Кейс 1 требует интерпретируемый ответ и различие расчётов, предсказаний и гипотез, запрещает вымышленные результаты и требует, чтобы сценарии совпадали с тем, что прототип умеет. Отдельного обязательного «негативного кейса» в тексте кейса 1 нет. Такой сценарий есть у кейса 2 другого трека: «неоднозначное происхождение» при недостатке данных. — [СИНТЕЗ_кейсы.pdf](file:///Users/forallxyz/DEARMO/%D0%A1%D0%98%D0%9D%D0%A2%D0%95%D0%97_%D0%BA%D0%B5%D0%B9%D1%81%D1%8B.pdf)
- Prompt-to-Pill после ввода задачи «automatically execute the entire workflow» и показывает фазы по очереди. Уточнение, отказ и кейс без ответа в README не описаны. — [README](https://raw.githubusercontent.com/ChatMED/Prompt-to-Pill/main/README.md)
- dock_assist просит белок и лиганды и выполняет пример «dock dopamine into SULT1A3». Это один успешный путь tool calling. README не описывает отказ смешивать assay и классификацию каждого числа. — [README](https://raw.githubusercontent.com/MauricioCafiero/dock_assist/main/README.md)
- delta хранит вызовы инструментов и provenance молекул и сравнивает мультиагентный режим с одним агентом на AKT1. Это аудит оптимизации, которая выдаёт молекулы и скоры, а не демонстрация отказа при несопоставимых данных. — [README](https://raw.githubusercontent.com/deltawave-tech/delta/main/README.md)
- AutoDock-Pipeline не подменяет провал ProLIF нулями, не выдаёт скрининг без сообщения, воспроизводится ли кристаллографическая поза, и прямо отделяет скор Vina от экспериментальной аффинности. Это шлюз внутри пайплайна без диалога на естественном языке. — [README](https://raw.githubusercontent.com/RahulRoktim/AutoDock-Pipeline/main/README.md)
- PocketDock на каждой задаче считает RDKit-дескрипторы и комбинирует вероятность кармана с аффинностью. README не описывает вопрос пользователю и не помечает числа как гипотезы. — [README](https://raw.githubusercontent.com/gozsari/PocketDock/main/README.md)
- Шаблоны ролей DEARMO уже говорят планировщику запросить уточнение при неоднозначности и всем агентам не выдумывать значения, разделяя записи базы, расчёты, прогнозы и гипотезы. Это локальные черновики промптов, не поведение чужих README. — [planner.txt](file:///Users/forallxyz/DEARMO/prompts/planner.txt); [review_agent.txt](file:///Users/forallxyz/DEARMO/prompts/review_agent.txt)

### Inferences

- Если движок чужой, жюри труднее свести демо к «обёртке», когда на экране происходит решение, которого PocketDock не принимает: система не выпускает число. Четыре приёма усиливают друг друга. Самый сильный одиночный — кейс без уверенного ответа, потому что успешный скор уже умеют PocketDock, AutoDock-Pipeline, Prompt-to-Pill и dock_assist.
- Уточняющий вопрос и отказ смешивать assay стоит показывать на малых молекулах, где RDKit уже отвечает честно. Негативный докинг тогда не выглядит отключенной кнопкой: химия посчитана и подписана, а сродство не выдумано.
- Называть этот негативный сценарий требованием кейса 1 нельзя. Это способ выполнить запрет на вымысел и одновременно не повторить успешный докинг из README. Кейс 2 использовать как цитату требования для кейса 1 некорректно.

### Gaps

- Не найден README из этой семёрки, где автор прямо пишет, что система задаёт уточняющий вопрос или отказывается усреднять IC50 с Ki. Отсутствие фразы в README не доказывает, что поведение невозможно в коде.
- Не сравнивались статьи Prompt-to-Pill (bioRxiv, DOI в README) и delta (arXiv 2508.03444) сверх того, что продублировано в README. В статье может быть негативный пример, которого нет в README.

## Какая тройка сценариев повышает шанс при одном работающем RDKit: три химии, три докинга или гибрид?

### Takeaway

Выбирать гибрид: два сценария на уже работающем RDKit и третий — докинг как честный отказ или уточнение, без выдуманного скора. Три белок-лигандных сценария за двое суток при `NotImplemented` по агентам и API проигрывают и по исполнимости, и по новизне: успешный P2Rank+Vina уже описан у PocketDock, AutoDock-Pipeline и Prompt-to-Pill. Три только химических сценария кейс формально допускает, но их легче всего назвать калькулятором RDKit, особенно если третий — снова дескрипторы. Поиск аналогов из локального черновика сейчас нечем считать: ChEMBL, PubChem и ранжирование не реализованы.

### Cited Findings

- В DEARMO работает разбор SMILES и шесть дескрипторов RDKit с происхождением `calculation`, версией RDKit и оговоркой, что дескрипторы не подтверждают активность. Нормализация солей, зарядов и таутомеров не делается. Ошибка пустого или невалидного SMILES возвращается как `status="error"`. — [README проекта](file:///Users/forallxyz/DEARMO/README.md); [rdkit_tools.py](file:///Users/forallxyz/DEARMO/tools/rdkit_tools.py)
- Не реализованы: `chembl_client`, `pubchem_client`, `ranking`, `chemistry_agent`, `data_agent`, `review_agent`, `planner`, `workflow`. Интерфейс пишет, что LLM, агенты и базы не подключены. — [README проекта](file:///Users/forallxyz/DEARMO/README.md); [app.py](file:///Users/forallxyz/DEARMO/app.py)
- Кейс требует нескольких агентов, LLM, внешний инструмент или базу, три живых сценария и запрет вымышленных результатов. Одно направление разрешено. — [СИНТЕЗ_кейсы.pdf](file:///Users/forallxyz/DEARMO/%D0%A1%D0%98%D0%9D%D0%A2%D0%95%D0%97_%D0%BA%D0%B5%D0%B9%D1%81%D1%8B.pdf)
- Готовый успешный докинг без LLM уже выложен: PocketDock, пример EGFR/erlotinib; AutoDock-Pipeline, CID → Vina. Успешный естественный язык плюс Vina уже выложен: dock_assist и Prompt-to-Pill, причём второй ещё и с P2Rank. — [PocketDock README](https://raw.githubusercontent.com/gozsari/PocketDock/main/README.md); [AutoDock-Pipeline README](https://raw.githubusercontent.com/RahulRoktim/AutoDock-Pipeline/main/README.md); [dock_assist README](https://raw.githubusercontent.com/MauricioCafiero/dock_assist/main/README.md); [Prompt-to-Pill README](https://raw.githubusercontent.com/ChatMED/Prompt-to-Pill/main/README.md)
- Локальные критерии сравнения, сопоставимость assays и паспорт мишени не заполнены: в `docs/selection_criteria.md` и `docs/target_description.md` стоят TODO. — [selection_criteria.md](file:///Users/forallxyz/DEARMO/docs/selection_criteria.md); [target_description.md](file:///Users/forallxyz/DEARMO/docs/target_description.md)

### Inferences

- Тройка, которую труднее всего свести к PocketDock или к калькулятору RDKit, собирается из того, что уже считается, и из того, что система сознательно не считает.
- Сценарий A, химия, исполнимый сейчас. Запрос со SMILES вроде этанола `CCO`. Ответ: шесть дескрипторов, у каждого подпись «расчёт / RDKit», канонический SMILES, версия, ограничение «это не активность». Это нужно как доказательство инструмента, но в одиночку это калькулятор.
- Сценарий B, химия, тоже исполнимый без баз. Запрос смешивает уровни: «у кого лучше logP и кто сильнее ингибирует» или даёт имя без SMILES. Ожидаемое поведение: уточняющий вопрос, если нет структуры; если есть только дескрипторы — отказ ранжировать активность и запрет подставлять assay. Числа logP остаются расчётами, фраза про активность не появляется. Это не умеет ни таблица PocketDock, ни функция RDKit сама по себе.
- Сценарий C, докинг как граница, а не как скрининг. Запрос «пристыкуй лиганд и назови скор». Пока Vina и P2Rank не подключены, фактический результат: статус `needs_clarification` или `no_data`, ни одного сродства, перечисление, каких данных не хватает (структура, карман, определение скора). Если SMILES валиден, дескрипторы из A можно приложить отдельно и подписать как расчёт. Такой сценарий отражает реальную возможность и не копирует демо EGFR у PocketDock.
- Три полноценных докинга хуже этой тройки: за двое суток их нельзя честно посчитать текущим кодом, а выдуманный ккал/моль нарушает кейс. Даже успешный один комплекс выглядит как тонкая обёртка над PocketDock или Prompt-to-Pill.
- Три только химических сценария хуже гибрида по различимости. Черновик «поиск аналогов» упирается в нереализованные PubChem/ChEMBL. Третий расчёт дескрипторов не добавляет координации. Гибрид оставляет два надёжных прогона и один сюжет, который нельзя отмахнуть фразой «это PocketDock» или «это RDKit».
- LLM и несколько агентов кейс всё равно требует. Гибрид не заменяет оркестрацию: планировщик должен выбрать A/B/C, химический агент вызвать RDKit, проверяющий агент не пропустить число без вида свидетельства. Это тонкий слой над уже готовой моделью `ResearchResult`, и он реалистичен за двое суток именно потому, что научный шаг один.

### Gaps

- Не измерялось время установки AutoDock Vina и P2Rank в окружении команды. Вывод «докинг не успеть» следует из статуса репозитория DEARMO и календаря, а не из неудачной попытки сборки.
- Нет протокола жюри с весами критериев. Рекомендация гибрида — вывод из текста кейса и из дыр README, не из опубликованной рубрики оценок.
