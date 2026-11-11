# Архитектура DEARMO

## Паттерн

Иерархический конвейер: **GigaChat-планировщик** (с keyword fallback) → **детерминированные handlers** → **review_agent** (программная проверка Evidence).

## Роли агентов

| Агент | Назначение |
|-------|------------|
| planner | Сценарий и сущности (JSON) |
| data_agent | BindingDB, PubChem, RAG |
| structure_agent | PDB метаданные |
| docking_agent | Загрузка/запуск докинга |
| chemistry_agent | RDKit |
| analysis_agent | Сравнение карманов, графики |
| validation_agent | Геометрия позы |
| review_agent | Контракт Evidence |

## Потоки данных

1. NL-запрос → `ResearchRequest` → `WorkflowState`.
2. Инструменты пишут в `state.evidence` и `state.artifacts`.
3. `ResearchResult` + `agent_trace` печатает `python cli.py ask`.

Планировщик собирает **список** сценариев. `orchestration/executor.py` идёт по ним в порядке мишень → докинг → сравнение → валидация, чинит недостающий PDB из рекомендации и повторяет шаг один раз.

## Внешние зависимости

- **Обязательные:** RDKit, PubChem (с fallback), локальные кеши BindingDB/PDB.
- **Опциональные:** GigaChat API (`GIGACHAT_CREDENTIALS`), P2Rank/Vina binaries, PoseBusters.

## Каталоги артефактов

- `data/examples/acetazolamide/` — демо докинг.
- `outputs/` — пользовательские прогоны и графики (в .gitignore).
