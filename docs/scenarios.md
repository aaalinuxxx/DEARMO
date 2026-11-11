# Сценарии DEARMO (кейс СИНТЕЗ №1)

Формат: **запрос → действия агентов → источники → ожидаемый ответ → ограничения**.

## 1. target_for_ligand

- **Запрос:** Подбери белок для докинга ацетазоламида. Обоснуй экспериментальными данными и предложи PDB.
- **Агенты:** `planner` → `chemistry_agent` (SMILES, RDKit) → `data_agent` (BindingDB, PubChem, RAG) → `structure_agent` (PDB) → `review_agent`.
- **Источники:** `data/cache/bindingdb/`, PubChem PUG REST, `data/rag/`, `data/cache/pdb/`.
- **Ожидаемый ответ:** рекомендация CA II, PDB **3HS4**, Ki из BindingDB как `database`, дескрипторы как `calculation`.
- **Ограничения:** Ki ≠ Vina score; при только косвенных данных — `hypothesis`.

## 2. dock_pockets

- **Запрос:** Докинг ацетазоламида в 3HS4, два кармана P2Rank.
- **Агенты:** `planner` → `chemistry_agent` → `structure_agent` → `docking_agent` (load cache) → `review_agent`.
- **Источники:** `data/examples/acetazolamide/` (poses.csv, report.md).
- **Ожидаемый ответ:** лучший `vina_score` (calculation), ссылка на pose file.
- **Ограничения:** на защите — **кешированный прогон**; live — `DOCKING_LIVE=1` и CLI `tools/docking_pipeline.py`.

## 3. compare_docking

- **Запрос:** Сравни три лиганда в двух карманах; меняется ли порядок по Vina?
- **Агенты:** `planner` → `analysis_agent` → `chemistry_agent` (MW/logP по лигандам) → `review_agent`.
- **Источники:** `poses.csv`, matplotlib.
- **Ожидаемый ответ:** heatmap `outputs/analysis/vina_heatmap.png`, признак `order_changes`.
- **Ограничения:** MW/logP не зависят от кармана.

## 4. validate_pose

- **Запрос:** Проверка правдоподобности позы ацетазоламида.
- **Агенты:** `planner` → `structure_agent` → `validation_agent` → `review_agent`.
- **Источники:** PDBQT поза, опционально PoseBusters.
- **Ожидаемый ответ:** `%` пройденных проверок, RMSD к кристаллу при наличии в кеше.

## Контроль redocking

- **3PTB / BEN:** `data/examples/3ptb_ben/poses.csv`, RMSD ≈ 0.356 Å (golden test).
