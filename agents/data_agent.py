"""Агент данных: интеграция PubChem, BindingDB, ChEMBL, сводной таблицы кандидатов и литературного RAG с метаданными кэша."""
from __future__ import annotations

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.bindingdb_client import fetch_activities as fetch_bindingdb_activities
from tools.chembl_client import deduplicate_activities, fetch_activities as fetch_chembl_activities
from tools.literature_rag import retrieve
from tools.pubchem_client import resolve_compound, search_similar
try:
    from tools.ranking import build_candidates_comparison_table
except ImportError:
    build_candidates_comparison_table = None  # type: ignore[assignment]


def collect_target_data(state: WorkflowState) -> WorkflowState:
    """Собрать данные по соединению и активностям из PubChem, BindingDB и ChEMBL."""
    with state.agent_step("data_agent", "bindingdb_chembl_pubchem"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else None
        )
        if not ligand:
            state.add_error("Не указан лиганд. Уточните название соединения или SMILES.")
            return state

        pubchem = resolve_compound(ligand)
        state.observations["pubchem"] = pubchem

        if pubchem.get("ambiguous"):
            candidates = pubchem.get("candidates") or []
            state.observations["clarification_candidates"] = candidates
            state.add_error(
                pubchem.get("message")
                or f"Уточните название молекулы '{ligand}': найдено несколько различных структур."
            )
            return state

        if pubchem.get("cid"):
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=pubchem["cid"],
                    method="pubchem_cid",
                    source="PubChem",
                    context={
                        **pubchem,
                        "is_cached": bool(pubchem.get("is_cached", False)),
                        "fetched_at": pubchem.get("fetched_at"),
                    },
                )
            )
        if pubchem.get("connectivity_smiles"):
            state.observations["smiles"] = pubchem["connectivity_smiles"]

        # 1. Данные из BindingDB (с метаданными is_cached и fetched_at)
        bdb = fetch_bindingdb_activities(ligand)
        state.observations["bindingdb"] = bdb
        bdb_is_cached = bool(bdb.get("is_cached") or bdb.get("_from_cache") or bdb.get("_from_download"))
        bdb_fetched_at = bdb.get("fetched_at")

        bdb_activities = bdb.get("activities") or []
        ca_rows = [
            act
            for act in bdb_activities
            if "carbonic" in (act.get("target_name") or "").lower()
            or "carbonate dehydratase" in (act.get("target_name") or "").lower()
        ]
        selected_bdb = (ca_rows or bdb_activities)[:8]
        for act in selected_bdb:
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=act.get("value"),
                    unit=act.get("units"),
                    method=act.get("affinity_type", "affinity"),
                    source="BindingDB",
                    context={
                        "target_name": act.get("target_name"),
                        "organism": act.get("organism"),
                        "pmid": act.get("pmid"),
                        "link": act.get("link"),
                        "is_cached": bdb_is_cached,
                        "fetched_at": bdb_fetched_at,
                    },
                )
            )

        if bdb.get("recommended_target"):
            rt = dict(bdb["recommended_target"])
            rt.setdefault("is_cached", bdb_is_cached)
            rt.setdefault("fetched_at", bdb_fetched_at)
            state.observations["recommended_target"] = rt
            state.add_evidence(
                Evidence(
                    kind="database" if not rt.get("hypothesis") else "hypothesis",
                    value=rt.get("value"),
                    unit=rt.get("units"),
                    method=rt.get("affinity_type", "target_recommendation"),
                    source="BindingDB",
                    context=rt,
                )
            )

        # 2. Данные из ChEMBL (нормализованные в нМ)
        target_id = state.request.target_id or state.observations.get("target_hint") or "CHEMBL205"
        chembl_rows = fetch_chembl_activities(target_id=target_id, molecule_id=ligand)
        chembl_rows, chembl_lims = deduplicate_activities(chembl_rows)
        state.observations["chembl_activities"] = chembl_rows
        for lim in chembl_lims:
            if lim not in state.limitations:
                state.limitations.append(lim)

        for row in chembl_rows[:8]:
            if not row.get("comparable"):
                state.limitations.append(
                    f"Несопоставимое измерение ChEMBL ({row.get('standard_type')} в {row.get('standard_units')}) "
                    "помечено как comparable=False и не сравнивается напрямую с концентрационными Ki/IC50 в nM."
                )
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=row.get("value_nm") if row.get("comparable") else row.get("standard_value"),
                    unit=row.get("normalized_units") if row.get("comparable") else row.get("standard_units"),
                    method=f"chembl_{str(row.get('standard_type') or 'activity').lower()}",
                    source="ChEMBL",
                    context=row,
                )
            )

        # 3. Явная фиксация отсутствия данных об активности (не подменяем на 0)
        if not selected_bdb and not chembl_rows:
            state.limitations.append(
                f"Экспериментальные данные активности для '{ligand}' в BindingDB и ChEMBL не найдены "
                "(отсутствие данных зафиксировано явно как None, без подмены на 0)."
            )
            state.observations["activity_data_missing"] = True

        # 4. Если передано несколько молекул — строим сводную таблицу сравнения кандидатов
        if len(state.request.molecules) > 1 and build_candidates_comparison_table is not None:
            state.observations["candidates_table"] = build_candidates_comparison_table(
                state.request.molecules,
                target_id=target_id,
                constraints=state.request.constraints,
            )

        # 5. Опциональный поиск похожих соединений при запросе аналогов
        if state.observations.get("smiles") and (
            "аналог" in state.request.text.lower() or "похож" in state.request.text.lower()
        ):
            sim_res = search_similar(state.observations["smiles"], threshold=80, max_records=5)
            state.observations["similar_compounds"] = sim_res.get("analogs", [])
            for analog in sim_res.get("analogs", [])[:5]:
                state.add_evidence(
                    Evidence(
                        kind="database",
                        value=analog.get("tanimoto_similarity"),
                        unit="Tanimoto",
                        method="pubchem_analog",
                        source="PubChem",
                        context={
                            **analog,
                            "is_cached": bool(sim_res.get("is_cached", False)),
                            "fetched_at": sim_res.get("fetched_at"),
                        },
                    )
                )

        return state


def run_rag(state: WorkflowState) -> WorkflowState:
    """Извлечь релевантные фрагменты литературы из локального RAG."""
    with state.agent_step("data_agent", "literature_rag"):
        query = state.request.text
        hits = list(retrieve(query))[:3]
        for h in hits:
            state.rag_context.append(h)
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=h.get("text", "")[:200],
                    source="RAG",
                    method="semantic_search",
                    context={**h, "is_cached": True, "fetched_at": "local-rag-corpus"},
                )
            )
        lim_note = "Контекст RAG — выдержки из заранее загруженных источников, не полнотекстовый поиск"
        if lim_note not in state.limitations:
            state.limitations.append(lim_note)
        return state


attach_literature_rag = run_rag


def analyze(state: WorkflowState) -> WorkflowState:
    state = collect_target_data(state)
    return run_rag(state)
