"""Агент данных: BindingDB, PubChem, RAG."""
from __future__ import annotations

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.bindingdb_client import fetch_activities
from tools.literature_rag import retrieve
from tools.pubchem_client import resolve_compound


def collect_target_data(state: WorkflowState) -> WorkflowState:
    with state.agent_step("data_agent", "bindingdb_pubchem"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else None
        )
        if not ligand:
            state.add_error("Не указан лиганд.")
            return state

        pubchem = resolve_compound(ligand)
        if pubchem.get("cid"):
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=pubchem["cid"],
                    method="pubchem_cid",
                    source="PubChem",
                    context=pubchem,
                )
            )
        if pubchem.get("connectivity_smiles"):
            state.observations["smiles"] = pubchem["connectivity_smiles"]

        bdb = fetch_activities(ligand)
        state.observations["bindingdb"] = bdb
        ca_rows = [
            act for act in (bdb.get("activities") or [])
            if "carbonic" in (act.get("target_name") or "").lower()
            or "carbonate dehydratase" in (act.get("target_name") or "").lower()
        ]
        for act in (ca_rows or (bdb.get("activities") or []))[:8]:
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
                    },
                )
            )
        if bdb.get("recommended_target"):
            rt = bdb["recommended_target"]
            state.observations["recommended_target"] = rt
            state.add_evidence(
                Evidence(
                    kind="database" if not rt.get("hypothesis") else "hypothesis",
                    value=rt.get("name"),
                    method="recommended_target",
                    source="BindingDB",
                    context=rt,
                )
            )
    return state


def attach_literature_rag(state: WorkflowState) -> WorkflowState:
    with state.agent_step("data_agent", "rag_retrieve"):
        terms = list(state.request.molecules or [])
        if state.observations.get("target_hint"):
            terms.append(state.observations["target_hint"])
        chunks = retrieve(terms)
        state.rag_context = chunks
        for ch in chunks:
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=ch.get("title") or ch.get("chunk_id"),
                    method="literature_chunk",
                    source=ch.get("source", "rag"),
                    context={"text": (ch.get("text") or "")[:500], "chunk_id": ch.get("chunk_id")},
                )
            )
        if chunks:
            state.limitations.append("Контекст RAG — выдержки из заранее загруженных источников, не полнотекстовый поиск.")
    return state


def collect_data(state: WorkflowState) -> WorkflowState:
    return collect_target_data(state)
