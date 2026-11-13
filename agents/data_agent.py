"""Агент данных: интеграция PubChem, BindingDB, ChEMBL, сводной таблицы кандидатов и литературного RAG с метаданными кэша."""
from __future__ import annotations

import re

from models.result import Evidence
from orchestration.state import WorkflowState
from tools.bindingdb_client import fetch_activities as fetch_bindingdb_activities
from tools.chembl_client import (
    deduplicate_activities,
    fetch_activities as fetch_chembl_activities,
    reset_transport_status,
    transport_status,
)
from tools.literature_rag import retrieve
from tools.pubchem_client import resolve_compound, search_similar
try:
    from tools.ranking import build_candidates_comparison_table
except ImportError:
    build_candidates_comparison_table = None  # type: ignore[assignment]


def collect_target_data(state: WorkflowState) -> WorkflowState:
    """Собрать сырые данные по соединению: PubChem и измерения BindingDB.

    Мишень здесь не выбирается. Выбор — отдельный шаг (`target_agent`), и он
    идёт после загрузки, чтобы решение принималось по полной выдаче измерений,
    а не по той её части, которую успели отфильтровать.
    """
    with state.agent_step("data_agent", "collect_activities"):
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
        bdb = fetch_bindingdb_activities(ligand, smiles=state.observations.get("smiles"))
        state.observations["bindingdb"] = bdb
        state.trace(
            "data_agent",
            f"bindingdb:{len(bdb.get('activities') or [])} измерений",
            detail=bdb.get("source"),
        )
        return state


def report_activities(state: WorkflowState) -> WorkflowState:
    """Доложить измерения выбранной мишени и дополнить их ChEMBL."""
    with state.agent_step("data_agent", "report_activities"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else None
        )
        bdb = state.observations.get("bindingdb") or {}
        bdb_is_cached = bool(bdb.get("is_cached") or bdb.get("_from_cache") or bdb.get("_from_download"))
        bdb_fetched_at = bdb.get("fetched_at")
        bdb_activities = bdb.get("activities") or []
        recommended = bdb.get("recommended_target") if isinstance(bdb.get("recommended_target"), dict) else None

        chosen_target = state.observations.get("chosen_target") or {}
        selected_bdb = _select_activity_rows(
            bdb_activities,
            target_hint=state.observations.get("target_hint") or state.request.target_id,
            organism_hint=chosen_target.get("organism"),
            recommended_name=(recommended or {}).get("name"),
        )
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

        rt: dict | None = None
        if bdb.get("recommended_target"):
            rt = dict(bdb["recommended_target"])
            rt.setdefault("is_cached", bdb_is_cached)
            rt.setdefault("fetched_at", bdb_fetched_at)
            # Без target_name факт показывался как «Ki 3.1 нМ / None»: имя мишени
            # в этой записи лежит в поле name.
            rt.setdefault("target_name", rt.get("name"))
            state.observations["recommended_target"] = rt
            # Мишень выбирает target_agent. Эвристика клиента BindingDB («лучшая
            # измеренная строка») рядом с его выбором давала в одном ответе две
            # разные мишени: CA3 из клиента и CA2 из агента. Факт остаётся в
            # наблюдениях как запасной вариант, но в Evidence не попадает.
            if state.observations.get("chosen_target"):
                rt = None
        if rt is not None:
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

        # 2. ChEMBL: названная мишень или активности молекулы. CA II не подставляется по умолчанию.
        target_id = state.request.target_id or state.observations.get("target_hint") or ""
        state.observations["chembl_query_mode"] = "target" if target_id else "molecule"
        reset_transport_status()
        chembl_rows = fetch_chembl_activities(target_id=target_id, molecule_id=ligand)
        chembl_transport = transport_status()
        if not chembl_rows and chembl_transport.get("reachable") is False:
            note = (
                f"ChEMBL сейчас недоступен ({chembl_transport.get('detail')}): его активности в ответ не вошли. "
                "Это не значит, что измерений нет — источник не ответил."
            )
            if note not in state.limitations:
                state.limitations.append(note)
            state.observations["chembl_unavailable"] = True
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

        # 3. Отсутствие данных фиксируется явно и с объяснением, что именно нашлось.
        wanted_target = str(state.observations.get("target_hint") or state.request.target_id or "").strip()
        if not selected_bdb and not chembl_rows:
            state.limitations.append(
                f"Экспериментальные данные активности для '{ligand}' в BindingDB и ChEMBL не найдены "
                "(отсутствие данных зафиксировано явно как None, без подмены на 0)."
            )
            state.observations["activity_data_missing"] = True
        elif wanted_target and not _rows_for_target(selected_bdb, wanted_target):
            # Кейс требует показать, что именно нашлось, когда точной пары нет.
            found = sorted({
                f"{row.get('target_name')} ({row.get('organism')})"
                for row in bdb_activities
                if row.get("target_name")
            })
            state.observations["exact_pair_missing"] = True
            state.limitations.append(
                f"Измерений для пары «{ligand} — {wanted_target}» не нашлось. "
                + (
                    "Для этого соединения в BindingDB есть измерения по другим мишеням: "
                    + ", ".join(found[:6])
                    + ("…" if len(found) > 6 else "")
                    + "."
                    if found
                    else "Измерений по этому соединению в BindingDB нет вовсе."
                )
            )

        # 4. Если передано несколько молекул — строим сводную таблицу сравнения кандидатов
        if len(state.request.molecules) > 1 and build_candidates_comparison_table is not None:
            table_kwargs = {"constraints": state.request.constraints}
            if target_id:
                table_kwargs["target_id"] = target_id
            state.observations["candidates_table"] = build_candidates_comparison_table(
                state.request.molecules,
                **table_kwargs,
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
        hits = list(retrieve(_rag_terms(state), limit=3))
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


def _rows_for_target(rows: list[dict], target_name: str) -> list[dict]:
    wanted = target_name.strip().lower()
    if not wanted:
        return list(rows)
    return [
        row
        for row in rows
        if wanted in (row.get("target_name") or "").lower()
        or (row.get("target_name") or "").lower() in wanted
    ]


def _select_activity_rows(
    activities: list[dict],
    *,
    target_hint: str | None,
    recommended_name: str | None,
    organism_hint: str | None = None,
) -> list[dict]:
    """Строки той мишени, которую назвали. Иначе — строки рекомендованной мишени, иначе сильнейшие измерения.

    Организм входит в условие отбора. Без него «Carbonic anhydrase 3» для
    человека собирала и бычьи измерения, и в ответ про человеческий белок
    попадала константа Ki 83400 нМ из другого организма — число настоящее,
    но относящееся не к той мишени.
    """
    hint = str(target_hint or "").strip().lower()
    organism = str(organism_hint or "").strip().lower()
    if hint:
        matched = [row for row in activities if hint in (row.get("target_name") or "").lower()]
        if matched and organism:
            same_organism = [
                row for row in matched if (row.get("organism") or "").strip().lower() == organism
            ]
            if same_organism:
                matched = same_organism
        if matched:
            return matched[:8]
    name = str(recommended_name or "").strip().lower()
    if name:
        matched = [
            row
            for row in activities
            if name in (row.get("target_name") or "").lower()
            or (row.get("target_name") or "").lower() in name
        ]
        if matched:
            return matched[:8]
    numeric: list[tuple[float, dict]] = []
    for row in activities:
        if row.get("relation") not in (None, "=", ""):
            continue
        try:
            numeric.append((float(row["value"]), row))
        except (TypeError, ValueError, KeyError):
            continue
    numeric.sort(key=lambda item: item[0])
    if numeric:
        return [row for _, row in numeric[:8]]
    return activities[:8]


def _rag_terms(state: WorkflowState) -> list[str]:
    query = state.request.text.lower()
    terms = re.findall(r"[a-zа-яё0-9]{4,}", query)
    for extra in (
        state.observations.get("ligand_query"),
        state.observations.get("target_hint"),
        state.observations.get("pdb_id"),
    ):
        if extra:
            token = str(extra).lower()
            if token not in terms:
                terms.append(token)
    return terms


attach_literature_rag = run_rag


def analyze(state: WorkflowState) -> WorkflowState:
    """Полная цепочка данных одним вызовом: сбор, выбор мишени, отчёт, литература."""
    from agents import target_agent  # локальный импорт: target_agent не зависит от data_agent

    state = collect_target_data(state)
    if state.errors:
        return state
    state = target_agent.select_target(state)
    state = report_activities(state)
    return run_rag(state)
