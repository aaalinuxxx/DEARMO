"""Агент структур: найти кристалл выбранной мишени и описать его состав.

Поиск идёт по accession UniProt, а не по словам: «carbonic anhydrase» в
полнотекстовом индексе RCSB находит сотни записей разных изоформ и организмов,
и выбор из такой выдачи обосновать нечем. Уровень поиска (A–D) возвращается
вместе с результатом и попадает в ответ: он говорит, насколько найденная
структура соответствует задаче.

Модель выбирает только номер варианта из списка, собранного RCSB. Ни один
идентификатор не берётся из её памяти.
"""
from __future__ import annotations

import re
from dataclasses import replace
from typing import Any

from models.result import Evidence
from orchestration.llm_client import decide
from orchestration.state import WorkflowState
from tools.pdb_client import (
    find_chem_comp_ids,
    get_entity_details,
    get_structure_metadata,
    search_structures,
    significant_cofactors,
)

_MAX_CANDIDATES = 5

# Чем выше уровень поиска, тем точнее структура отвечает задаче.
_TIER_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, None: 9}


def _tier_rank(tier: str | None) -> int:
    return _TIER_ORDER.get(tier, 9)


def _candidate_options(pdb_ids: list[str], ligand_comp_ids: list[str]) -> list[dict[str, Any]]:
    """Метаданные RCSB по каждому кандидату — основа для осознанного выбора."""
    wanted = {comp.upper() for comp in ligand_comp_ids}
    options: list[dict[str, Any]] = []
    for pdb_id in pdb_ids[:_MAX_CANDIDATES]:
        meta = get_structure_metadata(pdb_id)
        details = get_entity_details(pdb_id)
        cofactors = significant_cofactors(details)
        chains = sorted({ch for item in details.get("chains") or [] for ch in item.get("auth_chains") or []})
        options.append(
            {
                "pdb_id": pdb_id,
                "title": (meta.get("title") or "")[:110],
                "method": meta.get("method"),
                "resolution_A": meta.get("resolution_angstrom"),
                "chains": ", ".join(chains) or None,
                "metal_ions": ", ".join(
                    het["comp_id"] for het in cofactors if het.get("role") == "metal_ion"
                )
                or None,
                # Явные «да»/«нет» вместо True/False: булево поле модель
                # пересказывала как «структура содержит запрошенный лиганд»
                # даже тогда, когда стояло False.
                "содержит_запрошенный_лиганд": (
                    "да" if wanted & {het["comp_id"] for het in cofactors} else "нет"
                ),
                "_meta": meta,
                "_details": details,
            }
        )
    return options


# Утверждение «в структуре есть запрошенный лиганд» — то, чем модель охотнее
# всего подменяет данные: Qwen писала «выбрана структура, содержащая лиганд» для
# варианта, у которого прямо указано «нет». Обоснование сверяется с фактом.
_CLAIMS_LIGAND = re.compile(
    r"содерж\w*\s+(?:[\w,]+\s+){0,3}?(?:запрош\w*\s+|искомый\s+|этот\s+|нужный\s+)?лиганд",
    re.IGNORECASE,
)


def _checked_rationale(state: WorkflowState, rationale: str, chosen: dict[str, Any]) -> str:
    """Обоснование, не противоречащее карточке выбранной структуры."""
    if chosen.get("содержит_запрошенный_лиганд") == "да":
        return rationale
    if not _CLAIMS_LIGAND.search(rationale or ""):
        return rationale

    state.trace("structure_agent", "rationale_rejected", detail="заявлен лиганд, которого в структуре нет")
    state.limitations.append(
        f"Обоснование выбора {chosen.get('pdb_id')} от модели отклонено: в нём утверждалось "
        "наличие запрошенного лиганда, хотя в карточке структуры его нет. Выбор сохранён, "
        "обоснование заменено на фактическое."
    )
    resolution = chosen.get("resolution_A")
    return (
        f"Структура {chosen.get('pdb_id')} выбрана по разрешению "
        f"{resolution if resolution is not None else 'не указано'} Å; "
        "самого запрошенного лиганда в ней нет."
    )


def _choose(state: WorkflowState, options: list[dict[str, Any]], ligand: str | None) -> dict[str, Any]:
    """Выбор структуры: решает LLM, при её отсутствии — лучшее разрешение."""
    public = [{k: v for k, v in opt.items() if not k.startswith("_")} for opt in options]
    decision = decide(
        task=f"Выбрать кристаллическую структуру из PDB для докинга лиганда «{ligand or 'не указан'}»",
        options=public,
        criteria=(
            "структура, где содержит_запрошенный_лиганд=да, предпочтительнее; "
            "рентгеновская структура предпочтительнее; меньшее разрешение в ангстремах лучше; "
            "наличие функционального иона металла в кармане — плюс"
        ),
        context=state.request.text[:300],
    )
    if decision is not None:
        chosen = options[decision.index]
        rationale = _checked_rationale(state, decision.rationale, chosen)
        state.observations["structure_choice_source"] = decision.source
        state.observations["structure_choice_rationale"] = rationale
        decision = replace(decision, rationale=rationale)
        state.trace(
            "structure_agent",
            f"llm_choice:{chosen['pdb_id']}",
            detail=f"{decision.source}: {decision.rationale[:120]}",
        )
        state.add_evidence(
            Evidence(
                kind="hypothesis",
                value=chosen["pdb_id"],
                method="structure_choice_rationale",
                source=decision.source,
                context={
                    "rationale": decision.rationale,
                    "considered": [opt["pdb_id"] for opt in public],
                    "note": "Обоснование выбора сформулировала LLM. Сам список кандидатов получен из RCSB.",
                },
            )
        )
        return chosen

    def _rank(option: dict[str, Any]) -> tuple[int, float]:
        value = option.get("resolution_A")
        resolution = float(value) if isinstance(value, (int, float)) else 99.0
        return (0 if option.get("содержит_запрошенный_лиганд") == "да" else 1, resolution)

    chosen = min(options, key=_rank)
    state.observations["structure_choice_source"] = "heuristic"
    state.trace(
        "structure_agent",
        f"heuristic_choice:{chosen['pdb_id']}",
        detail="комплекс с этим лигандом, затем лучшее разрешение",
    )
    state.limitations.append(
        "LLM недоступна: структура выбрана эвристикой (комплекс с этим лигандом, затем "
        "лучшее разрешение), без оценки формулировки запроса."
    )
    return chosen


def _find_candidates(state: WorkflowState, ligand: str | None) -> dict[str, Any]:
    """Спуск по уровням A→D. Уровень фиксируется в наблюдениях и в evidence."""
    accession = state.observations.get("uniprot_accession")
    smiles = state.observations.get("smiles")
    inchikey = (state.observations.get("pubchem") or {}).get("inchikey") or (
        state.observations.get("pubchem_info") or {}
    ).get("inchikey")

    comp_ids = find_chem_comp_ids(smiles, inchikey=inchikey) if smiles else []
    if comp_ids:
        state.observations["ligand_comp_ids"] = comp_ids
        state.trace("structure_agent", f"chem_comp:{','.join(comp_ids)}", detail="RCSB chemical search")

    target_hint = state.observations.get("target_hint") or (
        state.observations.get("recommended_target") or {}
    ).get("name")
    text_query = " ".join(part for part in (str(ligand or ""), str(target_hint or "")) if part).strip()

    # Кандидатов accession может быть несколько (слитая мишень вроде BCR-ABL).
    # Берём тот, у которого поиск дал более высокий уровень: комплекс с этим
    # лигандом есть только у настоящей мишени лекарства.
    accessions = [a for a in (state.observations.get("uniprot_candidates") or []) if a]
    if accession and accession not in accessions:
        accessions.insert(0, accession)

    found = None
    for candidate in accessions or [None]:
        attempt = search_structures(
            uniprot_accession=candidate,
            ligand_comp_ids=comp_ids,
            text_query=text_query or None,
        )
        if found is None or _tier_rank(attempt["tier"]) < _tier_rank(found["tier"]):
            found = attempt
            found["accession"] = candidate
        if found["tier"] == "A":
            break

    if found and found.get("accession"):
        winner = found["accession"]
        if winner != accession:
            state.observations["uniprot_accession"] = winner
            state.trace(
                "structure_agent",
                f"uniprot_disambiguated:{winner}",
                detail=f"уровень {found['tier']}",
            )
        if len(accessions) > 1:
            # Факт об accession откладывался до этого момента: теперь известно,
            # какая из записей UniProt подтверждена кристаллической структурой.
            chosen_record = next(
                (r for r in (state.observations.get("uniprot_records") or []) if r.get("accession") == winner),
                state.observations.get("uniprot") or {},
            )
            state.observations["uniprot"] = chosen_record or {"accession": winner}
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=winner,
                    method="uniprot_accession",
                    source="UniProt",
                    context={
                        **(chosen_record or {}),
                        "target_name": chosen_record.get("protein_name") if chosen_record else None,
                        "disambiguated_by": f"наличие структуры уровня {found['tier']}",
                        "considered": accessions,
                    },
                )
            )
    state.observations["structure_search"] = {
        "tier": found["tier"],
        "tier_note": found["tier_note"],
        "query": found["query_description"],
        "candidates": found["pdb_ids"],
    }
    return found


def _report_composition(state: WorkflowState, pdb_id: str, details: dict[str, Any]) -> None:
    """Цепи и значимые кофакторы выбранной структуры как отдельные факты."""
    if not details.get("available"):
        state.limitations.append(
            f"Состав структуры {pdb_id} получить не удалось: цепь и кофакторы не указаны."
        )
        return

    accession = state.observations.get("uniprot_accession")
    for chain in details.get("chains") or []:
        auth = chain.get("auth_chains") or []
        if not auth:
            continue
        # Цепь интересна та, что соответствует выбранной мишени.
        if accession and chain.get("uniprot_accessions") and accession not in chain["uniprot_accessions"]:
            continue
        known = state.observations.setdefault("target_chains", [])
        # pdb_select вызывается и докингом тоже: без проверки цепи удваивались.
        known.extend(ch for ch in auth if ch not in known)
        state.add_evidence(
            Evidence(
                kind="database",
                value=", ".join(auth),
                method="structure_chain",
                source="RCSB PDB",
                context={
                    "pdb_id": pdb_id,
                    "entity_name": chain.get("name"),
                    "uniprot_accessions": chain.get("uniprot_accessions"),
                },
            )
        )

    for het in significant_cofactors(details):
        state.add_evidence(
            Evidence(
                kind="database",
                value=het["comp_id"],
                method="metal_ion" if het.get("role") == "metal_ion" else "bound_ligand",
                source="RCSB PDB",
                context={
                    "pdb_id": pdb_id,
                    "name": het.get("name"),
                    "chains": het.get("auth_chains"),
                    "role": het.get("role"),
                },
            )
        )


def _resolve_crystal_ligand(state: WorkflowState, pdb_id: str, details: dict[str, Any]) -> None:
    """Код лиганда, который реально лежит в выбранной структуре, — для контроля докинга.

    RCSB сопоставляет молекуле несколько компонентов (для ацетазоламида AZM, VFM,
    VGS…), а в конкретной структуре есть лишь один из них. Контроль по RMSD
    обязан брать именно тот, что присутствует в этой записи, иначе эталонных
    координат не найти. Результат не влияет на сам расчёт: докинг идёт вслепую.
    """
    smiles = state.observations.get("smiles")
    if not smiles:
        return
    comp_ids = state.observations.get("ligand_comp_ids")
    if not comp_ids:
        pubchem = state.observations.get("pubchem") or state.observations.get("pubchem_info") or {}
        comp_ids = find_chem_comp_ids(smiles, inchikey=pubchem.get("inchikey"))
        if comp_ids:
            state.observations["ligand_comp_ids"] = comp_ids
    present = {het.get("comp_id") for het in details.get("heteroatoms") or []}
    matching = [comp for comp in comp_ids or [] if comp in present]
    if matching:
        state.observations["crystal_comp_id"] = matching[0]
        state.trace("structure_agent", f"crystal_ligand:{matching[0]}", detail=f"в структуре {pdb_id}")
    elif comp_ids:
        state.limitations.append(
            f"В структуре {pdb_id} нет кристаллического лиганда, соответствующего запрошенной "
            "молекуле: контрольный RMSD до кристалла посчитать не с чем."
        )


def select_structure(state: WorkflowState) -> WorkflowState:
    with state.agent_step("structure_agent", "pdb_select"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else None
        )
        pdb_id = state.observations.get("pdb_id") or state.request.pdb_id
        details: dict[str, Any] | None = None

        if not pdb_id:
            found = _find_candidates(state, str(ligand) if ligand else None)
            candidate_ids = found["pdb_ids"]
            if candidate_ids:
                options = _candidate_options(candidate_ids, state.observations.get("ligand_comp_ids") or [])
                state.observations["pdb_candidates"] = [
                    {k: v for k, v in opt.items() if not k.startswith("_")} for opt in options
                ]
                chosen = _choose(state, options, str(ligand) if ligand else None)
                pdb_id = chosen["pdb_id"]
                details = chosen.get("_details")
                state.add_evidence(
                    Evidence(
                        kind="database",
                        value=found["tier"],
                        method="structure_search_tier",
                        source="RCSB PDB",
                        context={
                            "pdb_id": pdb_id,
                            "meaning": found["tier_note"],
                            "query": found["query_description"],
                            "candidates": candidate_ids,
                        },
                    )
                )
                if found["tier"] in ("C", "D"):
                    state.limitations.append(
                        f"Уровень поиска {found['tier']}: {found['tier_note']}. "
                        f"Запрос: {found['query_description']}."
                    )

        if not pdb_id:
            state.add_error("Не выбран PDB ID. Укажите структуру явно или назовите мишень.")
            return state

        pdb_id = str(pdb_id).strip().upper()
        state.observations["pdb_id"] = pdb_id
        state.request.pdb_id = pdb_id
        meta = get_structure_metadata(pdb_id)
        state.observations["pdb_meta"] = meta
        is_cached = bool(meta.get("is_cached") or meta.get("_from_cache") or meta.get("_from_download"))
        fetched_at = meta.get("fetched_at")

        state.add_evidence(
            Evidence(
                kind="database",
                value=pdb_id,
                method="pdb_structure",
                source="RCSB PDB",
                context={**meta, "is_cached": is_cached, "fetched_at": fetched_at},
            )
        )
        if meta.get("resolution_angstrom"):
            state.add_evidence(
                Evidence(
                    kind="database",
                    value=meta["resolution_angstrom"],
                    unit="Å",
                    method="resolution",
                    source="RCSB PDB",
                    context={"pdb_id": pdb_id, "is_cached": is_cached, "fetched_at": fetched_at},
                )
            )

        if details is None:
            details = get_entity_details(pdb_id)
        state.observations["pdb_entities"] = details
        _report_composition(state, pdb_id, details)
        _resolve_crystal_ligand(state, pdb_id, details)
        chains = state.observations.get("target_chains") or []
        if chains:
            # Одна цепь для рецептора: несколько копий цепи дают дубли карманов.
            state.observations["pdb_chain"] = str(chains[0])
        return state
