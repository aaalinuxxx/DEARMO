"""Агент мишени: какой белок рекомендовать и как он называется в UniProt.

Разделение ответственности с `data_agent` намеренное. Данные собирает один
агент, решение принимает другой: так на схеме видно, где в системе находится
точка выбора, и её легко проверить отдельно от загрузки измерений.

Решение принимает LLM, но только как индекс в списке, который собрали
измерения BindingDB. Accession модель не видит вовсе — его резолвит UniProt.
"""
from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from models.result import Evidence
from orchestration.llm_client import decide
from orchestration.state import WorkflowState
from tools.pdb_client import count_complexes, find_chem_comp_ids
from tools.uniprot_client import resolve_target_parts

_MAX_OPTIONS = 8

# Организмы, которые BindingDB называет человеческими.
_HUMAN = ("human", "homo sapiens")


def _is_human(organism: str | None) -> bool:
    lowered = (organism or "").strip().lower()
    return any(marker in lowered for marker in _HUMAN)


# Хвост «[A65S,N67Q,…]» или «[19-324]» в названии мишени BindingDB — это мутант
# или фрагмент белка, а не другой белок. Без отделения такого хвоста карбоангидраза II
# в выдаче ацетазоламида числилась «белком с мутациями» и не сопоставлялась с UniProt.
_VARIANT_TAIL = re.compile(r"\s*\[([^\]]{1,120})\]\s*$")


def split_variant(name: str) -> tuple[str, str | None]:
    """«Carbonic anhydrase 2 [A65S,N67Q]» -> («Carbonic anhydrase 2», «A65S,N67Q»)."""
    match = _VARIANT_TAIL.search(name or "")
    if not match:
        return (name or "").strip(), None
    return name[: match.start()].strip(), match.group(1).strip()


def group_targets(activities: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Сгруппировать измерения по мишени: кандидаты для осознанного выбора.

    Числа здесь — реальные значения из BindingDB; LLM их не меняет, она только
    указывает, какая из мишеней отвечает на вопрос пользователя.
    """
    groups: dict[tuple[str, str], dict[str, Any]] = {}
    for row in activities:
        raw_name = (row.get("target_name") or "").strip()
        if not raw_name:
            continue
        name, variant = split_variant(raw_name)
        organism = (row.get("organism") or "").strip()
        key = (name.lower(), organism.lower())
        try:
            value = float(row["value"])
        except (TypeError, ValueError, KeyError):
            value = None
        item = groups.setdefault(
            key,
            {
                "target_name": name,
                "organism": organism,
                "measurements": 0,
                "best_affinity_nM": None,
                "affinity_type": row.get("affinity_type"),
                "link": row.get("link"),
                "variants": [],
            },
        )
        if variant and variant not in item["variants"]:
            item["variants"].append(variant)
        item["measurements"] += 1
        if value is not None and row.get("relation") in (None, "=", ""):
            current = item["best_affinity_nM"]
            if current is None or value < current:
                item["best_affinity_nM"] = value
                item["affinity_type"] = row.get("affinity_type")
                item["link"] = row.get("link") or item.get("link")
    # Сортировка сначала по организму, потом по силе связывания.
    #
    # Если ранжировать только по аффинности, список обрезается на ортологах из
    # бактерий, дрожжей и мыши, а человеческая мишень до модели не доходит. Для
    # ацетазоламида так терялась карбоангидраза II: в выдаче BindingDB она идёт
    # девятой. Все организмы остаются в списке, меняется только порядок.
    ranked = sorted(
        groups.values(),
        key=lambda item: (
            0 if _is_human(item["organism"]) else 1,
            item["best_affinity_nM"] if item["best_affinity_nM"] is not None else float("inf"),
            -item["measurements"],
        ),
    )
    return ranked[:_MAX_OPTIONS]


def _ligand_comp_ids(state: WorkflowState) -> list[str]:
    comp_ids = state.observations.get("ligand_comp_ids")
    if comp_ids:
        return list(comp_ids)
    smiles = state.observations.get("smiles")
    if not smiles:
        return []
    pubchem = state.observations.get("pubchem") or state.observations.get("pubchem_info") or {}
    comp_ids = find_chem_comp_ids(smiles, inchikey=pubchem.get("inchikey"))
    if comp_ids:
        state.observations["ligand_comp_ids"] = comp_ids
    return list(comp_ids)


def _structural_evidence(option: dict[str, Any], comp_ids: list[str]) -> dict[str, Any]:
    """UniProt и число кристаллических комплексов с лигандом для одного кандидата."""
    records = resolve_target_parts(option["target_name"], option.get("organism"))
    counts = [count_complexes(r["accession"], comp_ids) for r in records] if comp_ids else []
    known = [c for c in counts if c is not None]
    return {"records": records, "complexes": max(known) if known else None}


def enrich_with_structures(state: WorkflowState, options: list[dict[str, Any]]) -> None:
    """Дополнить кандидатов структурным свидетельством и переупорядочить.

    BindingDB REST отдаёт неполную выдачу: для ацетазоламида в ней нет дикого
    типа карбоангидразы II, только мутант, и по силе связывания побеждает другая
    изоформа. Число кристаллических комплексов белка именно с этим лигандом
    берётся из PDB, то есть из независимого эксперимента, и показывает, какой
    белок на самом деле изучали вместе с этой молекулой.
    """
    comp_ids = _ligand_comp_ids(state)
    if not comp_ids:
        return
    with ThreadPoolExecutor(max_workers=4) as pool:
        found = list(pool.map(lambda opt: _structural_evidence(opt, comp_ids), options))
    for option, info in zip(options, found):
        records = info["records"]
        option["_uniprot_records"] = records
        if records:
            option["uniprot"] = "/".join(r["accession"] for r in records)
        if info["complexes"] is not None:
            option["комплексов_в_PDB_с_этим_лигандом"] = info["complexes"]

    options.sort(
        key=lambda item: (
            0 if _is_human(item["organism"]) else 1,
            -(item.get("комплексов_в_PDB_с_этим_лигандом") or 0),
            item["best_affinity_nM"] if item["best_affinity_nM"] is not None else float("inf"),
            -item["measurements"],
        )
    )


def _note_variants(state: WorkflowState, chosen: dict[str, Any]) -> None:
    if chosen.get("variants"):
        shown = "; ".join(chosen["variants"][:3])
        state.limitations.append(
            f"Измерения BindingDB по «{chosen['target_name']}» относятся к мутантным вариантам или "
            f"фрагментам ({shown}), дикого типа в выдаче нет: значения аффинности не переносятся на "
            "природный белок без оговорки."
        )


# Утверждение «у выбранной мишени самое сильное связывание» модель делает охотно,
# даже когда это не так: для ацетазоламида Qwen написала «наименьшая константа (25 нМ)»
# при том, что у другого кандидата было 3.1 нМ. Число настоящее, неверно сравнение.
_CLAIMS_STRONGEST = re.compile(
    r"(наименьш\w+|минимальн\w+|наилучш\w+|самы\w+\s+силь\w+|наиболее\s+силь\w+|сильнейш\w+|"
    r"лучш\w+)\s+(?:\w+\s+){0,2}?(констант\w*|связывани\w*|аффинност\w*|Ki|Kd|IC50)",
    re.IGNORECASE,
)


def _checked_rationale(
    state: WorkflowState, rationale: str, chosen: dict[str, Any], options: list[dict[str, Any]]
) -> str:
    """Обоснование без неверных сравнений: утверждения о «сильнейшем связывании» сверяются с данными."""
    if not _CLAIMS_STRONGEST.search(rationale or ""):
        return rationale
    values = [o["best_affinity_nM"] for o in options if o.get("best_affinity_nM") is not None]
    mine = chosen.get("best_affinity_nM")
    if not values or mine is None or mine <= min(values):
        return rationale

    strongest = min(
        (o for o in options if o.get("best_affinity_nM") is not None),
        key=lambda o: o["best_affinity_nM"],
    )
    state.trace("target_agent", "rationale_rejected", detail="заявлено сильнейшее связывание, которого нет")
    state.limitations.append(
        f"Обоснование выбора «{chosen['target_name']}» от модели отклонено: в нём утверждалось "
        "сильнейшее связывание, но среди кандидатов оно у другой мишени. Выбор сохранён, "
        "обоснование заменено на фактическое."
    )
    reasons = []
    if _is_human(chosen.get("organism")):
        reasons.append("человеческий белок")
    complexes = chosen.get("комплексов_в_PDB_с_этим_лигандом")
    if complexes:
        reasons.append(f"{complexes} кристаллических комплексов с этим лигандом в PDB")
    basis = ", ".join(reasons) or "совокупность данных"
    return (
        f"«{chosen['target_name']}» выбрана по критериям: {basis}. Самое сильное связывание в выдаче "
        f"у «{strongest['target_name']}» ({strongest['best_affinity_nM']} нМ), но оно не единственный довод."
    )


def choose_target(state: WorkflowState, activities: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Какую мишень считать ответом: решает LLM по сгруппированным измерениям."""
    options = group_targets(activities)
    if not options:
        return None
    enrich_with_structures(state, options)
    state.observations["target_options"] = [
        {k: v for k, v in opt.items() if not k.startswith("_")} for opt in options
    ]

    public = []
    for opt in options:
        card = {k: v for k, v in opt.items() if k not in ("link", "variants") and not k.startswith("_")}
        if opt.get("variants"):
            # Мутант или фрагмент: об этом модель должна знать, а не читать скобки в названии.
            card["измерения_относятся_к"] = "мутантам или фрагментам белка, не к дикому типу"
        public.append(card)
    decision = decide(
        task="Выбрать белок-мишень, по которому отвечать пользователю",
        options=public,
        criteria=(
            "соответствие вопросу пользователя; человеческий белок предпочтительнее, "
            "если пользователь не спросил про другой организм; больше комплексов_в_PDB_с_этим_лигандом "
            "означает, что именно этот белок изучали с данной молекулой; сильнее связывание "
            "(меньше наномолей); больше независимых измерений"
        ),
        context=state.request.text[:300],
    )
    if decision is None:
        state.observations["target_choice_source"] = "heuristic"
        if len(options) > 1:
            state.limitations.append(
                "LLM недоступна: мишень выбрана эвристикой (человек, затем число кристаллических "
                f"комплексов с лигандом, затем сила связывания) из {len(options)} кандидатов, "
                "без учёта формулировки запроса."
            )
        _note_variants(state, options[0])
        return options[0]

    chosen = options[decision.index]
    _note_variants(state, chosen)
    rationale = _checked_rationale(state, decision.rationale, chosen, options)
    state.observations["target_choice_source"] = decision.source
    state.observations["target_choice_rationale"] = rationale
    state.trace(
        "target_agent",
        f"llm_target:{chosen['target_name'][:40]}",
        detail=f"{decision.source}: {rationale[:120]}",
    )
    state.add_evidence(
        Evidence(
            kind="hypothesis",
            value=chosen["target_name"],
            method="target_choice_rationale",
            source=decision.source,
            context={
                "rationale": rationale,
                "considered": [f"{o['target_name']} ({o['organism']})" for o in options],
                "note": "Выбор мишени обосновала LLM; измерения взяты из BindingDB.",
            },
        )
    )
    return chosen


def resolve_identity(state: WorkflowState, target: dict[str, Any]) -> dict[str, Any] | None:
    """Сопоставить выбранную мишень записи UniProt и зафиксировать accession."""
    records = target.get("_uniprot_records") or resolve_target_parts(
        target.get("target_name") or "", target.get("organism")
    )
    record = records[0] if records else None
    if len(records) > 1:
        # Слитая мишень: какая половина является мишенью лекарства, решит поиск
        # структур по комплексу с этим лигандом, а не порядок слов в названии.
        state.observations["uniprot_candidates"] = [r["accession"] for r in records]
        state.observations["uniprot_records"] = records
        state.limitations.append(
            f"Название мишени «{target.get('target_name')}» соответствует нескольким записям "
            f"UniProt ({', '.join(r['accession'] + ' ' + str(r.get('gene')) for r in records)}). "
            "Выбор между ними сделан по наличию кристаллического комплекса с этим лигандом."
        )
    if record is None:
        state.observations["uniprot_unresolved"] = target.get("target_name")
        state.limitations.append(
            f"Мишень «{target.get('target_name')}» не сопоставлена записи UniProt: "
            "accession не указан. Поиск структуры пойдёт по названию, что менее точно. "
            "Идентификатор не подбирается по памяти модели."
        )
        return None

    state.observations["uniprot"] = record
    state.observations["uniprot_accession"] = record["accession"]
    if len(records) > 1:
        # Какая из записей настоящая, решит поиск структур. Факт об accession
        # запишет он же: иначе в ответ попадёт отвергнутая половина слияния.
        state.trace("target_agent", "uniprot_pending", detail=", ".join(r["accession"] for r in records))
        return record
    state.trace("target_agent", f"uniprot:{record['accession']}", detail=record.get("name_match"))
    state.add_evidence(
        Evidence(
            kind="database",
            value=record["accession"],
            method="uniprot_accession",
            source="UniProt",
            context={
                **record,
                "target_name": record.get("protein_name") or target.get("target_name"),
                "organism": record.get("organism"),
                "queried_name": target.get("target_name"),
            },
        )
    )
    if record.get("name_match") == "partial":
        state.limitations.append(
            f"Название мишени из BindingDB («{target.get('target_name')}») совпало с записью "
            f"UniProt {record['accession']} («{record.get('protein_name')}») лишь частично."
        )
    if not record.get("reviewed"):
        state.limitations.append(
            f"Запись UniProt {record['accession']} не проверена вручную (TrEMBL, не Swiss-Prot)."
        )
    return record


def _matches(target_name: str, row_name: str) -> bool:
    left = target_name.strip().lower()
    right = (row_name or "").strip().lower()
    return bool(left and right and (left in right or right in left))


def _honour_named_target(
    state: WorkflowState, activities: list[dict[str, Any]]
) -> dict[str, Any] | None:
    """Учесть уже названную мишень — и отбросить её, если её выдумала модель.

    Мишень пользователя авторитетна: он спросил именно про неё, даже если
    измерений нет. А вот подсказка планировщика — это извлечение из текста
    языковой моделью, и она может оказаться белком, которого в запросе не было
    (наблюдали `ACO1_HUMAN` в запросе про ацетазоламид). Поэтому подсказку
    модели принимаем только тогда, когда по ней есть измерения; иначе выбор
    делается по данным, а отброшенная подсказка попадает в ограничения.
    """
    user_target = str(state.request.target_id or "").strip()
    if user_target:
        state.observations["target_choice_source"] = "user"
        return {"target_name": user_target, "organism": _organism_for(activities, user_target)}

    hint = str(state.observations.get("target_hint") or "").strip()
    if not hint:
        return None

    source = str(state.observations.get("target_hint_source") or "")
    if not source.startswith("ollama"):
        # Подсказка получена разбором текста по правилам — она дословно из запроса.
        state.observations["target_choice_source"] = "planner_keywords"
        return {"target_name": hint, "organism": _organism_for(activities, hint)}

    # Подсказка модели принимается как слово пользователя только если эти слова
    # есть в самом запросе. Иначе это знание модели, а не просьба: на запрос
    # «подбери мишень для ацетазоламида» Qwen писала target_hint «Carbonic anhydrase»,
    # под которое подходят и бактериальные ортологи, и система отвечала про
    # Vibrio cholerae, вообще не запуская выбор по данным.
    if hint.lower() not in state.request.text.lower():
        state.observations["ignored_target_hint"] = hint
        state.observations.pop("target_hint", None)
        state.trace("target_agent", f"ignore_hint:{hint[:40]}", detail="мишень не названа в запросе")
        state.limitations.append(
            f"Подсказка планировщика «{hint}» не принята: в запросе такой мишени нет. "
            "Мишень выбрана по данным измерений и структур."
        )
        return None

    matched = [row for row in activities if _matches(hint, row.get("target_name") or "")]
    if matched:
        # Совпадение по подстроке даёт и «Carbonic anhydrase» дрожжей на запрос
        # про «Carbonic anhydrase 2» человека. Организм выбирается по тем же
        # правилам, что и мишень: точное совпадение имени, затем человек.
        exact = [
            row
            for row in matched
            if (row.get("target_name") or "").strip().lower() == hint.strip().lower()
        ]
        best = group_targets(exact or matched)
        state.observations["target_choice_source"] = "planner_llm_confirmed"
        return {
            "target_name": best[0]["target_name"] if best else hint,
            "organism": best[0]["organism"] if best else matched[0].get("organism"),
        }

    if not activities:
        # Проверить нечем: измерений нет вовсе. Берём подсказку, но честно говорим об этом.
        state.observations["target_choice_source"] = "planner_llm_unverified"
        state.limitations.append(
            f"Мишень «{hint}» извлечена из запроса языковой моделью и не подтверждена "
            "измерениями: для этого соединения в BindingDB записей нет."
        )
        return {"target_name": hint, "organism": None}

    state.observations["rejected_target_hint"] = hint
    state.observations.pop("target_hint", None)
    state.trace("target_agent", f"reject_hint:{hint[:40]}", detail="нет измерений по этой мишени")
    state.limitations.append(
        f"Подсказка планировщика «{hint}» отброшена: ни одно из {len(activities)} измерений "
        "BindingDB по этому соединению не относится к такой мишени. Мишень выбрана по данным."
    )
    return None


def select_target(state: WorkflowState) -> WorkflowState:
    """Выбрать мишень по измерениям и дать ей проверяемый идентификатор."""
    with state.agent_step("target_agent", "select_target"):
        activities = (state.observations.get("bindingdb") or {}).get("activities") or []
        chosen = _honour_named_target(state, activities)
        if chosen is None:
            chosen = choose_target(state, activities)

        if not chosen:
            if state.observations.get("bindingdb"):
                state.limitations.append(
                    "Мишень не выбрана: измерений для этого соединения не нашлось, "
                    "а выдумывать белок система не станет."
                )
            return state

        state.observations["target_hint"] = chosen["target_name"]
        state.observations["chosen_target"] = chosen
        resolve_identity(state, chosen)
        return state


def _organism_for(activities: list[dict[str, Any]], target_name: str) -> str | None:
    """Организм для названной пользователем мишени — из измерений, если они есть."""
    wanted = target_name.strip().lower()
    for row in activities:
        name = (row.get("target_name") or "").lower()
        if wanted and (wanted in name or name in wanted):
            return row.get("organism")
    return None
