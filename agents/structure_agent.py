"""Агент структур: собирает кандидатов PDB и даёт LLM выбрать из них с обоснованием.

Метаданные (разрешение, метод, год, лиганды) приходят из RCSB. Модель выбирает
только номер варианта из этого списка — ни один идентификатор не берётся из её памяти.
"""
from __future__ import annotations

from typing import Any

from models.result import Evidence
from orchestration.llm_client import decide
from orchestration.state import WorkflowState
from tools.pdb_client import get_structure_metadata, search_by_ligand_target

_MAX_CANDIDATES = 5


def _candidate_options(pdb_ids: list[str]) -> list[dict[str, Any]]:
    """Метаданные RCSB по каждому кандидату — основа для осознанного выбора."""
    options: list[dict[str, Any]] = []
    for pdb_id in pdb_ids[:_MAX_CANDIDATES]:
        meta = get_structure_metadata(pdb_id)
        options.append(
            {
                "pdb_id": pdb_id,
                "title": (meta.get("title") or "")[:110],
                "method": meta.get("method"),
                "resolution_A": meta.get("resolution_angstrom"),
                "_meta": meta,
            }
        )
    return options


def _choose(state: WorkflowState, options: list[dict[str, Any]], ligand: str | None) -> dict[str, Any]:
    """Выбор структуры: решает LLM, при её отсутствии — лучшее разрешение."""
    public = [{k: v for k, v in opt.items() if not k.startswith("_")} for opt in options]
    decision = decide(
        task=f"Выбрать кристаллическую структуру из PDB для докинга лиганда «{ligand or 'не указан'}»",
        options=public,
        criteria=(
            "рентгеновская структура предпочтительнее; меньшее разрешение в ангстремах лучше; "
            "структура должна относиться к белку из запроса пользователя"
        ),
        context=state.request.text[:300],
    )
    if decision is not None:
        chosen = options[decision.index]
        state.observations["structure_choice_source"] = decision.source
        state.observations["structure_choice_rationale"] = decision.rationale
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

    def _resolution(option: dict[str, Any]) -> float:
        value = option.get("resolution_A")
        return float(value) if isinstance(value, (int, float)) else 99.0

    chosen = min(options, key=_resolution)
    state.observations["structure_choice_source"] = "heuristic"
    state.trace("structure_agent", f"heuristic_choice:{chosen['pdb_id']}", detail="лучшее разрешение")
    state.limitations.append(
        "LLM недоступна: структура выбрана эвристикой по лучшему разрешению, без оценки релевантности белка."
    )
    return chosen


def select_structure(state: WorkflowState) -> WorkflowState:
    with state.agent_step("structure_agent", "pdb_select"):
        ligand = state.observations.get("ligand_query") or (
            state.request.molecules[0] if state.request.molecules else None
        )
        pdb_id = state.observations.get("pdb_id") or state.request.pdb_id

        if not pdb_id:
            # Мишень из измерений важнее полнотекстового поиска: она опирается на данные.
            rec = state.observations.get("recommended_target") or {}
            target_hint = state.observations.get("target_hint") or rec.get("name")
            if rec.get("pdb_id"):
                pdb_id = rec["pdb_id"]
            else:
                found = search_by_ligand_target(str(ligand or ""), target_hint)
                candidate_ids = [item["pdb_id"] for item in found if item.get("pdb_id")]
                if candidate_ids:
                    options = _candidate_options(candidate_ids)
                    state.observations["pdb_candidates"] = [
                        {k: v for k, v in opt.items() if not k.startswith("_")} for opt in options
                    ]
                    chosen = _choose(state, options, str(ligand) if ligand else None)
                    pdb_id = chosen["pdb_id"]
                    state.limitations.append(
                        f"PDB {pdb_id} выбран из {len(options)} кандидатов полнотекстового поиска RCSB "
                        f"по запросу «{' '.join(part for part in (str(ligand or ''), target_hint or '') if part).strip()}». "
                        "Это не утверждение, что структура является лучшим комплексом."
                    )

        if not pdb_id:
            state.add_error("Не выбран PDB ID. Укажите структуру явно или назовите мишень.")
            return state

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
        return state
