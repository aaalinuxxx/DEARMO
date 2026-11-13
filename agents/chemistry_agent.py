"""Агент химии: разрешение структуры, очистка солей, RDKit-дескрипторы (включая QED) и фильтрация по критериям."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from models.docking import load_poses_csv
from models.result import Evidence
from orchestration.state import WorkflowState
from tools.pubchem_client import resolve_compound
from tools.rdkit_tools import calculate_descriptors, evaluate_filters


def resolve_ligand(state: WorkflowState) -> WorkflowState:
    """Разрешить структуру лиганда, проверить неоднозначность и рассчитать свойства RDKit."""
    with state.agent_step("chemistry_agent", "resolve_ligand"):
        ligand = state.observations.get("ligand_query") or state.observations.get("ligand_id")
        if not ligand and state.request.molecules:
            ligand = state.request.molecules[0]
        if not ligand:
            return state

        state.observations["ligand_query"] = ligand
        info = resolve_compound(ligand)
        state.observations["pubchem_info"] = info

        if info.get("ambiguous"):
            candidates = info.get("candidates") or []
            state.observations["clarification_candidates"] = candidates
            cand_titles = ", ".join(
                f"{c.get('title')} (CID {c.get('cid')})" for c in candidates[:5] if c.get("cid")
            )
            state.add_error(
                f"Уточните соединение: запрос '{ligand}' неоднозначен. Возможные варианты: {cand_titles or 'несколько структур'}."
            )
            return state

        smiles = info.get("connectivity_smiles") or state.observations.get("smiles")
        if not smiles:
            # Если переданная строка похожа на битый SMILES или неизвестный идентификатор
            desc_attempt = calculate_descriptors(ligand, constraints=state.request.constraints)
            if desc_attempt.status == "error":
                state.limitations.append(f"Не удалось получить SMILES для '{ligand}' в PubChem.")
                return state
            smiles = ligand

        state.observations["smiles"] = smiles
        desc = calculate_descriptors(smiles, constraints=state.request.constraints)
        if desc.status == "error":
            for err in desc.errors:
                state.add_error(f"Ошибка структуры '{smiles}': {err}")
            return state

        for item in desc.evidence:
            state.add_evidence(item)
            if "filter_evaluation" in item.context:
                state.observations["filter_evaluation"] = item.context["filter_evaluation"]
        for lim in desc.limitations:
            if lim not in state.limitations:
                state.limitations.append(lim)
        return state


def filter_candidate_molecules(
    molecules: list[str],
    constraints: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Проверить список молекул по правилам биотехнологов (Липински, Вебера, PAINS) и пользовательским ограничениям.

    Возвращает словарь с прошедшими молекулами (`passed_molecules`), отсеянными (`rejected_molecules`)
    и явными причинами исключения для каждого кандидата.
    """
    passed: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []

    for mol_query in molecules:
        info = resolve_compound(mol_query)
        smiles = info.get("connectivity_smiles") or mol_query
        res = calculate_descriptors(smiles, constraints=constraints)
        if res.status != "success":
            rejected.append(
                {
                    "query": mol_query,
                    "smiles": smiles,
                    "passed": False,
                    "failed_reasons": res.errors or ["Некорректный SMILES"],
                }
            )
            continue

        desc_map = {item.method: item.value for item in res.evidence if item.method}
        filt = (
            res.evidence[0].context.get("filter_evaluation")
            if res.evidence and "filter_evaluation" in res.evidence[0].context
            else evaluate_filters(desc_map, constraints=constraints)
        )
        entry = {
            "query": mol_query,
            "smiles": res.evidence[0].context.get("canonical_smiles", smiles) if res.evidence else smiles,
            "descriptors": desc_map,
            "passed": bool(filt.get("passed")),
            "failed_reasons": list(filt.get("failed_reasons", [])),
        }
        if entry["passed"]:
            passed.append(entry)
        else:
            rejected.append(entry)

    return {
        "passed_molecules": passed,
        "rejected_molecules": rejected,
        "all_rejected": len(passed) == 0 and len(rejected) > 0,
    }


def append_ligand_descriptors_for_poses(state: WorkflowState, csv_path: Path) -> WorkflowState:
    """Добавить дескрипторы RDKit (MW, logP, QED) для лигандов из таблицы поз докинга."""
    with state.agent_step("chemistry_agent", "ligand_descriptors"):
        rows = load_poses_csv(csv_path)
        ligands = sorted({r.ligand_id for r in rows})
        for lig in ligands:
            # Структуру лиганда всегда ищем в PubChem: локального словаря SMILES нет.
            info = resolve_compound(lig)
            smi = info.get("connectivity_smiles")
            if not smi:
                state.limitations.append(f"Структура лиганда '{lig}' не найдена в PubChem, дескрипторы не считались.")
                continue
            desc = calculate_descriptors(smi)
            for item in desc.evidence:
                if item.method in ("molecular_weight", "logp", "qed"):
                    tagged = Evidence(
                        kind=item.kind,
                        value=item.value,
                        unit=item.unit,
                        source=item.source,
                        method=f"{item.method}_{lig}",
                        context={**item.context, "ligand_id": lig},
                    )
                    state.add_evidence(tagged)
        state.limitations.append("Молекулярная масса, logP и QED относятся к лиганду и не меняются при смене кармана.")
        return state


def analyze(state: WorkflowState) -> WorkflowState:
    return resolve_ligand(state)
