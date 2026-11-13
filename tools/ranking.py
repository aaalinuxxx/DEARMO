"""Фильтрация, ранжирование и сводная таблица сравнения кандидатов по правилам биотехнологов."""
from __future__ import annotations

from typing import Any

from tools.chembl_client import fetch_activities as fetch_chembl_activities, normalize_to_nm
from tools.pubchem_client import resolve_compound
from tools.rdkit_tools import calculate_descriptors, evaluate_filters


def _extract_candidate_descriptors(cand: dict[str, Any], criteria: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Извлечь или рассчитать дескрипторы RDKit и проверить фильтры для одного кандидата."""
    smiles = cand.get("smiles") or cand.get("canonical_smiles") or cand.get("connectivity_smiles")
    mol_name = cand.get("molecule") or cand.get("name") or cand.get("ligand_id") or smiles or "unknown"

    if not smiles and isinstance(mol_name, str) and mol_name != "unknown":
        info = resolve_compound(mol_name)
        smiles = info.get("connectivity_smiles")

    if smiles:
        res = calculate_descriptors(smiles, constraints=criteria)
        if res.status == "error":
            return (
                {
                    "molecule": mol_name,
                    "smiles": smiles,
                    "mw": None,
                    "logp": None,
                    "tpsa": None,
                    "hbd": None,
                    "hba": None,
                    "rotatable_bonds": None,
                    "qed": None,
                    "filter_passed": False,
                },
                res.errors or ["Некорректный SMILES"],
            )
        vals = {e.method: e.value for e in res.evidence if e.method}
        canonical = res.evidence[0].context.get("canonical_smiles", smiles) if res.evidence else smiles
        filt = (
            res.evidence[0].context.get("filter_evaluation")
            if res.evidence and "filter_evaluation" in res.evidence[0].context
            else evaluate_filters(vals, constraints=criteria)
        )
        return (
            {
                "molecule": mol_name,
                "smiles": canonical,
                "mw": round(float(vals["molecular_weight"]), 2) if vals.get("molecular_weight") is not None else None,
                "logp": round(float(vals["logp"]), 3) if vals.get("logp") is not None else None,
                "tpsa": round(float(vals["tpsa"]), 2) if vals.get("tpsa") is not None else None,
                "hbd": int(vals["hbd"]) if vals.get("hbd") is not None else None,
                "hba": int(vals["hba"]) if vals.get("hba") is not None else None,
                "rotatable_bonds": int(vals["rotatable_bonds"]) if vals.get("rotatable_bonds") is not None else None,
                "qed": round(float(vals["qed"]), 4) if vals.get("qed") is not None else None,
                "filter_passed": bool(filt.get("passed", True)),
            },
            list(filt.get("failed_reasons", [])),
        )

    # Если SMILES не передан, но переданы готовые численные дескрипторы в словаре кандидата
    desc_map = {
        "molecular_weight": cand.get("mw") if cand.get("mw") is not None else cand.get("molecular_weight", 0.0),
        "logp": cand.get("logp", 0.0),
        "tpsa": cand.get("tpsa", 0.0),
        "hbd": cand.get("hbd", 0),
        "hba": cand.get("hba", 0),
        "rotatable_bonds": cand.get("rotatable_bonds", 0),
        "qed": cand.get("qed", 0.5),
    }
    filt = evaluate_filters(desc_map, constraints=criteria)
    return (
        {
            "molecule": mol_name,
            "smiles": None,
            "mw": cand.get("mw", cand.get("molecular_weight")),
            "logp": cand.get("logp"),
            "tpsa": cand.get("tpsa"),
            "hbd": cand.get("hbd"),
            "hba": cand.get("hba"),
            "rotatable_bonds": cand.get("rotatable_bonds"),
            "qed": cand.get("qed"),
            "filter_passed": bool(filt.get("passed", True)),
        },
        list(filt.get("failed_reasons", [])),
    )


def rank_candidates(
    candidates: list[dict[str, Any]],
    criteria: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Отфильтровать и проранжировать кандидатов с явным указанием причин исключения и источников аргументов.

    Правила:
    1. Пропуски в измерениях активности (`ki_nm=None`) сохраняются как `None` (никогда не заменяются на 0).
    2. Несопоставимые единицы (`%`, `ug/mL`) или разные организмы/методики помечаются в `comparability_warnings`
       и не участвуют в прямом численном ранжировании по Ki.
    3. Кандидаты со статусом `filter_status == "passed"` всегда идут выше `"rejected"`.
    4. Среди прошедших фильтры приоритет имеют соединения с подтверждённой сопоставимой активностью (меньший Ki/IC50 в нМ),
       затем лучшей оценкой докинга (`vina_score`, более отрицательный), затем большим `qed`.
    """
    crit = dict(criteria or {})
    target_organism_req = str(crit.get("target_organism") or "Homo sapiens").lower()
    ki_max_nm = crit.get("ki_max_nm")
    ic50_max_nm = crit.get("ic50_max_nm")

    evaluated: list[dict[str, Any]] = []

    for cand in candidates:
        base_info, failed_reasons = _extract_candidate_descriptors(cand, crit)
        comparability_warnings: list[str] = []
        evidence_sources: list[str] = []
        missing_data: list[str] = []

        if base_info.get("mw") is not None:
            evidence_sources.append("calculation:RDKit(MW,logP,TPSA,QED)")
        else:
            missing_data.append("descriptors")

        # Обработка экспериментальной активности (Ki / IC50)
        raw_ki = cand.get("ki_nm") if "ki_nm" in cand else cand.get("ki")
        raw_unit = cand.get("unit") or cand.get("units") or ("nM" if "ki_nm" in cand else None)
        comparable_flag = cand.get("comparable", True)

        ki_nm: float | None = None
        if raw_ki is not None:
            if raw_unit and raw_unit != "nM":
                norm_val, _, is_comp = normalize_to_nm(raw_ki, raw_unit)
                if is_comp:
                    ki_nm = norm_val
                else:
                    comparable_flag = False
                    comparability_warnings.append(
                        f"Несопоставимая единица измерения ({raw_ki} {raw_unit}) не переводится в nM"
                    )
            else:
                try:
                    ki_nm = float(raw_ki)
                except (ValueError, TypeError):
                    ki_nm = None

        organism = cand.get("target_organism") or cand.get("organism")
        if organism and target_organism_req and target_organism_req not in str(organism).lower():
            comparable_flag = False
            comparability_warnings.append(
                f"Измерение получено на другом организме ({organism} вместо {crit.get('target_organism', 'Homo sapiens')})"
            )

        if cand.get("conflicting_measurements"):
            comparability_warnings.append("Обнаружены противоречащие друг другу измерения активности")

        if not comparable_flag and ki_nm is not None:
            # Не используем несопоставимое измерение как валидный нМ-ранг
            ki_nm_for_sort = float("inf")
        elif ki_nm is not None:
            ki_nm_for_sort = ki_nm
            evidence_sources.append("database:BindingDB/ChEMBL(Ki_nM)")
        else:
            ki_nm_for_sort = float("inf")
            missing_data.append("experimental_Ki")

        # Проверка порогов активности из criteria
        if ki_max_nm is not None and ki_nm is not None and comparable_flag:
            if ki_nm > float(ki_max_nm):
                failed_reasons.append(f"Constraint: Ki ({ki_nm:.2f} nM) > ki_max_nm ({ki_max_nm} nM)")

        ic50_nm = cand.get("ic50_nm")
        if ic50_max_nm is not None and ic50_nm is not None and comparable_flag:
            if float(ic50_nm) > float(ic50_max_nm):
                failed_reasons.append(f"Constraint: IC50 ({ float(ic50_nm):.2f} nM) > ic50_max_nm ({ic50_max_nm} nM)")

        vina_score = cand.get("vina_score")
        if vina_score is not None:
            try:
                vina_val = float(vina_score)
                evidence_sources.append("prediction:AutoDock_Vina")
            except (ValueError, TypeError):
                vina_val = None
        else:
            vina_val = None
            missing_data.append("docking_pose")

        passed = len(failed_reasons) == 0 and bool(base_info.get("filter_passed", True))
        filter_status = "passed" if passed else "rejected"

        qed_val = base_info.get("qed")
        qed_sort = -float(qed_val) if qed_val is not None else 0.0
        vina_sort = vina_val if vina_val is not None else 0.0
        has_comparable_ki = 0 if (ki_nm is not None and comparable_flag) else 1

        row = {
            "molecule": base_info["molecule"],
            "smiles": base_info["smiles"],
            "mw": base_info["mw"],
            "logp": base_info["logp"],
            "tpsa": base_info["tpsa"],
            "hbd": base_info["hbd"],
            "hba": base_info["hba"],
            "rotatable_bonds": base_info["rotatable_bonds"],
            "qed": qed_val,
            "ki_nm": ki_nm if comparable_flag else None,
            "raw_activity": raw_ki,
            "raw_unit": raw_unit,
            "comparable": bool(comparable_flag and raw_ki is not None),
            "vina_score": vina_val,
            "filter_status": filter_status,
            "passed": passed,
            "failed_reasons": failed_reasons,
            "comparability_warnings": comparability_warnings,
            "evidence_sources": evidence_sources,
            "missing_data": missing_data,
            "preliminary_choice": ki_nm is None or not comparable_flag,
            "_sort_key": (
                0 if passed else 1,
                has_comparable_ki,
                ki_nm_for_sort,
                vina_sort,
                qed_sort,
                str(base_info["molecule"]),
            ),
        }
        evaluated.append(row)

    evaluated.sort(key=lambda r: r["_sort_key"])
    for idx, item in enumerate(evaluated, start=1):
        item.pop("_sort_key", None)
        item["rank"] = idx

    return evaluated


def build_candidates_comparison_table(
    molecules: list[str],
    *,
    target_id: str = "",
    constraints: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Собрать данные из PubChem, RDKit и ChEMBL для списка молекул и построить сводную таблицу сравнения."""
    raw_candidates: list[dict[str, Any]] = []
    for mol_query in molecules:
        info = resolve_compound(mol_query)
        smiles = info.get("connectivity_smiles") or (mol_query if not info.get("error") else None)
        acts = fetch_chembl_activities(target_id=target_id, molecule_id=mol_query)

        best_ki: float | None = None
        best_ic50: float | None = None
        comp_flag = False
        conflict_flag = False
        organism_val = None

        for act in acts:
            if act.get("conflicting_measurements"):
                conflict_flag = True
            if act.get("comparable") and act.get("value_nm") is not None:
                comp_flag = True
                organism_val = act.get("target_organism")
                stype = str(act.get("standard_type") or "").upper()
                val = float(act["value_nm"])
                if stype == "KI" and (best_ki is None or val < best_ki):
                    best_ki = val
                elif stype == "IC50" and (best_ic50 is None or val < best_ic50):
                    best_ic50 = val

        raw_candidates.append(
            {
                "molecule": info.get("title") or mol_query,
                "query": mol_query,
                "smiles": smiles or mol_query,
                "ki_nm": best_ki if best_ki is not None else best_ic50,
                "ic50_nm": best_ic50,
                "unit": "nM" if (best_ki is not None or best_ic50 is not None) else None,
                "comparable": comp_flag,
                "conflicting_measurements": conflict_flag,
                "target_organism": organism_val or "Homo sapiens",
            }
        )

    return rank_candidates(raw_candidates, criteria=constraints)
