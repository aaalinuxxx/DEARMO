"""Сценарные тесты и проверка всех 7 обязательных проблемных случаев из ТЗ (роли.pdf)."""
from __future__ import annotations

from unittest.mock import patch

import pytest

from agents import chemistry_agent, data_agent, planner, review_agent
from models.request import ResearchRequest
from models.result import Evidence
from orchestration.state import WorkflowState
from orchestration.workflow import run_workflow
from tools import chembl_client, rdkit_tools


# ---------------------------------------------------------------------------
# Проверка валидации входных параметров ResearchRequest (Пункт 1 доски)
# ---------------------------------------------------------------------------

def test_request_validates_pdb_id_format():
    """PDB ID должен состоять ровно из 4 алфавитно-цифровых символов."""
    req = ResearchRequest(text="Докинг ацетазоламида", pdb_id="3hs4")
    assert req.pdb_id == "3HS4"

    with pytest.raises(ValueError, match="pdb_id"):
        ResearchRequest(text="Докинг", pdb_id="INVALID_PDB")

    with pytest.raises(ValueError, match="pdb_id"):
        ResearchRequest(text="Докинг", pdb_id="12")


def test_request_validates_numeric_constraints():
    """Проверка типов и диапазонов числовых ограничений constraints."""
    req = ResearchRequest(
        text="Анализ молекулы",
        constraints={"mw_min": 100, "mw_max": 450, "logp_max": 3.5},
    )
    assert req.constraints["mw_max"] == 450

    with pytest.raises(ValueError, match="должно быть числом"):
        ResearchRequest(text="Анализ", constraints={"mw_max": "five hundred"})

    with pytest.raises(ValueError, match="Некорректный диапазон"):
        ResearchRequest(text="Анализ", constraints={"mw_min": 500, "mw_max": 200})


# ---------------------------------------------------------------------------
# 7 обязательных проблемных случаев из роли.pdf (раздел 5, стр. 4)
# ---------------------------------------------------------------------------

def test_problem_case_1_invalid_smiles():
    """Случай 1: Некорректный SMILES — возвращается структурированная ошибка без падения пайплайна."""
    res = rdkit_tools.calculate_descriptors("NOT_A_VALID_SMILES_###")
    assert res.status == "error"
    assert len(res.errors) >= 1
    assert res.evidence == []

    state = WorkflowState(
        request=ResearchRequest(
            text="Рассчитай свойства для молекулы",
            molecules=["BROKEN_SMILES_((("],
        )
    )
    with patch("agents.chemistry_agent.resolve_compound", return_value={"connectivity_smiles": "BROKEN_SMILES_((("}):
        state = chemistry_agent.resolve_ligand(state)
    assert any("ошибка структуры" in e.lower() for e in state.errors)


def test_problem_case_2_ambiguous_name_requests_clarification():
    """Случай 2: Неоднозначное название соединения или мишени — возвращается status='needs_clarification'."""
    req = ResearchRequest(
        text="Подбери белок для докинга сульфонамида и оцени его активность",
        molecules=["сульфонамид"],
    )
    result, state = run_workflow(req)
    assert result.status == "needs_clarification"
    assert any("уточн" in err.lower() for err in result.errors)


def test_problem_case_3_missing_activity_data_preserved_as_none():
    """Случай 3: Отсутствие данных об активности — фиксируется явно в limitations, не подменяется на 0."""
    state = WorkflowState(
        request=ResearchRequest(
            text="Найди активность для нового соединения",
            molecules=["CCO"],
        )
    )
    with (
        patch("agents.data_agent.fetch_bindingdb_activities", return_value={"activities": []}),
        patch("agents.data_agent.fetch_chembl_activities", return_value=[]),
    ):
        state = data_agent.collect_target_data(state)

    assert state.observations.get("activity_data_missing") is True
    assert any("не найдены" in lim.lower() and "без подмены на 0" in lim.lower() for lim in state.limitations)
    activity_ev = [e for e in state.evidence if e.source in ("BindingDB", "ChEMBL")]
    assert activity_ev == []


def test_problem_case_4_incomparable_measurements_flagged():
    """Случай 4: Несопоставимые измерения (например, % ингибирования или ug/mL против Ki в нМ)."""
    raw_incomparable = [
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "assay_chembl_id": "CHEMBL1",
            "assay_type": "B",
            "standard_type": "Inhibition",
            "standard_relation": "=",
            "standard_value": "85.0",
            "standard_units": "%",
        },
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Bos taurus",
            "assay_chembl_id": "CHEMBL2",
            "assay_type": "F",
            "standard_type": "IC50",
            "standard_relation": "=",
            "standard_value": "25.0",
            "standard_units": "ug/mL",
        },
    ]
    norm = [chembl_client._normalize_record(r) for r in raw_incomparable]
    assert all(r["comparable"] is False for r in norm)
    assert norm[0]["normalized_units"] == "%"
    assert norm[1]["normalized_units"] == "ug/mL"


def test_problem_case_5_conflicting_activity_records_not_averaged():
    """Случай 5: Противоречащие друг другу записи активности — сохраняются раздельно с пометкой расхождения."""
    raw_conflicting = [
        chembl_client._normalize_record(
            {
                "molecule_chembl_id": "CHEMBL20",
                "target_chembl_id": "CHEMBL205",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "CHEMBL100",
                "assay_type": "B",
                "standard_type": "Ki",
                "standard_relation": "=",
                "standard_value": "5.0",
                "standard_units": "nM",
            }
        ),
        chembl_client._normalize_record(
            {
                "molecule_chembl_id": "CHEMBL20",
                "target_chembl_id": "CHEMBL205",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "CHEMBL100",
                "assay_type": "B",
                "standard_type": "Ki",
                "standard_relation": "=",
                "standard_value": "15.0",
                "standard_units": "uM",
            }
        ),
    ]
    deduped, lims = chembl_client.deduplicate_activities(raw_conflicting, divergence_fold=10.0)
    assert len(deduped) == 2
    values_nm = sorted(r["value_nm"] for r in deduped)
    assert values_nm == [5.0, 15000.0]
    assert all(r["conflicting_measurements"] is True for r in deduped)
    assert any("противоречащие измерения" in lim.lower() for lim in lims)


def test_problem_case_6_no_molecules_matching_constraints():
    """Случай 6: Отсутствие молекул, удовлетворяющих ограничениям — пустой список и явные причины отсева."""
    candidates = [
        "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
        "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
    ]
    report = chemistry_agent.filter_candidate_molecules(
        candidates,
        constraints={"mw_max": 120.0, "tpsa_max": 50.0},
    )
    assert report["all_rejected"] is True
    assert report["passed_molecules"] == []
    assert len(report["rejected_molecules"]) == 2
    for item in report["rejected_molecules"]:
        assert item["passed"] is False
        assert len(item["failed_reasons"]) >= 2
        assert any("mw_max" in r for r in item["failed_reasons"])


def test_problem_case_7_unsupported_conclusion_flagged_by_planner_and_review():
    """Случай 7: Просьба сделать вывод, который данные не подтверждают (клиническая эффективность / токсичность, либо число с потолка)."""
    req = ResearchRequest(
        text="Подбери белок для докинга ацетазоламида и докажи, что он безопасно вылечит пациента в клинике без побочных эффектов in vivo."
    )
    state = WorkflowState(request=req)
    state = planner.plan(state)
    assert state.observations.get("unsupported_conclusion_requested") is True
    assert any("не позволяют делать выводы о клинической эффективности" in lim.lower() for lim in state.limitations)

    state.evidence.append(Evidence(kind="calculation", value=222.25, unit="g/mol", source="RDKit", method="mw"))
    state.explanation = "Молекулярная масса равна 222.25 г/моль, а выживаемость в клинике составит 98.75%."
    state = review_agent.review(state)
    assert state.review_passed is False
    assert any("98.75" in err for err in state.review_errors)