"""Маршрутизация workflow без GigaChat."""
from models.request import ResearchRequest
from orchestration.scenarios import classify_by_keywords
from orchestration.workflow import run_workflow


def test_keyword_classifier_targets():
    assert classify_by_keywords("Подбери белок для ацетазоламида") == "target_for_ligand"


def test_scenario1_target_acetazolamide():
    text = (
        "Подбери белок для докинга ацетазоламида. "
        "Обоснуй выбор экспериментальными данными и предложи структуру из PDB."
    )
    result, state = run_workflow(ResearchRequest(text=text))
    assert result.status == "success"
    assert state.scenario_id == "target_for_ligand"
    assert any(e.method == "pubchem_cid" for e in result.evidence)
    assert any(e.value == "3HS4" for e in result.evidence if e.method == "pdb_structure")


def test_scenario3_compare_pockets():
    text = "Сравни три лиганда в двух карманах. Покажи, меняется ли порядок по Vina."
    result, state = run_workflow(ResearchRequest(text=text))
    assert result.status == "success"
    assert state.scenario_id == "compare_docking"
    assert state.artifacts.get("heatmap_png")
    assert any(e.method == "ligand_ranking_by_pocket" for e in result.evidence)


def test_compound_request_runs_target_then_dock_then_pose():
    text = (
        "Подбери белок для ацетазоламида, затем выполни докинг в 3HS4 "
        "и проверь правдоподобность позы."
    )
    result, state = run_workflow(ResearchRequest(text=text))
    assert result.status == "success"
    assert state.observations["completed_steps"] == [
        "target_for_ligand",
        "dock_pockets",
        "validate_pose",
    ]


def test_scenario4_validate_pose():
    text = "Проверь правдоподобность позы ацетазоламида в 3HS4."
    result, _state = run_workflow(ResearchRequest(text=text))
    assert result.status == "success"
    assert any(e.method == "checks_passed_percent" for e in result.evidence)
