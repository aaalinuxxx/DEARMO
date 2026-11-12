"""Маршрут инструментов зависит от вопроса и не подменяет чужую мишень карбоангидразой."""
from models.request import ResearchRequest
from orchestration.state import WorkflowState
from orchestration.tool_router import choose_tool_names
from tools.bindingdb_client import _activities_from_payload, recommend_target
from tools import chembl_client, pdb_client


def _state(text: str) -> WorkflowState:
    return WorkflowState(request=ResearchRequest(text=text))


def test_property_question_does_not_open_activity_databases():
    names, source = choose_tool_names(
        _state("Посчитай logP и QED для аспирина"),
        "target_for_ligand",
    )
    assert names == ["rdkit_ligand"]
    assert source == "heuristic"


def test_target_question_uses_activity_then_structure():
    names, source = choose_tool_names(
        _state("Подбери белок для докинга аспирина и предложи PDB"),
        "target_for_ligand",
    )
    assert names == ["rdkit_ligand", "activity_sources", "pdb_select", "literature"]
    assert source == "heuristic"


def test_qwen_cannot_drop_required_sources_or_invent_tools(monkeypatch):
    monkeypatch.setattr(
        "orchestration.tool_router.ask_json",
        lambda *_args, **_kwargs: {"tools": ["literature", "vina_anywhere"]},
    )
    names, source = choose_tool_names(
        _state("Подбери белок для докинга аспирина"),
        "target_for_ligand",
    )
    assert names[:3] == ["rdkit_ligand", "activity_sources", "pdb_select"]
    assert "literature" in names
    assert "vina_anywhere" not in names
    assert source.startswith("ollama:")


def test_dock_step_stays_docking_even_if_first_scenario_was_target():
    state = _state("Подбери белок, затем выполни докинг")
    state.scenario_id = "target_for_ligand"
    names, _source = choose_tool_names(state, "dock_pockets")
    assert names == ["dock_pipeline"]


def test_bindingdb_recommendation_follows_the_measured_protein():
    payload = {
        "response": {
            "bdb.affinities": [
                {
                    "bdb.target": "Thrombin",
                    "bdb.species": "Human",
                    "bdb.affinity_type": "Ki",
                    "bdb.affinity": "12",
                }
            ]
        }
    }
    rows = _activities_from_payload(payload)
    recommended = recommend_target({"ligand": "argatroban", "activities": rows})
    assert rows[0]["target_name"] == "Thrombin"
    assert recommended["name"] == "Thrombin"
    assert recommended["pdb_id"] is None
    assert recommended["value"] == 12.0


def test_chembl_without_named_target_does_not_assume_ca2(monkeypatch):
    monkeypatch.setattr(
        chembl_client,
        "_fetch_by_molecule",
        lambda *_args, **_kwargs: [{"standard_type": "Ki", "query_mode": "molecule"}],
    )
    rows = chembl_client.fetch_activities("", "aspirin")
    assert rows[0]["query_mode"] == "molecule"


def test_local_index_wins_over_network_and_absence_falls_back(monkeypatch, tmp_path):
    """Локальный указатель структур используется, если он есть; иначе идёт поиск RCSB.

    Тест создаёт указатель сам: привязка «молекула -> PDB» не должна жить
    в непрослеживаемом файле кэша, иначе результат зависит от истории машины.
    """
    import json

    monkeypatch.setattr(pdb_client, "_CACHE_ROOT", tmp_path)
    (tmp_path / "recommendations.json").write_text(
        json.dumps({"somecompound": [{"pdb_id": "1ABC", "protein_name": "test protein"}]}),
        encoding="utf-8",
    )

    def fail_search(*_args, **_kwargs):
        raise AssertionError("при наличии локального указателя сеть не нужна")

    monkeypatch.setattr(pdb_client, "search_rcsb_entries", fail_search)
    assert pdb_client.search_by_ligand_target("somecompound")[0]["pdb_id"] == "1ABC"

    monkeypatch.setattr(
        pdb_client,
        "search_rcsb_entries",
        lambda *_args, **_kwargs: [{"pdb_id": "1OXR", "source": "rcsb_search"}],
    )
    assert pdb_client.search_by_ligand_target("unknown-compound")[0]["pdb_id"] == "1OXR"
