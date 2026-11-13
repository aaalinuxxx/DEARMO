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
    assert names == [
        "rdkit_ligand",
        "activity_sources",
        "target_select",
        "activity_report",
        "pdb_select",
        "literature",
    ]
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
    assert names[:5] == [
        "rdkit_ligand",
        "activity_sources",
        "target_select",
        "activity_report",
        "pdb_select",
    ]
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


def test_structure_search_descends_the_ladder(monkeypatch):
    """Уровень поиска определяется данными: с комплексом — A, без него — B, без accession — D."""
    calls: list[dict] = []

    def fake_request(url, *, body=None, **_kwargs):
        calls.append(body or {})
        # Непустой ответ только для уровня B: комплекса с этим лигандом нет.
        description = str(body)
        if "nonpolymer_comp_id" in description:
            return {"result_set": []}
        if "nonpolymer_entity_count" in description:
            return {"result_set": [{"identifier": "1ABC"}]}
        return {"result_set": []}

    monkeypatch.setattr(pdb_client, "request_json", fake_request)
    found = pdb_client.search_structures(uniprot_accession="P00918", ligand_comp_ids=["AZM"])
    assert found["tier"] == "B"
    assert found["pdb_ids"] == ["1ABC"]
    assert "A" not in str(found["tier"])


def test_full_text_search_is_the_last_resort(monkeypatch):
    """Без accession остаётся поиск словами, и это отмечается уровнем D."""

    def fake_request(url, *, body=None, **_kwargs):
        if "full_text" in str(body):
            return {"result_set": [{"identifier": "2XYZ"}]}
        return {"result_set": []}

    monkeypatch.setattr(pdb_client, "request_json", fake_request)
    found = pdb_client.search_structures(uniprot_accession=None, text_query="carbonic anhydrase")
    assert found["tier"] == "D"
    assert found["pdb_ids"] == ["2XYZ"]


def test_structure_search_without_anything_invents_nothing(monkeypatch):
    """Ни accession, ни текста — ни одной структуры, а не случайная запись."""
    monkeypatch.setattr(pdb_client, "request_json", lambda *a, **k: {"result_set": []})
    found = pdb_client.search_structures(uniprot_accession=None, text_query=None)
    assert found["pdb_ids"] == []
    assert found["tier"] is None
