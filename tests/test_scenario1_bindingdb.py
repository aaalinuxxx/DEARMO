"""Рекомендация мишени BindingDB опирается только на измерения этого соединения."""
from tools.bindingdb_client import fetch_activities, recommend_target


def test_bindingdb_cache_is_keyed_by_molecule():
    data = fetch_activities("acetazolamide")
    assert data.get("_from_download") is True
    assert len(data.get("activities") or []) >= 2


def test_recommendation_is_derived_not_hardcoded():
    """Для вымышленной молекулы рекомендация берётся из её же строк, без подстановки."""
    record = {
        "ligand": "test_compound",
        "activities": [
            {"target_name": "Kinase X", "organism": "Homo sapiens", "affinity_type": "Ki",
             "value": 50.0, "relation": "=", "units": "nM"},
            {"target_name": "Kinase Y", "organism": "Rattus norvegicus", "affinity_type": "Ki",
             "value": 5.0, "relation": "=", "units": "nM"},
        ],
    }
    rec = recommend_target(record)
    # Человеческий белок приоритетнее более сильного измерения у крысы.
    assert rec["name"] == "Kinase X"
    assert rec["pdb_id"] is None
    assert rec["hypothesis"] is False
