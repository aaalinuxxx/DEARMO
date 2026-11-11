from tools.bindingdb_client import fetch_activities


def test_bindingdb_cache_acetazolamide():
    data = fetch_activities("ацетазоламид")
    assert data.get("_from_download") is True
    assert len(data.get("activities") or []) >= 2
    assert data["recommended_target"]["pdb_id"] == "3HS4"
