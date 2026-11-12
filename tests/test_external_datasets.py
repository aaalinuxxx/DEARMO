"""Проверки на скачанных RCSB, BindingDB и PubChem."""
import json
from pathlib import Path

from tools.bindingdb_client import fetch_activities
from tools.pdb_client import get_structure_metadata
from tools.pubchem_client import resolve_compound

ROOT = Path(__file__).resolve().parent.parent / "data" / "external"


def test_downloaded_pdb_is_ca2_acetazolamide():
    pdb = (ROOT / "rcsb" / "3HS4.pdb").read_text(encoding="utf-8", errors="replace")
    assert "ACETAZOLAMIDE" in pdb.upper()
    assert "ATOM" in pdb
    meta = get_structure_metadata("3HS4")
    assert meta["resolution_angstrom"] == 1.1
    assert meta["pubmed_id"] == 19851004
    assert meta["_from_download"] is True


def test_bindingdb_recommendation_comes_from_measurements():
    """Мишень выводится из строк BindingDB, а не из зашитого соответствия молекуле."""
    record = fetch_activities("acetazolamide")
    assert record["_from_download"] is True
    assert len(record["activities"]) >= 20

    rec = record["recommended_target"]
    measured_names = {a.get("target_name") for a in record["activities"]}
    # Имя рекомендованной мишени обязано присутствовать среди измеренных строк.
    assert rec["name"] in measured_names
    # PDB не выдумывается: он либо пришёл из той же строки, либо отсутствует.
    assert rec.get("pdb_id") in {None, *(a.get("pdb_id") for a in record["activities"])}


def test_pubchem_download_cid():
    payload = json.loads((ROOT / "pubchem" / "acetazolamide.json").read_text(encoding="utf-8"))
    cid = payload["PropertyTable"]["Properties"][0]["CID"]
    assert cid == 1986
    resolved = resolve_compound("acetazolamide")
    assert resolved["cid"] == 1986
    assert resolved["source"] == "pubchem_download"
