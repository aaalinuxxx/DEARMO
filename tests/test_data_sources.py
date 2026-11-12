"""Тесты внешних источников данных: PubChem и ChEMBL (нормализация единиц в нМ, дедупликация, кэш, моки API)."""
from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from tools import chembl_client, pubchem_client


# ---------------------------------------------------------------------------
# Блок 1: ChEMBL — приведение единиц к нМ и обработка пропусков/несопоставимых
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "value,unit,expected_nm,expected_comp",
    [
        (12.5, "nM", 12.5, True),
        ("0.018", "uM", 18.0, True),
        (1.5, "µM", 1500.0, True),
        (0.002, "mM", 2000.0, True),
        (500, "pM", 0.5, True),
        (1e-8, "M", 10.0, True),
    ],
)
def test_chembl_normalize_to_nm_valid_units(value, unit, expected_nm, expected_comp):
    val_nm, norm_unit, comparable = chembl_client.normalize_to_nm(value, unit)
    assert comparable is expected_comp
    assert norm_unit == "nM"
    assert val_nm == pytest.approx(expected_nm, rel=1e-4)


@pytest.mark.parametrize("missing_val", [None, "", "   ", "N/A"])
def test_chembl_normalize_preserves_none_never_zero(missing_val):
    """Пропуски никогда не должны превращаться в 0.0 nM."""
    val_nm, _, comparable = chembl_client.normalize_to_nm(missing_val, "nM")
    assert val_nm is None
    assert comparable is False


@pytest.mark.parametrize(
    "value,incomparable_unit",
    [
        (85.0, "%"),
        (12.0, "ug/mL"),
        (5.0, "mg/kg"),
        (10.0, None),
    ],
)
def test_chembl_incomparable_units_flagged(value, incomparable_unit):
    """Несопоставимые единицы (%, ug/mL) не переводятся в нМ и помечаются comparable=False."""
    val_nm, norm_unit, comparable = chembl_client.normalize_to_nm(value, incomparable_unit)
    assert comparable is False
    assert norm_unit != "nM"
    assert val_nm == pytest.approx(float(value))


# ---------------------------------------------------------------------------
# Блок 2: ChEMBL — дедупликация и обнаружение противоречивых записей
# ---------------------------------------------------------------------------

def test_chembl_deduplicate_replicates_uses_median():
    """Близкие повторные измерения в одном ассее агрегируются медианой с сохранением исходного списка."""
    records = [
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "standard_type": "Ki",
            "standard_relation": "=",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "value_nm": 10.0,
            "normalized_units": "nM",
            "comparable": True,
            "document_chembl_id": "DOC1",
        },
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "standard_type": "Ki",
            "standard_relation": "=",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "value_nm": 12.0,
            "normalized_units": "nM",
            "comparable": True,
            "document_chembl_id": "DOC2",
        },
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "standard_type": "Ki",
            "standard_relation": "=",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "value_nm": 14.0,
            "normalized_units": "nM",
            "comparable": True,
            "document_chembl_id": "DOC2",
        },
    ]
    deduped, lims = chembl_client.deduplicate_activities(records)
    assert len(deduped) == 1
    assert deduped[0]["value_nm"] == pytest.approx(12.0)
    assert deduped[0]["replicates_count"] == 3
    assert deduped[0]["raw_values_nm"] == [10.0, 12.0, 14.0]
    assert deduped[0]["conflicting_measurements"] is False
    assert lims == []


def test_chembl_deduplicate_preserves_conflicting_records_without_silent_averaging():
    """При сильном расхождении (>=10x) записи НЕ усредняются молча, а помечаются conflicting_measurements=True."""
    records = [
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "standard_type": "Ki",
            "standard_relation": "=",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "value_nm": 10.0,
            "normalized_units": "nM",
            "comparable": True,
        },
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_organism": "Homo sapiens",
            "standard_type": "Ki",
            "standard_relation": "=",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "value_nm": 2500.0,
            "normalized_units": "nM",
            "comparable": True,
        },
    ]
    deduped, lims = chembl_client.deduplicate_activities(records, divergence_fold=10.0)
    assert len(deduped) == 2
    assert all(r["conflicting_measurements"] is True for r in deduped)
    assert any("противоречащие измерения" in lim.lower() or "расхождение" in lim.lower() for lim in lims)


def test_chembl_fetch_activities_with_mock_api(tmp_path, monkeypatch):
    """Проверка fetch_activities ChEMBL с моком API, нормализацией uM -> nM и сохранением кэша с датой."""
    monkeypatch.setattr(chembl_client, "_CACHE_DIR", tmp_path / "chembl_cache")
    mock_payload = {
        "activities": [
            {
                "molecule_chembl_id": "CHEMBL20",
                "target_chembl_id": "CHEMBL205",
                "target_pref_name": "Carbonic anhydrase II",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "CHEMBL658249",
                "assay_type": "B",
                "standard_type": "Ki",
                "standard_relation": "=",
                "standard_value": "0.012",
                "standard_units": "uM",
                "pchembl_value": "7.92",
                "document_chembl_id": "CHEMBL1145821",
                "data_validity_comment": None,
            }
        ]
    }
    with patch.object(chembl_client, "_get_json", return_value=mock_payload):
        rows = chembl_client.fetch_activities("CHEMBL205", "acetazolamide", use_cache=True)

    assert len(rows) == 1
    row = rows[0]
    assert row["molecule_chembl_id"] == "CHEMBL20"
    assert row["target_chembl_id"] == "CHEMBL205"
    assert row["value_nm"] == pytest.approx(12.0)
    assert row["normalized_units"] == "nM"
    assert row["comparable"] is True
    assert row["is_cached"] is False
    assert row["fetched_at"]

    # Второй вызов должен читаться из локального кэша с is_cached=True и сохранённой датой fetched_at
    with patch.object(chembl_client, "_get_json", return_value=None):
        cached_rows = chembl_client.fetch_activities("CHEMBL205", "acetazolamide", use_cache=True)
    assert len(cached_rows) == 1
    assert cached_rows[0]["is_cached"] is True
    assert cached_rows[0]["fetched_at"] == row["fetched_at"]


# ---------------------------------------------------------------------------
# Блок 3: PubChem — поиск по имени, SMILES (POST), CID, аналоги и неоднозначность
# ---------------------------------------------------------------------------

def test_pubchem_resolve_by_smiles_uses_post(tmp_path, monkeypatch):
    """Запрос по SMILES должен отправляться POST-запросом (чтобы спецсимволы не ломали URL)."""
    monkeypatch.setattr(pubchem_client, "_CACHE_DIR", tmp_path / "pubchem_cache")
    called = {}

    def fake_request_json(url, *, method="GET", form_data=None, timeout=10.0, max_retries=3):
        called["url"] = url
        called["method"] = method
        called["form_data"] = form_data
        return {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 1986,
                        "Title": "Acetazolamide",
                        "ConnectivitySMILES": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
                        "CanonicalSMILES": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
                        "MolecularWeight": 222.25,
                    }
                ]
            }
        }

    with patch.object(pubchem_client, "_request_json", side_effect=fake_request_json):
        res = pubchem_client.resolve_compound("CC(=O)Nc1nnc(S(N)(=O)=O)s1", use_cache=False)

    assert called["method"] == "POST"
    assert called["form_data"] == {"smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1"}
    assert res["cid"] == 1986
    assert res["query_mode"] == "smiles"
    assert res["ambiguous"] is False


def test_pubchem_ambiguous_name_returns_candidates():
    """Неоднозначное название (например, 'сульфонамид') возвращает флаг ambiguous=True и список кандидатов."""
    res = pubchem_client.resolve_compound("сульфонамид")
    assert res["ambiguous"] is True
    assert res["status"] == "needs_clarification"
    assert len(res["candidates"]) >= 2


def test_pubchem_search_similar_returns_analogs_with_tanimoto(tmp_path, monkeypatch):
    """Поиск похожих соединений возвращает список аналогов с рассчитанным Tanimoto similarity."""
    monkeypatch.setattr(pubchem_client, "_CACHE_DIR", tmp_path / "pubchem_sim_cache")
    with patch.object(pubchem_client, "_request_json", return_value=None):
        sim = pubchem_client.search_similar("CC(=O)Nc1nnc(S(N)(=O)=O)s1", threshold=80, max_records=4)

    assert sim["status"] == "success"
    assert len(sim["analogs"]) >= 2
    assert all("tanimoto_similarity" in a for a in sim["analogs"])
    assert sim["analogs"][0]["tanimoto_similarity"] >= sim["analogs"][-1]["tanimoto_similarity"]
    assert sim["fetched_at"]
