"""Проверка расчёта дескрипторов RDKit, QED, очистки солей, фильтрации и сходства Танимото."""
import pytest

from tools.rdkit_tools import (
    calculate_descriptors,
    evaluate_filters,
    parse_smiles,
    standardize_smiles,
    tanimoto_similarity,
)


def test_ethanol_descriptors():
    result = calculate_descriptors("CCO")
    assert result.status == "success"
    values = {item.method: item.value for item in result.evidence}
    assert values["molecular_weight"] == pytest.approx(46.069, abs=0.01)
    assert values["hbd"] == 1
    assert values["hba"] == 1
    assert "qed" in values
    assert 0.0 < values["qed"] <= 1.0
    assert all(item.context["rdkit_version"] for item in result.evidence)
    assert all(item.kind == "calculation" for item in result.evidence)


@pytest.mark.parametrize("smiles", ["", "   ", "not-a-smiles", "CCO ethanol", "C(C)(C"])
def test_invalid_input_has_no_results(smiles):
    result = calculate_descriptors(smiles)
    assert result.status == "error"
    assert result.errors
    assert result.evidence == []


def test_repeated_calculation_is_stable():
    first = calculate_descriptors("CCO")
    second = calculate_descriptors("CCO")
    assert [(e.method, e.value) for e in first.evidence] == [(e.method, e.value) for e in second.evidence]


def test_parse_smiles_rejects_empty_structure():
    with pytest.raises(ValueError):
        parse_smiles("")


def test_salt_stripping_and_limitation_note():
    """Соль ацетазоламида с [Na+] или Cl должна очищаться до родительской молекулы с записью в limitations."""
    salted = "CC(=O)Nc1nnc(S(N)(=O)=O)s1.[Na+]"
    parent = "CC(=O)Nc1nnc(S(N)(=O)=O)s1"

    std = standardize_smiles(salted, desalt=True)
    assert std["had_salt"] is True
    assert std["removed_fragments"]

    res_salted = calculate_descriptors(salted, desalt=True)
    res_parent = calculate_descriptors(parent, desalt=True)
    assert res_salted.status == "success"
    vals_salted = {e.method: e.value for e in res_salted.evidence}
    vals_parent = {e.method: e.value for e in res_parent.evidence}
    assert vals_salted["molecular_weight"] == pytest.approx(vals_parent["molecular_weight"], abs=0.05)
    assert vals_salted["qed"] == pytest.approx(vals_parent["qed"], abs=0.01)
    assert any("соль" in lim.lower() or "противоион" in lim.lower() for lim in res_salted.limitations)


def test_evaluate_filters_lipinski_veber_and_constraints():
    """Проверка прохождения фильтров для ацетазоламида и отсева при жёстких ограничениях."""
    res_ok = calculate_descriptors("CC(=O)Nc1nnc(S(N)(=O)=O)s1")
    filt_ok = res_ok.evidence[0].context["filter_evaluation"]
    assert filt_ok["passed"] is True
    assert filt_ok["lipinski_passed"] is True
    assert filt_ok["veber_passed"] is True
    assert filt_ok["failed_reasons"] == []

    # Задаём заведомо строгие ограничения, которые ацетазоламид (MW ~ 222.25, TPSA ~ 115) не проходит
    res_fail = calculate_descriptors(
        "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
        constraints={"mw_max": 150.0, "tpsa_max": 90.0, "qed_min": 0.9},
    )
    filt_fail = res_fail.evidence[0].context["filter_evaluation"]
    assert filt_fail["passed"] is False
    assert len(filt_fail["failed_reasons"]) >= 3
    assert any("mw_max" in r for r in filt_fail["failed_reasons"])
    assert any("tpsa_max" in r for r in filt_fail["failed_reasons"])
    assert any("qed_min" in r for r in filt_fail["failed_reasons"])
    assert any("не прошла фильтры" in lim.lower() for lim in res_fail.limitations)


def test_evaluate_filters_direct_lipinski_violation():
    heavy_desc = {
        "molecular_weight": 650.0,
        "logp": 6.2,
        "tpsa": 165.0,
        "hbd": 7,
        "hba": 12,
        "rotatable_bonds": 14,
        "qed": 0.15,
    }
    filt = evaluate_filters(heavy_desc)
    assert filt["passed"] is False
    assert filt["lipinski_passed"] is False
    assert filt["veber_passed"] is False
    assert len(filt["failed_reasons"]) >= 4


def test_tanimoto_similarity_identical_and_analogs():
    acetazolamide = "CC(=O)Nc1nnc(S(N)(=O)=O)s1"
    methazolamide = "CC(=O)N=c1sc(S(N)(=O)=O)nn1C"
    ethanol = "CCO"

    assert tanimoto_similarity(acetazolamide, acetazolamide) == pytest.approx(1.0)
    sim_analog = tanimoto_similarity(acetazolamide, methazolamide)
    sim_distant = tanimoto_similarity(acetazolamide, ethanol)
    assert 0.2 < sim_analog < 1.0
    assert sim_analog > sim_distant
