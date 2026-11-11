"""Проверка расчёта и ошибочных входов на первом этапе."""
import pytest
from tools.rdkit_tools import calculate_descriptors, parse_smiles


def test_ethanol_descriptors():
    result = calculate_descriptors("CCO")
    assert result.status == "success"
    values = {item.method: item.value for item in result.evidence}
    assert values["molecular_weight"] == pytest.approx(46.069, abs=0.01)
    assert values["hbd"] == 1
    assert values["hba"] == 1
    assert all(item.context["rdkit_version"] for item in result.evidence)
    assert all(item.kind == "calculation" for item in result.evidence)

@pytest.mark.parametrize("smiles", ["", "   ", "not-a-smiles", "CCO ethanol"])
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
