"""Состав выбранной структуры и устойчивость обоснования к выдумкам модели."""
from __future__ import annotations

from agents import review_agent, structure_agent
from models.request import ResearchRequest
from models.result import Evidence
from orchestration.state import WorkflowState
from tools import pdb_client


def _state() -> WorkflowState:
    return WorkflowState(request=ResearchRequest(text="Подбери структуру"))


def test_metal_ion_is_separated_from_crystallization_artifacts():
    """Zn — кофактор, глицерин и вода — артефакты кристаллизации."""
    details = {
        "heteroatoms": [
            {"comp_id": "ZN", "role": "metal_ion", "name": "ZINC ION"},
            {"comp_id": "GOL", "role": "crystallization_artifact", "name": "GLYCEROL"},
            {"comp_id": "AZM", "role": "bound_ligand", "name": "acetazolamide"},
        ]
    }
    kept = {het["comp_id"] for het in pdb_client.significant_cofactors(details)}
    assert kept == {"ZN", "AZM"}


def test_chain_and_cofactor_become_separate_facts():
    state = _state()
    state.observations["uniprot_accession"] = "P00918"
    details = {
        "available": True,
        "chains": [
            {"auth_chains": ["A"], "name": "Carbonic anhydrase 2", "uniprot_accessions": ["P00918"]},
            {"auth_chains": ["B"], "name": "Другой белок", "uniprot_accessions": ["P99999"]},
        ],
        "heteroatoms": [{"comp_id": "ZN", "role": "metal_ion", "name": "ZINC ION"}],
    }

    structure_agent._report_composition(state, "3HS4", details)

    chains = [e.value for e in state.evidence if e.method == "structure_chain"]
    ions = [e.value for e in state.evidence if e.method == "metal_ion"]
    assert chains == ["A"], "цепь постороннего белка в ответ не идёт"
    assert ions == ["ZN"]


def test_rationale_claiming_an_absent_ligand_is_replaced():
    """Qwen писала «структура содержит запрошенный лиганд» там, где его нет."""
    state = _state()
    chosen = {"pdb_id": "3UYQ", "resolution_A": 1.7, "содержит_запрошенный_лиганд": "нет"}

    result = structure_agent._checked_rationale(
        state,
        "Выбрана структура, так как она содержит запрошенный лиганд и имеет лучшее разрешение.",
        chosen,
    )

    assert "содержит запрошенный лиганд" not in result.lower()
    assert "3UYQ" in result
    assert any("отклонено" in lim for lim in state.limitations)


def test_truthful_rationale_survives_untouched():
    state = _state()
    chosen = {"pdb_id": "3HS4", "resolution_A": 1.1, "содержит_запрошенный_лиганд": "да"}
    text = "Структура содержит запрошенный лиганд и имеет лучшее разрешение."

    assert structure_agent._checked_rationale(state, text, chosen) == text
    assert state.limitations == []


def test_review_rejects_an_invented_pdb_id():
    """Код структуры сочинить так же легко, как число, — значит, его надо проверять."""
    evidence = [Evidence(kind="database", value="3HS4", method="pdb_structure", source="RCSB PDB")]
    issues = review_agent._check_identifiers("Рекомендуется структура 1ABC для докинга.", evidence)
    assert any("1ABC" in issue for issue in issues)


def test_review_accepts_identifiers_that_are_backed_by_evidence():
    evidence = [
        Evidence(kind="database", value="3HS4", method="pdb_structure", source="RCSB PDB"),
        Evidence(kind="database", value="P00918", method="uniprot_accession", source="UniProt"),
    ]
    text = "Структура 3HS4 относится к белку P00918."
    assert review_agent._check_identifiers(text, evidence) == []


def test_digits_inside_an_accession_are_not_read_as_a_number():
    """P07451 — идентификатор, а не число 07451 без источника."""
    found = [m.group(0) for m in review_agent._NUMBER_IN_TEXT.finditer("Белок P07451 и структура 1M17.")]
    assert found == []
    found = [m.group(0) for m in review_agent._NUMBER_IN_TEXT.finditer("Разрешение 1.7 Å, Ki 3.1 нМ.")]
    assert found == ["1.7", "3.1"]
