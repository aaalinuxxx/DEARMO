"""Мишень и её идентификатор: выбор по данным, accession только из UniProt."""
from __future__ import annotations

import pytest

from agents import target_agent
from models.request import ResearchRequest
from orchestration.state import WorkflowState
from tools import uniprot_client

_ROWS = [
    {"target_name": "Carbonic anhydrase", "organism": "Yeast", "affinity_type": "Ki",
     "value": 11.0, "relation": "=", "units": "nM"},
    {"target_name": "Carbonic anhydrase 2", "organism": "Human", "affinity_type": "Ki",
     "value": 25.0, "relation": "=", "units": "nM"},
    {"target_name": "Carbonic anhydrase 2", "organism": "Human", "affinity_type": "Ki",
     "value": 30.0, "relation": "=", "units": "nM"},
    {"target_name": "Carbonic anhydrase", "organism": "Vibrio cholerae", "affinity_type": "Ki",
     "value": 6.8, "relation": "=", "units": "nM"},
]


def _state(text: str = "Подбери белок-мишень") -> WorkflowState:
    return WorkflowState(request=ResearchRequest(text=text))


def test_human_target_is_not_pushed_out_by_a_stronger_ortholog():
    """Сильное связывание у бактериального ортолога не должно вытеснять человеческий белок."""
    options = target_agent.group_targets(_ROWS)
    names = [(o["target_name"], o["organism"]) for o in options]
    assert names[0] == ("Carbonic anhydrase 2", "Human")
    # Остальные организмы из списка не выбрасываются, только отодвигаются.
    assert ("Carbonic anhydrase", "Vibrio cholerae") in names


def test_grouping_keeps_the_best_measurement_per_target():
    options = target_agent.group_targets(_ROWS)
    human = next(o for o in options if o["organism"] == "Human")
    assert human["best_affinity_nM"] == 25.0
    assert human["measurements"] == 2


def test_llm_target_hint_without_measurements_is_rejected():
    """Мишень, которую выдумал планировщик, не принимается молча.

    Наблюдали на практике: по запросу про ацетазоламид Qwen извлекла мишень
    ACO1_HUMAN, которой в запросе не было. Такая подсказка обязана отпасть.
    """
    # Пользователь назвал мишень сам, но по ней нет ни одного измерения этого соединения.
    state = _state("Подбери лиганды для ACO1_HUMAN")
    state.observations["bindingdb"] = {"activities": _ROWS}
    state.observations["target_hint"] = "ACO1_HUMAN"
    state.observations["target_hint_source"] = "ollama:qwen2.5:7b-instruct"

    chosen = target_agent._honour_named_target(state, _ROWS)

    assert chosen is None, "подсказка без измерений не может быть принята"
    assert state.observations["rejected_target_hint"] == "ACO1_HUMAN"
    assert any("отброшена" in lim for lim in state.limitations)


def test_llm_target_hint_backed_by_measurements_is_kept():
    # Мишень названа самим пользователем, модель лишь её извлекла.
    state = _state("Подбери лиганды для Carbonic anhydrase 2")
    state.observations["target_hint"] = "Carbonic anhydrase 2"
    state.observations["target_hint_source"] = "ollama:qwen2.5:7b-instruct"

    chosen = target_agent._honour_named_target(state, _ROWS)

    assert chosen["target_name"] == "Carbonic anhydrase 2"
    assert chosen["organism"] == "Human"
    assert state.observations["target_choice_source"] == "planner_llm_confirmed"


def test_user_named_target_is_honoured_even_without_data():
    """Пользователь спросил про конкретный белок — отвечаем про него, а не про другой."""
    state = WorkflowState(request=ResearchRequest(text="активность по тромбину", target_id="Thrombin"))

    chosen = target_agent._honour_named_target(state, _ROWS)

    assert chosen["target_name"] == "Thrombin"
    assert state.observations["target_choice_source"] == "user"


def test_accession_is_never_guessed_when_uniprot_says_nothing(monkeypatch):
    """Нет подтверждения из UniProt — нет accession. Память модели не используется."""
    monkeypatch.setattr(target_agent, "resolve_target_parts", lambda *a, **k: [])
    state = _state()

    record = target_agent.resolve_identity(state, {"target_name": "Неизвестный белок", "organism": None})

    assert record is None
    assert "uniprot_accession" not in state.observations
    assert any("не сопоставлена записи UniProt" in lim for lim in state.limitations)


def test_uniprot_rejects_an_entry_whose_name_does_not_match(monkeypatch):
    """Выдача поиска, не совпавшая по имени, отбрасывается, а не принимается первой строкой."""
    monkeypatch.setattr(
        uniprot_client,
        "request_json",
        lambda *a, **k: {
            "results": [
                {
                    "primaryAccession": "P99999",
                    "uniProtkbId": "OTHER_HUMAN",
                    "proteinDescription": {
                        "recommendedName": {"fullName": {"value": "Completely different enzyme"}}
                    },
                    "organism": {"scientificName": "Homo sapiens"},
                }
            ]
        },
    )
    assert uniprot_client.resolve_target("Carbonic anhydrase 2", "Homo sapiens") is None


@pytest.mark.parametrize(
    "query,entry_name,expected",
    [
        ("Carbonic anhydrase 2", "Carbonic anhydrase 2", "exact"),
        # «domain-containing protein» — шум номенклатуры, он не должен мешать совпадению.
        ("Alpha-carbonic anhydrase domain-containing protein", "Alpha-carbonic anhydrase", "exact"),
        ("Carbonic anhydrase", "Carbonic anhydrase 2", "partial"),
        ("Thrombin", "Carbonic anhydrase 2", None),
    ],
)
def test_name_match_quality(query, entry_name, expected):
    entry = {
        "proteinDescription": {"recommendedName": {"fullName": {"value": entry_name}}},
    }
    assert uniprot_client._match_quality(query, entry) == expected


def test_fused_target_name_yields_both_halves(monkeypatch):
    """BCR-ABL: у слияния нет своей записи, а у половин есть. Выбор здесь не делается."""
    seen: list[str] = []

    def fake_resolve_one(name, organism, *, timeout):
        seen.append(name)
        if name == "Breakpoint cluster region protein":
            return {"accession": "P11274", "gene": "BCR", "name_match": "exact"}
        if name == "Tyrosine-protein kinase ABL1":
            return {"accession": "P00519", "gene": "ABL1", "name_match": "exact"}
        return None

    monkeypatch.setattr(uniprot_client, "_resolve_one", fake_resolve_one)
    found = uniprot_client.resolve_target_parts(
        "Breakpoint cluster region protein/Tyrosine-protein kinase ABL1", "Human"
    )

    assert [r["accession"] for r in found] == ["P11274", "P00519"]
    assert all(r["full_target_name"].startswith("Breakpoint") for r in found)


# ---- структурное свидетельство и мутанты ----

def test_mutation_tail_is_split_from_the_protein_name():
    """«[A65S,N67Q]» — мутант того же белка, а не другой белок: имя нужно для UniProt."""
    assert target_agent.split_variant("Carbonic anhydrase 2 [A65S,N67Q,E69T]") == (
        "Carbonic anhydrase 2",
        "A65S,N67Q,E69T",
    )
    assert target_agent.split_variant("Carbonic anhydrase 15 [19-324]") == ("Carbonic anhydrase 15", "19-324")
    assert target_agent.split_variant("Thrombin") == ("Thrombin", None)


def test_variant_rows_merge_into_the_wild_type_candidate():
    rows = [
        {"target_name": "Carbonic anhydrase 2 [A65S,N67Q]", "organism": "Human", "affinity_type": "Ki",
         "value": 25.0, "relation": "=", "units": "nM"},
        {"target_name": "Carbonic anhydrase 2", "organism": "Human", "affinity_type": "Ki",
         "value": 12.0, "relation": "=", "units": "nM"},
    ]
    options = target_agent.group_targets(rows)
    assert len(options) == 1
    assert options[0]["measurements"] == 2
    assert options[0]["variants"] == ["A65S,N67Q"]


def test_crystal_complexes_outrank_a_stronger_binder_without_structures(monkeypatch):
    """Для ацетазоламида в BindingDB сильнее CA3, но кристаллы у CA2: решают независимые данные."""
    options = [
        {"target_name": "Carbonic anhydrase 3", "organism": "Human", "best_affinity_nM": 3.1,
         "affinity_type": "Ki", "measurements": 2, "variants": [], "link": None},
        {"target_name": "Carbonic anhydrase 2", "organism": "Human", "best_affinity_nM": 25.0,
         "affinity_type": "Ki", "measurements": 1, "variants": ["A65S"], "link": None},
    ]
    accession = {"Carbonic anhydrase 3": "P07451", "Carbonic anhydrase 2": "P00918"}
    complexes = {"P07451": 0, "P00918": 14}
    monkeypatch.setattr(
        target_agent, "resolve_target_parts",
        lambda name, organism=None, **k: [{"accession": accession[name], "gene": "x"}],
    )
    monkeypatch.setattr(target_agent, "count_complexes", lambda acc, comps, **k: complexes[acc])
    state = _state()
    state.observations["ligand_comp_ids"] = ["AZM"]

    target_agent.enrich_with_structures(state, options)

    assert options[0]["target_name"] == "Carbonic anhydrase 2"
    assert options[0]["комплексов_в_PDB_с_этим_лигандом"] == 14


def test_unknown_complex_count_does_not_reorder_candidates(monkeypatch):
    """Нет связи с RCSB — значит неизвестно, а не ноль: прежний порядок сохраняется."""
    options = [
        {"target_name": "A", "organism": "Human", "best_affinity_nM": 1.0, "affinity_type": "Ki",
         "measurements": 1, "variants": [], "link": None},
        {"target_name": "B", "organism": "Human", "best_affinity_nM": 9.0, "affinity_type": "Ki",
         "measurements": 1, "variants": [], "link": None},
    ]
    monkeypatch.setattr(target_agent, "resolve_target_parts", lambda *a, **k: [])
    monkeypatch.setattr(target_agent, "count_complexes", lambda *a, **k: None)
    state = _state()
    state.observations["ligand_comp_ids"] = ["AZM"]
    target_agent.enrich_with_structures(state, options)
    assert [o["target_name"] for o in options] == ["A", "B"]
    assert all("комплексов_в_PDB_с_этим_лигандом" not in o for o in options)


def test_claim_of_strongest_binding_is_checked_against_the_numbers():
    """Qwen писала «наименьшая константа (25 нМ)», когда у другого кандидата было 3.1 нМ."""
    state = _state()
    chosen = {"target_name": "CA2", "organism": "Human", "best_affinity_nM": 25.0,
              "комплексов_в_PDB_с_этим_лигандом": 14}
    options = [chosen, {"target_name": "CA3", "organism": "Human", "best_affinity_nM": 3.1}]

    text = target_agent._checked_rationale(
        state, "Выбран из-за наименьшей константы связывания (25.0 нМ).", chosen, options
    )

    assert "наименьшей константы" not in text
    assert "3.1" in text and "14" in text
    assert any("отклонено" in lim for lim in state.limitations)


def test_true_claim_of_strongest_binding_is_left_alone():
    state = _state()
    chosen = {"target_name": "CA3", "organism": "Human", "best_affinity_nM": 3.1}
    options = [chosen, {"target_name": "CA2", "organism": "Human", "best_affinity_nM": 25.0}]
    text = "Выбран из-за наименьшей константы связывания."
    assert target_agent._checked_rationale(state, text, chosen, options) == text
    assert state.limitations == []


def test_mutant_caveat_travels_with_the_answer_text():
    """Оговорка про мутант стоит рядом с числом, а не только в списке ограничений."""
    from orchestration.workflow import _caveat_line, _sources_line
    from models.result import Evidence

    state = _state()
    state.limitations.append("Измерения BindingDB по «X» относятся к мутантным вариантам или фрагментам (A65S).")
    state.add_evidence(Evidence(kind="database", value=25.0, unit="nM", method="Ki", source="BindingDB",
                                context={"link": "https://www.bindingdb.org/x?monomerid=1"}))
    assert _caveat_line(state).startswith("Важно:")
    assert "bindingdb.org" in _sources_line(state)


def test_urls_added_by_code_do_not_trip_the_number_check():
    """monomerid=10880 в адресе — не числовое утверждение модели."""
    from agents import review_agent
    from models.result import Evidence

    state = _state()
    state.scenario_id = "target_for_ligand"
    state.add_evidence(Evidence(kind="database", value=25.0, unit="nM", method="Ki", source="BindingDB"))
    state.explanation = "Ki составляет 25.0 нМ.\n\nЗаписи: https://www.bindingdb.org/x?monomerid=10880"
    review_agent.review(state)
    assert state.review_passed, state.review_errors


def test_model_invented_hint_does_not_replace_data_driven_choice():
    """Регрессия: Qwen вписывала «Carbonic anhydrase» в запрос про ацетазоламид.

    Под это имя есть измерения у бактериальных ортологов, подсказка «подтверждалась»
    и подменяла выбор по данным: ответ приходил про Vibrio cholerae.
    """
    state = _state("Подбери белок-мишень для ацетазоламида и обоснуй выбор экспериментальными данными.")
    state.observations["target_hint"] = "Carbonic anhydrase"
    state.observations["target_hint_source"] = "ollama:qwen2.5:7b-instruct"
    rows = _ROWS + [
        {"target_name": "Carbonic anhydrase", "organism": "Vibrio cholerae", "affinity_type": "Ki",
         "value": 6.8, "relation": "=", "units": "nM"},
    ]

    chosen = target_agent._honour_named_target(state, rows)

    assert chosen is None, "выбор по данным должен остаться за target_agent"
    assert state.observations["ignored_target_hint"] == "Carbonic anhydrase"
    assert "target_hint" not in state.observations
    assert any("не принята" in lim for lim in state.limitations)


def test_acetazolamide_style_invented_hint_is_ignored_not_rejected():
    """Подсказка, которой нет в запросе, отбрасывается как «не принята», а не как «нет измерений»."""
    state = _state("Подбери белок-мишень для ацетазоламида")
    state.observations["target_hint"] = "ACO1_HUMAN"
    state.observations["target_hint_source"] = "ollama:qwen2.5:7b-instruct"
    assert target_agent._honour_named_target(state, _ROWS) is None
    assert state.observations["ignored_target_hint"] == "ACO1_HUMAN"
    assert "rejected_target_hint" not in state.observations
