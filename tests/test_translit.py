"""Русские названия молекул без словаря соединений."""
from __future__ import annotations

import pytest

from tools.translit import candidates, looks_russian, resolve_russian_name


def test_detects_cyrillic():
    assert looks_russian("эрлотиниба")
    assert not looks_russian("erlotinib")


@pytest.mark.parametrize(
    "word,expected",
    [
        ("эрлотиниба", "erlotinib"),
        ("иматиниб", "imatinib"),
        ("целекоксиб", "celecoxib"),
        ("аспирина", "aspirin"),
    ],
)
def test_candidate_list_contains_correct_spelling(word, expected):
    """Правильное написание обязано попасть в список вариантов, и не в самый хвост."""
    assert expected in candidates(word)


def test_resolution_requires_exact_pubchem_hit():
    """Вариант принимается только при точном совпадении: приблизительных нет."""
    seen: list[str] = []

    def resolver(name: str):
        seen.append(name)
        return {"cid": 5291, "title": "Imatinib"} if name == "imatinib" else {"cid": None}

    result = resolve_russian_name("иматиниба", resolver)
    assert result is not None
    assert result["cid"] == 5291
    assert result["translit_variant"] == "imatinib"


def test_unresolvable_word_is_not_guessed():
    """Если PubChem ничего не подтвердил, молекула не выдумывается."""
    assert resolve_russian_name("правдоподобность", lambda name: {"cid": None}) is None


def test_ambiguous_hit_is_rejected():
    """Неоднозначный ответ PubChem не принимается за найденную молекулу."""
    assert resolve_russian_name("вещество", lambda name: {"cid": 1, "ambiguous": True}) is None
