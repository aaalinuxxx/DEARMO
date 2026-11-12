"""Защита от возврата хардкода.

Кейс требует «отсутствия вымышленных результатов». Самый незаметный способ их
получить — зашитый в код справочник структур, CID или активностей: он молча
подставляет правдоподобные числа там, где источник ничего не ответил.
Эти тесты падают, если такой справочник появится снова.
"""
from __future__ import annotations

import re
from pathlib import Path
from unittest.mock import patch

import pytest

from tools import chembl_client, pubchem_client

_ROOT = Path(__file__).resolve().parent.parent
_CODE_DIRS = ("agents", "orchestration", "tools", "models")
# Скрипт разовой загрузки опорных данных вправе называть конкретные молекулы.
_ALLOWED_FILES = {"dataset_fetch.py"}

# Длинная строка, похожая на SMILES, в коде = зашитая структура молекулы.
_SMILES_LITERAL = re.compile(r'"[A-Za-z0-9@+\-\[\]()=#\\/%]{12,}"')
_SMILES_MARKERS = ("=O", "c1", "C1", "N1", "S(")


def _python_files() -> list[Path]:
    files: list[Path] = []
    for folder in _CODE_DIRS:
        files.extend(
            path
            for path in (_ROOT / folder).rglob("*.py")
            if path.name not in _ALLOWED_FILES and "__pycache__" not in str(path)
        )
    return files


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: str(p.relative_to(_ROOT)))
def test_no_hardcoded_smiles_in_source(path: Path) -> None:
    """В коде нет зашитых структур молекул: их источник — PubChem или пользователь."""
    offenders = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("#") or stripped.startswith('"""'):
            continue
        for match in _SMILES_LITERAL.finditer(line):
            literal = match.group(0)
            if any(marker in literal for marker in _SMILES_MARKERS):
                offenders.append(stripped[:100])
    assert not offenders, f"Зашитый SMILES в {path.name}: {offenders}"


def test_unknown_molecule_does_not_get_a_substituted_cid(tmp_path, monkeypatch) -> None:
    """Если PubChem не ответил, CID остаётся пустым, а не берётся у другой молекулы."""
    monkeypatch.setattr(pubchem_client, "_CACHE_DIR", tmp_path / "cache")
    monkeypatch.setattr(pubchem_client, "_PUBCHEM_DIR", tmp_path / "external")
    with patch.object(pubchem_client, "_request_json", return_value=None):
        result = pubchem_client.resolve_compound("nonexistent-compound-xyz")
    assert result.get("cid") is None
    assert result.get("error") == "no_data"


def test_chembl_without_api_returns_nothing_not_reference_values(monkeypatch, tmp_path) -> None:
    """Недоступный ChEMBL даёт пустой результат и признак недоступности, а не «эталонные» числа."""
    monkeypatch.setattr(chembl_client, "_CACHE_DIR", tmp_path / "chembl")
    chembl_client.reset_transport_status()
    with patch.object(chembl_client, "_get_json", return_value=None):
        rows = chembl_client.fetch_activities("CHEMBL205", "CHEMBL20", use_cache=False)
    assert rows == []


def test_target_id_is_not_defaulted_to_any_protein() -> None:
    """Пустая мишень остаётся пустой: подстановки «мишени по умолчанию» нет."""
    assert chembl_client.resolve_target_id("") == ""
    assert chembl_client.resolve_target_id("   ") == ""
