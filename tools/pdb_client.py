"""RCSB PDB: кеш метаданных и поиск по ключевым словам."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "pdb"
_EXTERNAL = Path(__file__).resolve().parent.parent / "data" / "external" / "rcsb"


def load_structure_cache(pdb_id: str) -> dict[str, Any] | None:
    path = _CACHE_ROOT / f"{pdb_id.upper()}.json"
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def _from_rcsb_entry(pdb_id: str) -> dict[str, Any] | None:
    path = _EXTERNAL / f"{pdb_id}.entry.json"
    pdb_path = _EXTERNAL / f"{pdb_id}.pdb"
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as f:
        entry = json.load(f)
    info = entry.get("rcsb_entry_info") or {}
    citation = entry.get("rcsb_primary_citation") or {}
    resolution = info.get("resolution_combined") or [None]
    return {
        "pdb_id": pdb_id,
        "title": (entry.get("struct") or {}).get("title"),
        "method": (info.get("experimental_method") or "").replace("X-ray", "X-RAY DIFFRACTION"),
        "resolution_angstrom": resolution[0] if resolution else None,
        "pubmed_id": (entry.get("rcsb_entry_container_identifiers") or {}).get("pubmed_id"),
        "doi": citation.get("pdbx_database_id_DOI"),
        "local_pdb": str(pdb_path) if pdb_path.is_file() else None,
        "source": "rcsb_download",
        "_from_download": True,
    }


def get_structure_metadata(pdb_id: str) -> dict[str, Any]:
    pdb_id = pdb_id.upper()[:4]
    downloaded = _from_rcsb_entry(pdb_id)
    if downloaded:
        return downloaded
    cached = load_structure_cache(pdb_id)
    if cached:
        cached["_from_cache"] = True
        return cached
    return {
        "pdb_id": pdb_id,
        "title": None,
        "source": "rcsb",
        "note": f"Нет кеша для {pdb_id}",
        "_from_cache": False,
    }


def search_by_ligand_target(ligand_key: str, target_hint: str | None) -> list[dict[str, Any]]:
    index_path = _CACHE_ROOT / "recommendations.json"
    if not index_path.is_file():
        return []
    with index_path.open(encoding="utf-8") as f:
        data = json.load(f)
    entries = data.get(ligand_key.lower(), [])
    if target_hint:
        hint = target_hint.lower()
        entries = [e for e in entries if hint in (e.get("protein_name") or "").lower()]
    return entries
