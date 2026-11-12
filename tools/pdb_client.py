"""RCSB PDB: кэш метаданных с датой получения (fetched_at), устойчивый REST-запрос (timeout + retry) и поиск структур."""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "pdb"
_EXTERNAL = Path(__file__).resolve().parent.parent / "data" / "external" / "rcsb"
_RCSB_REST = "https://data.rcsb.org/rest/v1/core/entry"


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _file_mtime_iso(path: Path) -> str:
    try:
        return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except OSError:
        return _now_iso()


def _get_json_with_retry(url: str, *, timeout: float = 10.0, max_retries: int = 3) -> dict[str, Any] | None:
    headers = {"Accept": "application/json"}
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=headers, method="GET")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504) and attempt + 1 < max_retries:
                time.sleep(0.2 * (attempt + 1))
                continue
            return None
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError):
            if attempt + 1 < max_retries:
                time.sleep(0.1 * (attempt + 1))
                continue
            return None
    return None


def load_structure_cache(pdb_id: str) -> dict[str, Any] | None:
    path = _CACHE_ROOT / f"{pdb_id.upper()}.json"
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            data.setdefault("fetched_at", _file_mtime_iso(path))
            data["is_cached"] = True
            data["_from_cache"] = True
            return data
    except (OSError, json.JSONDecodeError):
        return None
    return None


def save_structure_cache(pdb_id: str, metadata: dict[str, Any]) -> None:
    try:
        _CACHE_ROOT.mkdir(parents=True, exist_ok=True)
        path = _CACHE_ROOT / f"{pdb_id.upper()}.json"
        payload = dict(metadata)
        payload.setdefault("fetched_at", _now_iso())
        with path.open("w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
    except OSError:
        pass


def _parse_rcsb_entry(pdb_id: str, entry: dict[str, Any], *, pdb_path: Path | None = None) -> dict[str, Any]:
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
        "local_pdb": str(pdb_path) if (pdb_path and pdb_path.is_file()) else None,
    }


def _from_rcsb_entry(pdb_id: str) -> dict[str, Any] | None:
    path = _EXTERNAL / f"{pdb_id}.entry.json"
    pdb_path = _EXTERNAL / f"{pdb_id}.pdb"
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            entry = json.load(f)
        parsed = _parse_rcsb_entry(pdb_id, entry, pdb_path=pdb_path)
        parsed.update(
            {
                "source": "rcsb_download",
                "_from_download": True,
                "_from_cache": False,
                "is_cached": True,
                "fetched_at": _file_mtime_iso(path),
            }
        )
        return parsed
    except (OSError, json.JSONDecodeError):
        return None


def get_structure_metadata(
    pdb_id: str,
    *,
    use_cache: bool = True,
    timeout: float = 10.0,
) -> dict[str, Any]:
    """Получить метаданные структуры PDB (из предзагрузки, кэша или RCSB Data API с таймаутом и ретраями)."""
    clean_id = pdb_id.strip().upper()[:4]

    if use_cache:
        downloaded = _from_rcsb_entry(clean_id)
        if downloaded:
            return downloaded
        cached = load_structure_cache(clean_id)
        if cached:
            return cached

    # Live-запрос к RCSB Data API с таймаутом 10с и 3 попытками
    url = f"{_RCSB_REST}/{clean_id}"
    live_entry = _get_json_with_retry(url, timeout=timeout, max_retries=3)
    if live_entry:
        fetched_at = _now_iso()
        parsed = _parse_rcsb_entry(clean_id, live_entry)
        parsed.update(
            {
                "source": "rcsb_live",
                "_from_download": False,
                "_from_cache": False,
                "is_cached": False,
                "fetched_at": fetched_at,
            }
        )
        save_structure_cache(clean_id, parsed)
        return parsed

    return {
        "pdb_id": clean_id,
        "title": None,
        "source": "rcsb",
        "note": f"Нет локального кэша и не удалось получить данные RCSB для {clean_id}",
        "_from_cache": False,
        "_from_download": False,
        "is_cached": False,
        "fetched_at": _now_iso(),
    }


def search_by_ligand_target(ligand_key: str, target_hint: str | None = None) -> list[dict[str, Any]]:
    index_path = _CACHE_ROOT / "recommendations.json"
    if not index_path.is_file():
        return []
    try:
        with index_path.open(encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return []
    entries = data.get(ligand_key.lower(), [])
    if target_hint:
        hint = target_hint.lower()
        entries = [e for e in entries if hint in (e.get("protein_name") or "").lower()]
    return entries
