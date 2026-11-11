"""Скачивание узких научных наборов: RCSB, BindingDB, PubChem."""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "data" / "external"
ACETAZOLAMIDE_SMILES = "CC(=O)Nc1nnc(S(N)(=O)=O)s1"
PDB_IDS = ("3HS4", "3PTB")


def _get(url: str, timeout: float = 60.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "DEARMO/1.0 (research)"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def fetch_pdb(pdb_id: str) -> Path:
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    dest = ROOT / "rcsb" / f"{pdb_id}.pdb"
    _write(dest, _get(url))
    meta_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    meta_dest = ROOT / "rcsb" / f"{pdb_id}.entry.json"
    _write(meta_dest, _get(meta_url))
    return dest


def fetch_pubchem_acetazolamide() -> Path:
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetazolamide/"
        "property/Title,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,IUPACName/JSON"
    )
    dest = ROOT / "pubchem" / "acetazolamide.json"
    _write(dest, _get(url))
    return dest


def _bindingdb_rows(payload: dict | list) -> list:
    if isinstance(payload, list):
        return payload
    for key in ("bdb.affinities", "affinities"):
        if isinstance(payload.get(key), list):
            return payload[key]
    for val in payload.values():
        if isinstance(val, dict) and isinstance(val.get("bdb.affinities"), list):
            return val["bdb.affinities"]
        if isinstance(val, list):
            return val
    return []


def _parse_affinity(raw) -> tuple[str, float | None]:
    text = str(raw).strip().replace(",", "")
    relation = "="
    if text[:1] in (">", "<"):
        relation = text[0]
        text = text[1:].strip()
    try:
        return relation, float(text)
    except ValueError:
        return relation, None


def _normalize_bindingdb(payload: dict | list) -> dict:
    activities = []
    for row in _bindingdb_rows(payload):
        if not isinstance(row, dict):
            continue
        relation, value = _parse_affinity(row.get("bdb.affinity") or row.get("affinity"))
        if value is None:
            continue
        activities.append(
            {
                "target_name": row.get("bdb.target") or row.get("target"),
                "organism": row.get("bdb.species") or row.get("species"),
                "affinity_type": row.get("bdb.affinity_type") or row.get("affinity_type") or "Ki",
                "value": value,
                "relation": relation,
                "units": "nM",
                "monomer_id": row.get("bdb.monomerid"),
            }
        )
    hit = None
    if isinstance(payload, dict):
        inner = next(iter(payload.values()), {})
        if isinstance(inner, dict):
            hit = inner.get("bdb.hit")
    return {
        "ligand": "acetazolamide",
        "source": "bindingdb",
        "hit_count": hit,
        "activities": activities,
    }


def fetch_bindingdb_acetazolamide() -> Path:
    smiles = urllib.parse.quote(ACETAZOLAMIDE_SMILES, safe="")
    url = (
        "https://bindingdb.org/rest/getTargetByCompound"
        f"?smiles={smiles}&cutoff=100000&response=application/json"
    )
    raw_path = ROOT / "bindingdb" / "acetazolamide_raw.json"
    body = _get(url, timeout=120)
    _write(raw_path, body)
    payload = json.loads(body.decode("utf-8", errors="replace"))
    normalized = _normalize_bindingdb(payload)
    normalized["downloaded_at"] = datetime.now(timezone.utc).isoformat()
    normalized["url"] = url
    dest = ROOT / "bindingdb" / "acetazolamide.json"
    dest.write_text(json.dumps(normalized, ensure_ascii=False, indent=2), encoding="utf-8")
    return dest


def fetch_all() -> dict[str, str]:
    written = {}
    for pdb_id in PDB_IDS:
        written[pdb_id] = str(fetch_pdb(pdb_id))
    written["pubchem"] = str(fetch_pubchem_acetazolamide())
    written["bindingdb"] = str(fetch_bindingdb_acetazolamide())
    manifest = ROOT / "MANIFEST.json"
    manifest.write_text(
        json.dumps({"files": written, "fetched_at": datetime.now(timezone.utc).isoformat()}, indent=2),
        encoding="utf-8",
    )
    return written


if __name__ == "__main__":
    paths = fetch_all()
    for key, path in paths.items():
        print(f"{key}: {path}")
