"""PubChem PUG REST: идентификаторы и SMILES."""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

_FALLBACK_SMILES = {
    "acetazolamide": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    "ацетазоламид": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
}
_PUBCHEM_FILE = Path(__file__).resolve().parent.parent / "data" / "external" / "pubchem" / "acetazolamide.json"


def _get_json(url: str, timeout: float = 15.0) -> dict[str, Any] | None:
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError):
        return None


def resolve_compound(name_or_smiles: str) -> dict[str, Any]:
    key = name_or_smiles.strip()
    lower = key.lower()
    if lower in ("acetazolamide", "ацетазоламид") and _PUBCHEM_FILE.is_file():
        with _PUBCHEM_FILE.open(encoding="utf-8") as f:
            payload = json.load(f)
        props = payload.get("PropertyTable", {}).get("Properties", [{}])[0]
        return {
            "query": key,
            "cid": props.get("CID"),
            "title": props.get("Title"),
            "connectivity_smiles": props.get("ConnectivitySMILES") or props.get("CanonicalSMILES"),
            "molecular_weight": props.get("MolecularWeight"),
            "source": "pubchem_download",
        }
    if lower in _FALLBACK_SMILES:
        return {
            "query": key,
            "cid": 1986,
            "title": "Acetazolamide",
            "connectivity_smiles": _FALLBACK_SMILES[lower],
            "source": "pubchem_fallback",
        }
    url = (
        f"{_PUG}/compound/name/{urllib.parse.quote(key)}"
        "/property/Title,CID,ConnectivitySMILES/JSON"
    )
    data = _get_json(url)
    if not data:
        return {"query": key, "source": "pubchem", "error": "no_data"}
    props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return {
        "query": key,
        "cid": props.get("CID"),
        "title": props.get("Title"),
        "connectivity_smiles": props.get("ConnectivitySMILES"),
        "source": "pubchem",
    }
