"""BindingDB: чтение кеша и опциональный live API."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "bindingdb"
_EXTERNAL = Path(__file__).resolve().parent.parent / "data" / "external" / "bindingdb" / "acetazolamide.json"


def cache_path(ligand_key: str) -> Path:
    return _CACHE_ROOT / f"{ligand_key.lower()}.json"


def load_cached_activities(ligand_key: str) -> dict[str, Any] | None:
    path = cache_path(ligand_key)
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as f:
        return json.load(f)


_ALIASES = {
    "ацетазоламид": "acetazolamide",
    "acetazolamid": "acetazolamide",
}


def recommend_target(record: dict[str, Any]) -> dict[str, Any] | None:
    """Мишень для докинга — CA II по кристаллу 3HS4, не минимум Ki по всем изоформам."""
    rows = record.get("activities") or []
    human_ca2 = [
        a for a in rows
        if "carbonic anhydrase 2" in (a.get("target_name") or "").lower()
        and "human" in (a.get("organism") or "").lower()
        and a.get("relation") == "="
    ]
    if not rows and record.get("recommended_target"):
        return record["recommended_target"]
    mutant = next((a for a in human_ca2 if "[" in (a.get("target_name") or "")), None)
    chosen = mutant or (human_ca2[0] if human_ca2 else None)
    hypothesis = chosen is None or "[" in (chosen.get("target_name") or "")
    return {
        "name": "Carbonic anhydrase 2",
        "organism": "Homo sapiens",
        "pdb_id": "3HS4",
        "hypothesis": hypothesis,
        "affinity_type": None if chosen is None else chosen.get("affinity_type"),
        "value": None if chosen is None else chosen.get("value"),
        "units": None if chosen is None else chosen.get("units"),
        "measured_target_name": None if chosen is None else chosen.get("target_name"),
        "rationale": (
            "Структура 3HS4 — комплекс человеческой CA II с ацетазоламидом (RCSB). "
            "В выгрузке BindingDB нет Ki дикого типа CA II; ближайшая строка — мутант, "
            "поэтому аффинность не переносится на WT как факт."
            if hypothesis
            else "Есть измерение по человеческой CA II без пометки мутации."
        ),
    }


def fetch_activities(ligand_name: str, *, use_cache: bool = True) -> dict[str, Any]:
    raw = ligand_name.lower().replace(" ", "_")
    key = _ALIASES.get(raw, raw)
    if use_cache and key == "acetazolamide" and _EXTERNAL.is_file():
        with _EXTERNAL.open(encoding="utf-8") as f:
            record = json.load(f)
        record["_from_cache"] = False
        record["_from_download"] = True
        rec = recommend_target(record)
        if rec:
            record["recommended_target"] = rec
        return record
    if use_cache:
        cached = load_cached_activities(key)
        if cached:
            cached["_from_cache"] = True
            return cached
    return {
        "ligand": ligand_name,
        "activities": [],
        "source": "bindingdb",
        "note": "Нет локального кеша; live API не подключён в MVP.",
        "_from_cache": False,
    }
