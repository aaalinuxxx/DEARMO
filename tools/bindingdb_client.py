"""BindingDB: чтение кэша с метаданными даты (fetched_at), устойчивые HTTP-запросы (timeout + retry) и выбор мишени."""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "bindingdb"
_EXTERNAL = Path(__file__).resolve().parent.parent / "data" / "external" / "bindingdb" / "acetazolamide.json"

_ALIASES = {
    "ацетазоламид": "acetazolamide",
    "acetazolamid": "acetazolamide",
    "cc(=o)nc1nnc(s(n)(=o)=o)s1": "acetazolamide",
    "метазоламид": "methazolamide",
    "этоксзоламид": "ethoxzolamide",
    "бензамидин": "benzamidine",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _file_mtime_iso(path: Path) -> str:
    try:
        return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except OSError:
        return _now_iso()


def cache_path(ligand_key: str) -> Path:
    return _CACHE_ROOT / f"{ligand_key.lower()}.json"


def load_cached_activities(ligand_key: str) -> dict[str, Any] | None:
    path = cache_path(ligand_key)
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


def save_cached_activities(ligand_key: str, record: dict[str, Any]) -> None:
    try:
        _CACHE_ROOT.mkdir(parents=True, exist_ok=True)
        path = cache_path(ligand_key)
        payload = dict(record)
        payload.setdefault("fetched_at", _now_iso())
        with path.open("w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
    except OSError:
        pass


def _get_json_with_retry(url: str, *, timeout: float = 10.0, max_retries: int = 3) -> Any | None:
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
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError, ValueError):
            if attempt + 1 < max_retries:
                time.sleep(0.1 * (attempt + 1))
                continue
            return None
    return None


def recommend_target(record: dict[str, Any]) -> dict[str, Any] | None:
    """Мишень для докинга — CA II по кристаллу 3HS4, не минимум Ki по всем изоформам."""
    rows = record.get("activities") or []
    human_ca2 = [
        a
        for a in rows
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


def fetch_activities(
    ligand_name: str,
    *,
    use_cache: bool = True,
    timeout: float = 10.0,
) -> dict[str, Any]:
    """Получить записи активности из BindingDB (локальная предзагрузка, файловый кэш или REST API с повторами)."""
    raw = ligand_name.strip().lower().replace(" ", "_")
    key = _ALIASES.get(raw, raw)

    # 1. Предзагруженный эталонный файл репозитория (сохраняем обратную совместимость _from_download=True)
    if use_cache and key == "acetazolamide" and _EXTERNAL.is_file():
        with _EXTERNAL.open(encoding="utf-8") as f:
            record = json.load(f)
        fetched_at = record.get("fetched_at") or _file_mtime_iso(_EXTERNAL)
        record["_from_cache"] = False
        record["_from_download"] = True
        record["is_cached"] = True
        record["fetched_at"] = fetched_at
        rec = recommend_target(record)
        if rec:
            record["recommended_target"] = rec
        return record

    # 2. Локальный кэш в data/cache/bindingdb/
    if use_cache:
        cached = load_cached_activities(key)
        if cached:
            return cached

    # 3. Попытка обращения к BindingDB REST API (по UniProt P00918 / CA II человека) с таймаутом и ретраями
    api_url = (
        "https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot"
        f"?uniprot=P00918&cutoff=10000&code=0&response=application/json"
    )
    live_data = _get_json_with_retry(api_url, timeout=timeout, max_retries=3)
    if isinstance(live_data, dict) and live_data.get("affinities"):
        fetched_at = _now_iso()
        matched = [
            row
            for row in live_data.get("affinities", [])
            if key in str(row.get("monomerid") or row.get("smiles") or "").lower()
        ]
        if matched:
            record = {
                "ligand": ligand_name,
                "activities": matched,
                "source": "bindingdb_live",
                "_from_cache": False,
                "_from_download": False,
                "is_cached": False,
                "fetched_at": fetched_at,
            }
            save_cached_activities(key, record)
            return record

    return {
        "ligand": ligand_name,
        "activities": [],
        "source": "bindingdb",
        "note": "Нет локального кеша; данные из live API BindingDB недоступны.",
        "_from_cache": False,
        "_from_download": False,
        "is_cached": False,
        "fetched_at": _now_iso(),
    }
