"""BindingDB: чтение кэша с метаданными даты (fetched_at), устойчивые HTTP-запросы (timeout + retry) и выбор мишени."""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "bindingdb"
_EXTERNAL_DIR = Path(__file__).resolve().parent.parent / "data" / "external" / "bindingdb"

def _cache_key(ligand_name: str, smiles: str | None = None) -> str:
    """Ключ кэша без словаря известных молекул.

    Для латинского названия это его слаг, иначе — хэш SMILES: он одинаков для
    «ацетазоламид», «acetazolamide» и самой структуры, и не требует списка синонимов.
    """
    raw = (ligand_name or "").strip().lower().replace(" ", "_")
    ascii_slug = re.sub(r"[^a-z0-9_\-]", "", raw)
    if ascii_slug and ascii_slug == raw:
        return ascii_slug
    if smiles:
        return "smi_" + hashlib.sha1(smiles.strip().encode("utf-8")).hexdigest()[:12]
    return ascii_slug or "unknown"


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


def _parse_affinity(raw: Any) -> tuple[str, float | None]:
    text = str(raw).strip().replace(",", "")
    relation = "="
    if text[:1] in (">", "<"):
        relation = text[0]
        text = text[1:].strip()
    try:
        return relation, float(text)
    except ValueError:
        return relation, None


def _activities_from_payload(payload: Any) -> list[dict[str, Any]]:
    """Разобрать ответ getTargetByCompound: список есть и под разными обёртками."""
    buckets: list[list] = []

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key in ("bdb.affinities", "affinities"):
                if isinstance(node.get(key), list):
                    buckets.append(node[key])
                    return
            for value in node.values():
                walk(value)
        elif isinstance(node, list) and node and isinstance(node[0], dict):
            buckets.append(node)

    walk(payload)
    activities: list[dict[str, Any]] = []
    for bucket in buckets:
        for row in bucket:
            if not isinstance(row, dict):
                continue
            relation, value = _parse_affinity(row.get("bdb.affinity") or row.get("affinity"))
            if value is None:
                continue
            monomer_id = row.get("bdb.monomerid") or row.get("monomerid")
            activities.append(
                {
                    "target_name": row.get("bdb.target") or row.get("target"),
                    "organism": row.get("bdb.species") or row.get("species"),
                    "affinity_type": row.get("bdb.affinity_type") or row.get("affinity_type") or "Ki",
                    "value": value,
                    "relation": relation,
                    "units": "nM",
                    "monomer_id": monomer_id,
                    "link": record_url(monomer_id),
                    "pdb_id": row.get("pdb_id"),
                }
            )
    return activities


def record_url(monomer_id: Any) -> str | None:
    """Ссылка на саму запись BindingDB.

    REST-ответ не содержит ни PMID, ни ссылки, поэтому единственный способ дать
    пользователю проверяемый источник измерения — собрать адрес записи по её
    идентификатору. Статья приходит отдельно, из цитирования структуры PDB.
    """
    if monomer_id in (None, "", 0):
        return None
    return (
        "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/"
        f"MolStructure.jsp?monomerid={monomer_id}"
    )


def _looks_like_smiles(text: str) -> bool:
    return any(ch in text for ch in "=#()[]")


def recommend_target(record: dict[str, Any]) -> dict[str, Any] | None:
    """Мишень выводится из измерений этого соединения — одинаково для любой молекулы."""
    return _recommend_from_measurements(record)


def _recommend_from_measurements(record: dict[str, Any]) -> dict[str, Any] | None:
    """Лучшая измеренная строка: человек, затем Ki/Kd, затем меньшее значение. PDB не выдумывается."""
    ranked: list[tuple[tuple, dict[str, Any]]] = []
    for row in record.get("activities") or []:
        if row.get("relation") not in (None, "=", ""):
            continue
        try:
            value = float(row["value"])
        except (TypeError, ValueError, KeyError):
            continue
        organism = (row.get("organism") or "").lower()
        human = 0 if "human" in organism or "homo sapiens" in organism else 1
        kind = str(row.get("affinity_type") or "").upper()
        potency = 0 if kind in {"KI", "KD"} else 1
        ranked.append(((human, potency, value), row))
    if not ranked:
        return None
    _, chosen = sorted(ranked, key=lambda item: item[0])[0]
    return {
        "name": chosen.get("target_name"),
        "organism": chosen.get("organism"),
        "pdb_id": chosen.get("pdb_id"),
        "hypothesis": False,
        "affinity_type": chosen.get("affinity_type"),
        "value": chosen.get("value"),
        "units": chosen.get("units"),
        "measured_target_name": chosen.get("target_name"),
        "link": chosen.get("link") or record_url(chosen.get("monomer_id")),
        "rationale": (
            "Мишень взята из строки BindingDB этого соединения: сначала человек, "
            "затем Ki или Kd, затем меньшее значение. Идентификатор PDB из другой молекулы не подставляется."
        ),
    }


def _ensure_links(record: dict[str, Any]) -> dict[str, Any]:
    """Достроить ссылки на записи у строк, загруженных до их появления.

    Предзагруженный набор в data/external и ранние файлы кэша писались без поля
    link, и из-за этого «ссылка на запись» — прямое требование кейса — пропадала
    именно в демонстрационном прогоне. Идентификатор записи в этих файлах есть,
    так что адрес восстанавливается без обращения к сети.
    """
    for row in record.get("activities") or []:
        if not row.get("link"):
            row["link"] = record_url(row.get("monomer_id"))
    recommended = record.get("recommended_target")
    if isinstance(recommended, dict) and not recommended.get("link"):
        recommended["link"] = record_url(recommended.get("monomer_id"))
    return record


def fetch_activities(
    ligand_name: str,
    *,
    smiles: str | None = None,
    use_cache: bool = True,
    timeout: float = 12.0,
) -> dict[str, Any]:
    """Получить записи активности из BindingDB (локальная предзагрузка, файловый кэш или REST API с повторами)."""
    key = _cache_key(ligand_name, smiles)

    # 1. Предзагруженный эталонный файл репозитория (сохраняем обратную совместимость _from_download=True)
    preloaded = _EXTERNAL_DIR / f"{key}.json"
    if use_cache and preloaded.is_file():
        with preloaded.open(encoding="utf-8") as f:
            record = json.load(f)
        fetched_at = record.get("fetched_at") or _file_mtime_iso(preloaded)
        record["_from_cache"] = False
        record["_from_download"] = True
        record["is_cached"] = True
        record["fetched_at"] = fetched_at
        rec = recommend_target(_ensure_links(record))
        if rec:
            record["recommended_target"] = rec
        return _ensure_links(record)

    # 2. Локальный кэш в data/cache/bindingdb/
    if use_cache:
        cached = load_cached_activities(key)
        if cached:
            if not cached.get("recommended_target"):
                rec = recommend_target(_ensure_links(cached))
                if rec:
                    cached["recommended_target"] = rec
            return _ensure_links(cached)

    # 3. BindingDB по SMILES самого соединения, без привязки к заранее выбранной мишени.
    query_smiles = (smiles or "").strip()
    if not query_smiles and _looks_like_smiles(ligand_name):
        query_smiles = ligand_name.strip()
    if query_smiles:
        api_url = (
            "https://bindingdb.org/rest/getTargetByCompound"
            f"?smiles={urllib.parse.quote(query_smiles, safe='')}&cutoff=100000&response=application/json"
        )
        live_data = _get_json_with_retry(api_url, timeout=timeout, max_retries=2)
        activities = _activities_from_payload(live_data) if live_data else []
        if activities:
            record = {
                "ligand": ligand_name,
                "smiles": query_smiles,
                "activities": activities,
                "source": "bindingdb_live",
                "_from_cache": False,
                "_from_download": False,
                "is_cached": False,
                "fetched_at": _now_iso(),
            }
            rec = recommend_target(record)
            if rec:
                record["recommended_target"] = rec
            save_cached_activities(key, record)
            return _ensure_links(record)

    return {
        "ligand": ligand_name,
        "activities": [],
        "source": "bindingdb",
        "note": (
            "Нет локального кеша BindingDB. Живой запрос выполняется по SMILES соединения "
            "(getTargetByCompound) и сейчас ничего не вернул."
            if query_smiles
            else "Нет локального кеша BindingDB и нет SMILES, чтобы запросить активности соединения."
        ),
        "_from_cache": False,
        "_from_download": False,
        "is_cached": False,
        "fetched_at": _now_iso(),
    }
