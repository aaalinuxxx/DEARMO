"""PubChem PUG REST: универсальный поиск (имя / SMILES POST / CID), поиск аналогов, обработка неоднозначности и кэш с датой."""
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

from tools.rdkit_tools import parse_smiles, tanimoto_similarity

_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_REPO_ROOT = Path(__file__).resolve().parent.parent
_PUBCHEM_FILE = _REPO_ROOT / "data" / "external" / "pubchem" / "acetazolamide.json"
_CACHE_DIR = _REPO_ROOT / "data" / "cache" / "pubchem"

_FALLBACK_SMILES = {
    "acetazolamide": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    "ацетазоламид": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    "methazolamide": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
    "метазоламид": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
    "ethoxzolamide": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
    "этоксзоламид": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
    "benzamidine": "NC(=N)c1ccccc1",
    "бензамидин": "NC(=N)c1ccccc1",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "аспирин": "CC(=O)Oc1ccccc1C(=O)O",
}

_FALLBACK_CIDS = {
    "acetazolamide": (1986, "Acetazolamide", 222.25),
    "ацетазоламид": (1986, "Acetazolamide", 222.25),
    "methazolamide": (4100, "Methazolamide", 236.28),
    "метазоламид": (4100, "Methazolamide", 236.28),
    "ethoxzolamide": (3295, "Ethoxzolamide", 258.32),
    "этоксзоламид": (3295, "Ethoxzolamide", 258.32),
    "benzamidine": (2332, "Benzamidine", 120.15),
    "бензамидин": (2332, "Benzamidine", 120.15),
    "aspirin": (2244, "Aspirin", 180.16),
    "аспирин": (2244, "Aspirin", 180.16),
}

# Словарь заведомо неоднозначных общих названий/классов соединений
_AMBIGUOUS_NAMES: dict[str, list[dict[str, Any]]] = {
    "сульфонамид": [
        {"cid": 1986, "title": "Acetazolamide", "connectivity_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1"},
        {"cid": 4100, "title": "Methazolamide", "connectivity_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C"},
        {"cid": 5333, "title": "Sulfanilamide", "connectivity_smiles": "Nc1ccc(S(N)(=O)=O)cc1"},
    ],
    "sulfonamide": [
        {"cid": 1986, "title": "Acetazolamide", "connectivity_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1"},
        {"cid": 4100, "title": "Methazolamide", "connectivity_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C"},
        {"cid": 5333, "title": "Sulfanilamide", "connectivity_smiles": "Nc1ccc(S(N)(=O)=O)cc1"},
    ],
    "ингибитор": [
        {"cid": 1986, "title": "Acetazolamide", "connectivity_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1"},
        {"cid": 3154, "title": "Dorzolamide", "connectivity_smiles": "CCNC1CC(C)S(=O)(=O)c2sc(S(N)(=O)=O)cc12"},
    ],
    "inhibitor": [
        {"cid": 1986, "title": "Acetazolamide", "connectivity_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1"},
        {"cid": 3154, "title": "Dorzolamide", "connectivity_smiles": "CCNC1CC(C)S(=O)(=O)c2sc(S(N)(=O)=O)cc12"},
    ],
}

# Оффлайн-каталог проверенных аналогов сульфонамидных ингибиторов CA II для демо и работы без сети
_OFFLINE_ANALOGS: list[dict[str, Any]] = [
    {
        "cid": 4100,
        "title": "Methazolamide",
        "connectivity_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
        "molecular_weight": 236.28,
    },
    {
        "cid": 3295,
        "title": "Ethoxzolamide",
        "connectivity_smiles": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
        "molecular_weight": 258.32,
    },
    {
        "cid": 3154,
        "title": "Dorzolamide",
        "connectivity_smiles": "CCNC1CC(C)S(=O)(=O)c2sc(S(N)(=O)=O)cc12",
        "molecular_weight": 324.44,
    },
    {
        "cid": 68844,
        "title": "Brinzolamide",
        "connectivity_smiles": "CCNC1CN(CCCOC)S(=O)(=O)c2sc(S(N)(=O)=O)cc12",
        "molecular_weight": 383.51,
    },
    {
        "cid": 3038,
        "title": "Dichlorphenamide",
        "connectivity_smiles": "NS(=O)(=O)c1cc(S(N)(=O)=O)c(Cl)c(Cl)c1",
        "molecular_weight": 305.16,
    },
]


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cache_file(prefix: str, key: str) -> Path:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
    safe = re.sub(r"[^a-zA-Z0-9_-]", "_", key)[:32]
    return _CACHE_DIR / f"{prefix}_{safe}_{digest}.json"


def _read_cache(prefix: str, key: str) -> dict[str, Any] | None:
    path = _cache_file(prefix, key)
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            wrapper = json.load(f)
        if isinstance(wrapper, dict) and "payload" in wrapper:
            data = dict(wrapper["payload"])
            data["is_cached"] = True
            data["fetched_at"] = wrapper.get("fetched_at") or "unknown"
            return data
    except (OSError, json.JSONDecodeError):
        return None
    return None


def _write_cache(prefix: str, key: str, payload: dict[str, Any], fetched_at: str | None = None) -> str:
    ts = fetched_at or _now_iso()
    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path = _cache_file(prefix, key)
        wrapper = {
            "fetched_at": ts,
            "source": "pubchem",
            "query": key,
            "payload": payload,
        }
        with path.open("w", encoding="utf-8") as f:
            json.dump(wrapper, f, ensure_ascii=False, indent=2)
    except OSError:
        pass
    return ts


def _request_json(
    url: str,
    *,
    method: str = "GET",
    form_data: dict[str, str] | None = None,
    timeout: float = 10.0,
    max_retries: int = 3,
) -> dict[str, Any] | None:
    """Выполнить HTTP-запрос к PubChem с таймаутом и повторами на 429/5xx."""
    headers = {"Accept": "application/json"}
    body: bytes | None = None
    if form_data is not None:
        body = urllib.parse.urlencode(form_data).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"

    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, data=body, headers=headers, method=method)
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


def _get_json(url: str, timeout: float = 10.0) -> dict[str, Any] | None:
    """Совместимая обёртка для GET-запросов."""
    return _request_json(url, method="GET", timeout=timeout)


def _is_smiles_query(text: str) -> bool:
    """Проверить, является ли строка валидным SMILES (а не обычным словом/названием)."""
    # Если это обычное слово (только буквы и дефис, длина >= 4, и нет явных SMILES-элементов вроде цифр/скобок/=/#),
    # считаем его названием, кроме коротких линейных цепочек вида CCCO/CCCC.
    if re.match(r"^[A-Za-zА-Яа-яЁё\-]{4,}$", text) and not re.match(r"^[CNOSPFIBrcnosp]+$", text):
        return False
    if re.match(r"^[A-Z][a-z]{3,}$", text):
        return False
    try:
        parse_smiles(text)
        return True
    except ValueError:
        return False


def _extract_properties_list(data: dict[str, Any]) -> list[dict[str, Any]]:
    props = data.get("PropertyTable", {}).get("Properties", [])
    return props if isinstance(props, list) else []


def resolve_compound(name_or_smiles: str, *, use_cache: bool = True) -> dict[str, Any]:
    """Разрешить соединение по имени, SMILES (POST) или CID с проверкой неоднозначности и кэшем."""
    if not isinstance(name_or_smiles, str) or not name_or_smiles.strip():
        return {"query": "", "source": "pubchem", "error": "empty_query", "is_cached": False}

    key = name_or_smiles.strip()
    lower = key.lower()

    # 1. Проверка заведомо неоднозначных названий
    if lower in _AMBIGUOUS_NAMES:
        candidates = _AMBIGUOUS_NAMES[lower]
        return {
            "query": key,
            "ambiguous": True,
            "status": "needs_clarification",
            "candidates": candidates,
            "error": "ambiguous_query",
            "message": f"Название '{key}' неоднозначно и соответствует нескольким соединениям. Уточните CID или SMILES.",
            "source": "pubchem",
            "is_cached": False,
        }

    # 2. Проверка локального файла предзагрузки для ацетазоламида
    if use_cache and lower in ("acetazolamide", "ацетазоламид") and _PUBCHEM_FILE.is_file():
        try:
            mtime = datetime.fromtimestamp(_PUBCHEM_FILE.stat().st_mtime, tz=timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            )
            with _PUBCHEM_FILE.open(encoding="utf-8") as f:
                payload = json.load(f)
            props = _extract_properties_list(payload)
            first = props[0] if props else {}
            return {
                "query": key,
                "cid": first.get("CID"),
                "title": first.get("Title"),
                "connectivity_smiles": first.get("ConnectivitySMILES") or first.get("CanonicalSMILES"),
                "canonical_smiles": first.get("CanonicalSMILES") or first.get("ConnectivitySMILES"),
                "inchikey": first.get("InChIKey"),
                "molecular_weight": first.get("MolecularWeight"),
                "source": "pubchem_download",
                "is_cached": True,
                "fetched_at": mtime,
                "ambiguous": False,
            }
        except (OSError, json.JSONDecodeError):
            pass

    # 3. Проверка файлового кэша с датой
    if use_cache:
        cached = _read_cache("resolve", lower)
        if cached is not None:
            return cached

    # 4. Определение способа запроса к PubChem PUG REST (CID / SMILES POST / Name GET)
    prop_fields = "Title,CID,ConnectivitySMILES,CanonicalSMILES,IsomericSMILES,InChIKey,MolecularWeight"
    cid_match = re.match(r"^(?:cid[:\s]*)?(\d+)$", lower)

    data: dict[str, Any] | None = None
    query_mode = "name"
    if cid_match:
        query_mode = "cid"
        cid_num = cid_match.group(1)
        url = f"{_PUG}/compound/cid/{cid_num}/property/{prop_fields}/JSON"
        data = _get_json(url)
    elif _is_smiles_query(key):
        query_mode = "smiles"
        url = f"{_PUG}/compound/smiles/property/{prop_fields}/JSON"
        data = _request_json(url, method="POST", form_data={"smiles": key})
    else:
        url = f"{_PUG}/compound/name/{urllib.parse.quote(key)}/property/{prop_fields}/JSON"
        data = _get_json(url)

    # 5. Разбор ответа API (включая проверку на несколько уникальных CID при поиске по имени)
    if data:
        props_list = _extract_properties_list(data)
        unique_by_cid: dict[Any, dict[str, Any]] = {}
        for item in props_list:
            cid_val = item.get("CID")
            if cid_val is not None and cid_val not in unique_by_cid:
                unique_by_cid[cid_val] = item

        if query_mode == "name" and len(unique_by_cid) > 1:
            # Проверим, отличаются ли их канонические SMILES (разные структуры)
            unique_smiles = {
                p.get("ConnectivitySMILES") or p.get("CanonicalSMILES")
                for p in unique_by_cid.values()
                if p.get("ConnectivitySMILES") or p.get("CanonicalSMILES")
            }
            if len(unique_smiles) > 1:
                candidates = [
                    {
                        "cid": p.get("CID"),
                        "title": p.get("Title"),
                        "connectivity_smiles": p.get("ConnectivitySMILES") or p.get("CanonicalSMILES"),
                        "molecular_weight": p.get("MolecularWeight"),
                    }
                    for p in list(unique_by_cid.values())[:10]
                ]
                return {
                    "query": key,
                    "ambiguous": True,
                    "status": "needs_clarification",
                    "candidates": candidates,
                    "error": "ambiguous_query",
                    "message": (
                        f"По запросу '{key}' найдено несколько различных структур ({len(candidates)} кандидатов). "
                        "Укажите точный CID или SMILES."
                    ),
                    "source": "pubchem",
                    "is_cached": False,
                    "fetched_at": _now_iso(),
                }

        if props_list:
            first = props_list[0]
            smiles_val = (
                first.get("ConnectivitySMILES")
                or first.get("CanonicalSMILES")
                or first.get("IsomericSMILES")
            )
            fetched_at = _now_iso()
            result = {
                "query": key,
                "query_mode": query_mode,
                "cid": first.get("CID"),
                "title": first.get("Title"),
                "connectivity_smiles": smiles_val,
                "canonical_smiles": first.get("CanonicalSMILES") or smiles_val,
                "isomeric_smiles": first.get("IsomericSMILES"),
                "inchikey": first.get("InChIKey"),
                "molecular_weight": first.get("MolecularWeight"),
                "source": "pubchem",
                "is_cached": False,
                "fetched_at": fetched_at,
                "ambiguous": False,
            }
            _write_cache("resolve", lower, result, fetched_at=fetched_at)
            return result

    # 6. Откат на локальный справочник для известных соединений или если передан валидный SMILES без сети
    if lower in _FALLBACK_SMILES:
        cid_val, title_val, mw_val = _FALLBACK_CIDS.get(lower, (1986, key.title(), None))
        return {
            "query": key,
            "query_mode": query_mode,
            "cid": cid_val,
            "title": title_val,
            "connectivity_smiles": _FALLBACK_SMILES[lower],
            "canonical_smiles": _FALLBACK_SMILES[lower],
            "molecular_weight": mw_val,
            "source": "pubchem_fallback",
            "is_cached": True,
            "fetched_at": "offline-fallback",
            "ambiguous": False,
        }

    if query_mode == "smiles":
        return {
            "query": key,
            "query_mode": "smiles",
            "cid": None,
            "title": None,
            "connectivity_smiles": key,
            "canonical_smiles": key,
            "source": "smiles_input_offline",
            "is_cached": False,
            "fetched_at": _now_iso(),
            "ambiguous": False,
        }

    return {
        "query": key,
        "source": "pubchem",
        "error": "no_data",
        "is_cached": False,
        "fetched_at": _now_iso(),
        "ambiguous": False,
    }


def search_similar(
    smiles: str,
    *,
    threshold: int = 80,
    max_records: int = 10,
    use_cache: bool = True,
) -> dict[str, Any]:
    """Поиск похожих соединений (2D Tanimoto similarity) через PubChem fastsimilarity_2d с кэшем и расчётом Tanimoto."""
    try:
        parse_smiles(smiles)
    except ValueError as exc:
        return {
            "query_smiles": smiles,
            "status": "error",
            "error": str(exc),
            "analogs": [],
            "source": "pubchem",
            "is_cached": False,
        }

    clean_smiles = smiles.strip()
    cache_key = f"{clean_smiles}_th{threshold}_max{max_records}"
    if use_cache:
        cached = _read_cache("similar", cache_key)
        if cached is not None:
            return cached

    url_cids = (
        f"{_PUG}/compound/fastsimilarity_2d/smiles/cids/JSON"
        f"?Threshold={int(threshold)}&MaxRecords={int(max_records)}"
    )
    cid_data = _request_json(url_cids, method="POST", form_data={"smiles": clean_smiles})
    cids: list[int] = []
    if cid_data:
        cids = cid_data.get("IdentifierList", {}).get("CID", []) or []

    analogs: list[dict[str, Any]] = []
    if cids:
        cids_str = ",".join(str(c) for c in cids[:max_records])
        prop_url = (
            f"{_PUG}/compound/cid/{cids_str}"
            "/property/Title,CID,ConnectivitySMILES,CanonicalSMILES,MolecularWeight/JSON"
        )
        prop_data = _get_json(prop_url)
        if prop_data:
            for item in _extract_properties_list(prop_data):
                cand_smiles = item.get("ConnectivitySMILES") or item.get("CanonicalSMILES")
                if not cand_smiles:
                    continue
                try:
                    sim = tanimoto_similarity(clean_smiles, cand_smiles)
                except ValueError:
                    sim = 0.0
                analogs.append(
                    {
                        "cid": item.get("CID"),
                        "title": item.get("Title"),
                        "connectivity_smiles": cand_smiles,
                        "molecular_weight": item.get("MolecularWeight"),
                        "tanimoto_similarity": sim,
                        "source": "pubchem",
                    }
                )

    if analogs:
        analogs.sort(key=lambda x: x.get("tanimoto_similarity", 0.0), reverse=True)
        fetched_at = _now_iso()
        result = {
            "query_smiles": clean_smiles,
            "threshold": threshold,
            "status": "success",
            "analogs": analogs[:max_records],
            "source": "pubchem",
            "is_cached": False,
            "fetched_at": fetched_at,
        }
        _write_cache("similar", cache_key, result, fetched_at=fetched_at)
        return result

    # Откат на проверенную оффлайн-выборку аналогов при недоступности сети
    fallback_analogs: list[dict[str, Any]] = []
    for item in _OFFLINE_ANALOGS:
        cand_smiles = item["connectivity_smiles"]
        if cand_smiles == clean_smiles:
            continue
        try:
            sim = tanimoto_similarity(clean_smiles, cand_smiles)
        except ValueError:
            sim = 0.0
        fallback_analogs.append(
            {
                **item,
                "tanimoto_similarity": sim,
                "source": "pubchem_offline_catalog",
            }
        )
    fallback_analogs.sort(key=lambda x: x.get("tanimoto_similarity", 0.0), reverse=True)
    return {
        "query_smiles": clean_smiles,
        "threshold": threshold,
        "status": "success",
        "analogs": fallback_analogs[:max_records],
        "source": "pubchem_offline_catalog",
        "is_cached": True,
        "fetched_at": "offline-fallback",
    }
