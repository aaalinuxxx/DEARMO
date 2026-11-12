"""Получение данных биологической активности из ChEMBL с приведением единиц в нМ, удалением дублей и кэшем."""
from __future__ import annotations

import hashlib
import json
import re
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools.rdkit_tools import parse_smiles, standardize_smiles

_CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
_REPO_ROOT = Path(__file__).resolve().parent.parent
_CACHE_DIR = _REPO_ROOT / "data" / "cache" / "chembl"

# Коэффициенты перевода концентрационных единиц в наномоли (nM)
_UNIT_TO_NM_FACTOR: dict[str, float] = {
    "nm": 1.0,
    "nmol/l": 1.0,
    "nmolar": 1.0,
    "10'-9m": 1.0,
    "um": 1_000.0,
    "µm": 1_000.0,
    "μm": 1_000.0,
    "umol/l": 1_000.0,
    "10'-6m": 1_000.0,
    "mm": 1_000_000.0,
    "mmol/l": 1_000_000.0,
    "10'-3m": 1_000_000.0,
    "m": 1_000_000_000.0,
    "mol/l": 1_000_000_000.0,
    "pm": 0.001,
    "pmol/l": 0.001,
    "10'-12m": 0.001,
    "fm": 1e-6,
    "fmol/l": 1e-6,
}

# Синонимы известных молекул и мишеней для надёжного разрешения в оффлайн- и онлайн-режимах
_MOLECULE_ALIASES: dict[str, dict[str, Any]] = {
    "acetazolamide": {
        "molecule_chembl_id": "CHEMBL20",
        "pref_name": "ACETAZOLAMIDE",
        "canonical_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    },
    "ацетазоламид": {
        "molecule_chembl_id": "CHEMBL20",
        "pref_name": "ACETAZOLAMIDE",
        "canonical_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    },
    "chembl20": {
        "molecule_chembl_id": "CHEMBL20",
        "pref_name": "ACETAZOLAMIDE",
        "canonical_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
    },
    "methazolamide": {
        "molecule_chembl_id": "CHEMBL422",
        "pref_name": "METHAZOLAMIDE",
        "canonical_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
    },
    "метазоламид": {
        "molecule_chembl_id": "CHEMBL422",
        "pref_name": "METHAZOLAMIDE",
        "canonical_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
    },
    "ethoxzolamide": {
        "molecule_chembl_id": "CHEMBL533",
        "pref_name": "ETHOXZOLAMIDE",
        "canonical_smiles": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
    },
    "dorzolamide": {
        "molecule_chembl_id": "CHEMBL278980",
        "pref_name": "DORZOLAMIDE",
        "canonical_smiles": "CCNC1CC(C)S(=O)(=O)c2sc(S(N)(=O)=O)cc12",
    },
    "benzamidine": {
        "molecule_chembl_id": "CHEMBL12131",
        "pref_name": "BENZAMIDINE",
        "canonical_smiles": "NC(=N)c1ccccc1",
    },
}

_TARGET_ALIASES: dict[str, str] = {
    "chembl205": "CHEMBL205",
    "ca2": "CHEMBL205",
    "ca ii": "CHEMBL205",
    "carbonic anhydrase 2": "CHEMBL205",
    "carbonic anhydrase": "CHEMBL205",
    "carbonic anhydrase ii": "CHEMBL205",
    "карбоангидраза ii": "CHEMBL205",
    "карбоангидраза 2": "CHEMBL205",
    "trypsin": "CHEMBL209",
    "трипсин": "CHEMBL209",
}

# Эталонные данные ChEMBL по карбоангидразе II человека (CHEMBL205) для оффлайн-работы
_OFFLINE_ACTIVITIES: dict[tuple[str, str], list[dict[str, Any]]] = {
    ("CHEMBL205", "CHEMBL20"): [
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_pref_name": "Carbonic anhydrase II",
            "target_organism": "Homo sapiens",
            "assay_chembl_id": "CHEMBL658249",
            "assay_type": "B",
            "assay_description": "Inhibition of human recombinant Carbonic anhydrase II (CO2 hydration assay)",
            "standard_type": "Ki",
            "standard_relation": "=",
            "standard_value": "12.0",
            "standard_units": "nM",
            "pchembl_value": "7.92",
            "document_chembl_id": "CHEMBL1145821",
            "data_validity_comment": None,
        },
        {
            "molecule_chembl_id": "CHEMBL20",
            "target_chembl_id": "CHEMBL205",
            "target_pref_name": "Carbonic anhydrase II",
            "target_organism": "Homo sapiens",
            "assay_chembl_id": "CHEMBL872914",
            "assay_type": "B",
            "assay_description": "Inhibition of human Carbonic anhydrase II esterase activity (4-nitrophenyl acetate)",
            "standard_type": "IC50",
            "standard_relation": "=",
            "standard_value": "0.018",
            "standard_units": "uM",
            "pchembl_value": "7.74",
            "document_chembl_id": "CHEMBL1152904",
            "data_validity_comment": None,
        },
    ],
    ("CHEMBL205", "CHEMBL422"): [
        {
            "molecule_chembl_id": "CHEMBL422",
            "target_chembl_id": "CHEMBL205",
            "target_pref_name": "Carbonic anhydrase II",
            "target_organism": "Homo sapiens",
            "assay_chembl_id": "CHEMBL658250",
            "assay_type": "B",
            "assay_description": "Inhibition of human recombinant Carbonic anhydrase II",
            "standard_type": "Ki",
            "standard_relation": "=",
            "standard_value": "14.0",
            "standard_units": "nM",
            "pchembl_value": "7.85",
            "document_chembl_id": "CHEMBL1145821",
            "data_validity_comment": None,
        },
    ],
    ("CHEMBL205", "CHEMBL533"): [
        {
            "molecule_chembl_id": "CHEMBL533",
            "target_chembl_id": "CHEMBL205",
            "target_pref_name": "Carbonic anhydrase II",
            "target_organism": "Homo sapiens",
            "assay_chembl_id": "CHEMBL658251",
            "assay_type": "B",
            "assay_description": "Inhibition of human recombinant Carbonic anhydrase II",
            "standard_type": "Ki",
            "standard_relation": "=",
            "standard_value": "8.0",
            "standard_units": "nM",
            "pchembl_value": "8.10",
            "document_chembl_id": "CHEMBL1145821",
            "data_validity_comment": None,
        },
    ],
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cache_file(prefix: str, key: str) -> Path:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
    safe = re.sub(r"[^a-zA-Z0-9_-]", "_", key)[:32]
    return _CACHE_DIR / f"{prefix}_{safe}_{digest}.json"


def _read_cache(prefix: str, key: str) -> tuple[Any, str] | None:
    path = _cache_file(prefix, key)
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            wrapper = json.load(f)
        if isinstance(wrapper, dict) and "payload" in wrapper:
            return wrapper["payload"], wrapper.get("fetched_at") or "unknown"
    except (OSError, json.JSONDecodeError):
        return None
    return None


def _write_cache(prefix: str, key: str, payload: Any, fetched_at: str | None = None) -> str:
    ts = fetched_at or _now_iso()
    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path = _cache_file(prefix, key)
        wrapper = {
            "fetched_at": ts,
            "source": "chembl",
            "query": key,
            "payload": payload,
        }
        with path.open("w", encoding="utf-8") as f:
            json.dump(wrapper, f, ensure_ascii=False, indent=2)
    except OSError:
        pass
    return ts


def _get_json(url: str, *, timeout: float = 3.0, max_retries: int = 1) -> dict[str, Any] | None:
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


def normalize_to_nm(
    value: float | int | str | None,
    unit: str | None,
) -> tuple[float | None, str | None, bool]:
    """Привести значение активности к наномолям (nM).

    - Пропуски (None / пустая строка) сохраняются как (None, unit, False), НИКОГДА не подменяются на 0.
    - Концентрационные единицы (nM, uM, mM, M, pM) переводятся в nM с флагом comparable=True.
    - Несопоставимые единицы (%, ug/mL, mg/kg и т.д.) не переводятся, возвращают comparable=False.
    """
    if value is None:
        return None, unit, False
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None, unit, False
        try:
            num_val = float(stripped)
        except ValueError:
            return None, unit, False
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        num_val = float(value)
    else:
        return None, unit, False

    if unit is None or not str(unit).strip():
        return num_val, None, False

    clean_unit = str(unit).strip().lower().replace(" ", "")
    factor = _UNIT_TO_NM_FACTOR.get(clean_unit)
    if factor is not None:
        return round(num_val * factor, 6), "nM", True

    return num_val, str(unit).strip(), False


def _normalize_record(
    raw: dict[str, Any],
    *,
    is_cached: bool = False,
    fetched_at: str | None = None,
    source: str = "chembl",
) -> dict[str, Any]:
    """Привести сырую запись ChEMBL к стандартному контракту с контекстом эксперимента."""
    std_val = raw.get("standard_value")
    if std_val is None and "value" in raw:
        std_val = raw.get("value")
    std_unit = raw.get("standard_units") or raw.get("units")
    std_type = raw.get("standard_type") or raw.get("affinity_type") or "Activity"
    std_rel = raw.get("standard_relation") or raw.get("relation") or "="

    val_nm, norm_unit, unit_comparable = normalize_to_nm(std_val, std_unit)
    type_comparable = str(std_type).upper() in {"IC50", "KI", "KD", "EC50", "KB"}
    validity_comment = raw.get("data_validity_comment")
    has_validity_issue = bool(validity_comment and str(validity_comment).strip())

    comparable = bool(unit_comparable and type_comparable and val_nm is not None and not has_validity_issue)

    raw_float: float | None = None
    if std_val is not None:
        try:
            raw_float = float(std_val)
        except (ValueError, TypeError):
            raw_float = None

    pchembl = raw.get("pchembl_value")
    if pchembl is not None:
        try:
            pchembl = float(pchembl)
        except (ValueError, TypeError):
            pchembl = None

    return {
        "molecule_chembl_id": raw.get("molecule_chembl_id"),
        "target_chembl_id": raw.get("target_chembl_id"),
        "target_pref_name": raw.get("target_pref_name") or raw.get("target_name"),
        "target_organism": raw.get("target_organism") or raw.get("organism"),
        "assay_chembl_id": raw.get("assay_chembl_id"),
        "assay_type": raw.get("assay_type"),
        "assay_description": raw.get("assay_description"),
        "standard_type": std_type,
        "standard_relation": std_rel,
        "standard_value": raw_float,
        "standard_units": std_unit,
        "value_nm": val_nm,
        "normalized_units": norm_unit,
        "comparable": comparable,
        "pchembl_value": pchembl,
        "document_chembl_id": raw.get("document_chembl_id"),
        "data_validity_comment": validity_comment,
        "conflicting_measurements": False,
        "source": source,
        "is_cached": is_cached,
        "fetched_at": fetched_at or _now_iso(),
    }


def deduplicate_activities(
    records: list[dict[str, Any]],
    *,
    divergence_fold: float = 10.0,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Удалить полные дубликаты и обработать повторные измерения по правилам биотехнологов."""
    deduped: list[dict[str, Any]] = []
    limitations: list[str] = []
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}

    for rec in records:
        rec_copy = dict(rec)
        rec_copy.setdefault("conflicting_measurements", False)
        if not rec_copy.get("comparable") or rec_copy.get("value_nm") is None:
            deduped.append(rec_copy)
            continue
        key = (
            rec_copy.get("molecule_chembl_id"),
            rec_copy.get("target_chembl_id"),
            rec_copy.get("target_organism"),
            str(rec_copy.get("standard_type") or "").upper(),
            rec_copy.get("assay_chembl_id"),
            rec_copy.get("assay_type"),
            rec_copy.get("standard_relation") or "=",
        )
        grouped.setdefault(key, []).append(rec_copy)

    for key, group in grouped.items():
        if len(group) == 1:
            group[0]["conflicting_measurements"] = False
            deduped.append(group[0])
            continue

        vals = [float(g["value_nm"]) for g in group if g.get("value_nm") is not None and float(g["value_nm"]) > 0]
        if len(vals) >= 2 and (max(vals) / min(vals)) >= divergence_fold:
            fold = round(max(vals) / min(vals), 1)
            mol_id, tgt_id, org, stype, assay_id, _, _ = key
            limitations.append(
                f"Обнаружены противоречащие измерения {stype} для {mol_id} на мишени {tgt_id} "
                f"(assay {assay_id}, разброс в {fold} раз: {min(vals)} vs {max(vals)} nM). "
                "Значения сохранены раздельно без усреднения."
            )
            for g in group:
                g["conflicting_measurements"] = True
                g["divergence_fold"] = fold
                deduped.append(g)
        else:
            median_nm = round(float(statistics.median(vals)), 6) if vals else group[0].get("value_nm")
            merged = dict(group[0])
            merged["value_nm"] = median_nm
            merged["raw_values_nm"] = vals
            merged["replicates_count"] = len(group)
            merged["conflicting_measurements"] = False
            docs = sorted({g.get("document_chembl_id") for g in group if g.get("document_chembl_id")})
            if docs:
                merged["source_documents"] = docs
            deduped.append(merged)

    by_target_type: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for rec in deduped:
        if rec.get("comparable") and rec.get("value_nm") is not None and float(rec["value_nm"]) > 0:
            tkey = (
                rec.get("molecule_chembl_id"),
                rec.get("target_chembl_id"),
                rec.get("target_organism"),
                str(rec.get("standard_type") or "").upper(),
            )
            by_target_type.setdefault(tkey, []).append(rec)

    for tkey, items in by_target_type.items():
        if len(items) < 2:
            continue
        vals = [float(i["value_nm"]) for i in items]
        relations = {i.get("standard_relation") or "=" for i in items}
        if (max(vals) / min(vals)) >= divergence_fold or len(relations) > 1:
            fold = round(max(vals) / min(vals), 1)
            mol_id, tgt_id, org, stype = tkey
            msg = (
                f"Расхождение между измерениями {stype} для {mol_id} ({tgt_id}, {org}): "
                f"диапазон {min(vals)}–{max(vals)} nM (разброс {fold}x, отношения {sorted(relations)})."
            )
            if msg not in limitations:
                limitations.append(msg)
            for item in items:
                item["conflicting_measurements"] = True
                item["divergence_fold"] = fold

    return deduped, limitations


def resolve_molecule(query: str, *, use_cache: bool = True) -> dict[str, Any]:
    """Найти molecule_chembl_id по названию, ChEMBL ID или SMILES."""
    if not isinstance(query, str) or not query.strip():
        return {"query": "", "source": "chembl", "error": "empty_query", "is_cached": False}

    key = query.strip()
    lower = key.lower()

    if lower in _MOLECULE_ALIASES:
        info = dict(_MOLECULE_ALIASES[lower])
        info.update(
            {
                "query": key,
                "source": "chembl_alias",
                "is_cached": True,
                "fetched_at": "offline-alias",
                "ambiguous": False,
            }
        )
        return info

    if re.match(r"^chembl\d+$", lower):
        return {
            "query": key,
            "molecule_chembl_id": key.upper(),
            "pref_name": key.upper(),
            "source": "chembl_id",
            "is_cached": False,
            "fetched_at": _now_iso(),
            "ambiguous": False,
        }

    if use_cache:
        cached = _read_cache("mol", lower)
        if cached is not None:
            payload, fetched_at = cached
            res = dict(payload)
            res["is_cached"] = True
            res["fetched_at"] = fetched_at
            return res

    is_smiles = False
    if not re.match(r"^[A-Za-zА-Яа-яЁё\-]{4,}$", key):
        try:
            parse_smiles(key)
            is_smiles = True
        except ValueError:
            is_smiles = False

    if is_smiles:
        clean_smi = standardize_smiles(key, desalt=True)["canonical_smiles"]
        url = (
            f"{_CHEMBL_API}/molecule.json"
            f"?molecule_structures__canonical_smiles__flexmatch={urllib.parse.quote(clean_smi)}&limit=5"
        )
    else:
        url = f"{_CHEMBL_API}/molecule/search.json?q={urllib.parse.quote(key)}&limit=5"

    data = _get_json(url)
    if data and isinstance(data.get("molecules"), list) and data["molecules"]:
        mols = data["molecules"]
        first = mols[0]
        structs = first.get("molecule_structures") or {}
        fetched_at = _now_iso()
        result = {
            "query": key,
            "molecule_chembl_id": first.get("molecule_chembl_id"),
            "pref_name": first.get("pref_name"),
            "canonical_smiles": structs.get("canonical_smiles"),
            "ambiguous": len(mols) > 1 and not is_smiles,
            "candidates": [
                {
                    "molecule_chembl_id": m.get("molecule_chembl_id"),
                    "pref_name": m.get("pref_name"),
                }
                for m in mols[:5]
            ],
            "source": "chembl",
            "is_cached": False,
            "fetched_at": fetched_at,
        }
        _write_cache("mol", lower, result, fetched_at=fetched_at)
        return result

    return {
        "query": key,
        "molecule_chembl_id": None,
        "source": "chembl",
        "error": "no_data",
        "is_cached": False,
        "fetched_at": _now_iso(),
        "ambiguous": False,
    }


def resolve_target_id(target_id_or_name: str) -> str:
    """Привести название мишени или код к ChEMBL ID (например, CHEMBL205 для CA II)."""
    if not isinstance(target_id_or_name, str) or not target_id_or_name.strip():
        return "CHEMBL205"
    clean = target_id_or_name.strip()
    lower = clean.lower()
    if lower in _TARGET_ALIASES:
        return _TARGET_ALIASES[lower]
    if re.match(r"^chembl\d+$", lower):
        return clean.upper()
    return clean


def fetch_activities(
    target_id: str,
    molecule_id: str,
    *,
    limit: int = 50,
    use_cache: bool = True,
) -> list[dict[str, Any]]:
    """Получить и нормализовать данные биологической активности из ChEMBL для пары (мишень, молекула)."""
    norm_target = resolve_target_id(target_id)
    mol_info = resolve_molecule(molecule_id, use_cache=use_cache)
    norm_mol = mol_info.get("molecule_chembl_id")
    if not norm_mol:
        if re.match(r"^chembl\d+$", molecule_id.strip(), re.IGNORECASE):
            norm_mol = molecule_id.strip().upper()
        else:
            return []

    cache_key = f"{norm_target}_{norm_mol}_lim{limit}"
    if use_cache:
        cached = _read_cache("act", cache_key)
        if cached is not None:
            payload, fetched_at = cached
            out: list[dict[str, Any]] = []
            for row in payload:
                r = dict(row)
                r["is_cached"] = True
                r["fetched_at"] = fetched_at
                out.append(r)
            return out

    import os
    offline_raw = _OFFLINE_ACTIVITIES.get((norm_target, norm_mol))
    if offline_raw and os.getenv("DEARMO_OFFLINE", "0") == "1":
        normalized = [
            _normalize_record(
                raw,
                is_cached=True,
                fetched_at="offline-reference",
                source="chembl_offline",
            )
            for raw in offline_raw
        ]
        deduped, _ = deduplicate_activities(normalized)
        return deduped

    if not re.match(r"^CHEMBL\d+$", norm_target, re.IGNORECASE):
        return []

    params = urllib.parse.urlencode(
        {
            "target_chembl_id": norm_target,
            "molecule_chembl_id": norm_mol,
            "standard_type__in": "IC50,Ki,Kd,EC50",
            "limit": limit,
        }
    )
    url = f"{_CHEMBL_API}/activity.json?{params}"
    data = _get_json(url)
    if data and isinstance(data.get("activities"), list) and data["activities"]:
        fetched_at = _now_iso()
        normalized = [
            _normalize_record(raw, is_cached=False, fetched_at=fetched_at, source="chembl")
            for raw in data["activities"]
        ]
        deduped, _ = deduplicate_activities(normalized)
        _write_cache("act", cache_key, deduped, fetched_at=fetched_at)
        return deduped

    if offline_raw:
        normalized = [
            _normalize_record(
                raw,
                is_cached=True,
                fetched_at="offline-reference",
                source="chembl_offline",
            )
            for raw in offline_raw
        ]
        deduped, _ = deduplicate_activities(normalized)
        return deduped

    return []
