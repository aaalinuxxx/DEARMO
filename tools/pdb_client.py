"""RCSB PDB: метаданные структур, состав комплекса и поиск по мишени.

Поиск структуры устроен лестницей из четырёх уровней (см. `search_structures`).
Уровень, на котором нашлись кандидаты, сам по себе является результатом: он
говорит, насколько структура соответствует задаче, и попадает в ответ как факт,
а не как оговорка в тексте.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools.http_json import Transport, request_json

_CACHE_ROOT = Path(__file__).resolve().parent.parent / "data" / "cache" / "pdb"
_EXTERNAL = Path(__file__).resolve().parent.parent / "data" / "external" / "rcsb"
_RCSB_REST = "https://data.rcsb.org/rest/v1/core/entry"
_RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_CORE = "https://data.rcsb.org/rest/v1/core"

_TRANSPORT = Transport()

# Кофакторы и ионы, без которых карман не тот: их наличие проверяется отдельно.
# Список намеренно не содержит привязки к конкретному белку — это просто
# разметка «ион металла», и решение о значимости принимает агент по контексту.
_METAL_IONS = frozenset(
    {"ZN", "MG", "MN", "FE", "FE2", "CU", "CU1", "CA", "NI", "CO", "CD", "K", "NA"}
)
# Частые артефакты кристаллизации: присутствуют почти везде и ничего не значат.
_CRYSTALLIZATION_ARTIFACTS = frozenset(
    {"HOH", "GOL", "EDO", "PEG", "PG4", "SO4", "PO4", "CL", "ACT", "DMS", "MPD", "TRS", "IMD"}
)


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _file_mtime_iso(path: Path) -> str:
    try:
        return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except OSError:
        return _now_iso()


def transport_status() -> dict[str, Any]:
    return _TRANSPORT.as_dict()


def _get_json_with_retry(url: str, *, timeout: float = 10.0, max_retries: int = 3) -> dict[str, Any] | None:
    return request_json(url, timeout=timeout, max_retries=max_retries, transport=_TRANSPORT)


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


def _search(query: dict[str, Any], *, limit: int, timeout: float, sort_by_resolution: bool) -> list[str]:
    """Вернуть идентификаторы структур по готовому запросу Search API."""
    options: dict[str, Any] = {"paginate": {"start": 0, "rows": limit}}
    if sort_by_resolution:
        options["sort"] = [
            {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
        ]
    payload = request_json(
        _RCSB_SEARCH,
        body={"query": query, "return_type": "entry", "request_options": options},
        timeout=timeout,
        max_retries=2,
        transport=_TRANSPORT,
    )
    if not payload:
        return []
    found: list[str] = []
    for hit in payload.get("result_set") or []:
        pdb_id = str(hit.get("identifier") or "").upper()[:4]
        if len(pdb_id) == 4 and pdb_id not in found:
            found.append(pdb_id)
    return found


def _attr(attribute: str, operator: str, value: Any) -> dict[str, Any]:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def _uniprot_node(accession: str) -> dict[str, Any]:
    return _attr(
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "exact_match",
        accession.strip().upper(),
    )


def chem_comp_inchikey(comp_id: str, *, timeout: float = 10.0) -> str | None:
    data = _get_json_with_retry(f"{_RCSB_CORE}/chemcomp/{comp_id.strip().upper()}", timeout=timeout)
    if not data:
        return None
    return (data.get("rcsb_chem_comp_descriptor") or {}).get("InChIKey")


def _same_skeleton(left: str | None, right: str | None) -> bool:
    """Сравнение по связности: первый блок InChIKey одинаков у стереоизомеров.

    Кристалл может содержать другой стереоизомер или соль того же соединения —
    это та же молекула для целей выбора структуры. А вот другой скелет это уже
    другое вещество, и подставлять его нельзя.
    """
    if not left or not right:
        return False
    return left.split("-")[0] == right.split("-")[0]


def find_chem_comp_ids(
    smiles: str,
    *,
    inchikey: str | None = None,
    limit: int = 5,
    timeout: float = 15.0,
) -> list[str]:
    """Трёхбуквенные коды компонентов PDB, соответствующие структуре лиганда.

    Нужны, чтобы спросить у RCSB «структуры этого белка именно с этим лигандом».
    Сопоставление графовое, по самой молекуле, а не по названию; `graph-relaxed`
    возвращает и похожие соединения, поэтому при известном InChIKey выдача
    дополнительно сверяется поатомно и остаются только совпадения по скелету.
    """
    clean = str(smiles or "").strip()
    if not clean:
        return []
    query = {
        "type": "terminal",
        "service": "chemical",
        "parameters": {
            "value": clean,
            "type": "descriptor",
            "descriptor_type": "SMILES",
            "match_type": "graph-relaxed",
        },
    }
    payload = request_json(
        _RCSB_SEARCH,
        body={
            "query": query,
            "return_type": "mol_definition",
            "request_options": {"paginate": {"start": 0, "rows": limit}},
        },
        timeout=timeout,
        max_retries=2,
        transport=_TRANSPORT,
    )
    if not payload:
        return []
    found = [
        str(hit.get("identifier") or "").upper()
        for hit in payload.get("result_set") or []
        if hit.get("identifier")
    ]
    if not inchikey:
        # Без эталонного InChIKey проверить нечем: берём только первое, самое
        # близкое совпадение, чтобы не тянуть в запрос похожие соединения.
        return found[:1]
    return [comp for comp in found if _same_skeleton(chem_comp_inchikey(comp), inchikey)]


# Уровни поиска: от самого строгого соответствия задаче к самому слабому.
# Текст — то, что увидит пользователь рядом с выбранной структурой.
SEARCH_TIERS: dict[str, str] = {
    "A": "комплекс этого белка именно с этим лигандом",
    "B": "структура этого белка с каким-то связанным лигандом (карман в holo-форме)",
    "C": "структура этого белка без связанного лиганда (apo; геометрия кармана может отличаться)",
    "D": "полнотекстовый поиск по названию: мишень не удалось сопоставить с UniProt",
}


def search_structures(
    *,
    uniprot_accession: str | None = None,
    ligand_comp_ids: list[str] | None = None,
    text_query: str | None = None,
    limit: int = 5,
    timeout: float = 15.0,
) -> dict[str, Any]:
    """Спуск по уровням A→D до первого непустого результата.

    Возвращает `{"tier", "tier_note", "pdb_ids", "query_description"}`. Уровень
    выбирается данными, а не настройкой: если у белка есть кристалл с этим
    лигандом, мы его найдём, а если нет — честно скажем, что берём apo-форму.
    """
    accession = (uniprot_accession or "").strip().upper()
    comp_ids = [c for c in (ligand_comp_ids or []) if c]

    if accession:
        if comp_ids:
            # Уровень A: любой из подходящих компонентов связан с этим белком.
            ligand_node: dict[str, Any] = (
                _attr(
                    "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                    "exact_match",
                    comp_ids[0],
                )
                if len(comp_ids) == 1
                else {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        _attr(
                            "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                            "exact_match",
                            comp_id,
                        )
                        for comp_id in comp_ids
                    ],
                }
            )
            found = _search(
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [_uniprot_node(accession), ligand_node],
                },
                limit=limit,
                timeout=timeout,
                sort_by_resolution=True,
            )
            if found:
                return {
                    "tier": "A",
                    "tier_note": SEARCH_TIERS["A"],
                    "pdb_ids": found,
                    "query_description": f"UniProt {accession} + компонент {'/'.join(comp_ids)}",
                }

        # Уровень B: holo-структура с любым лигандом, рентген, по разрешению.
        found = _search(
            {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    _uniprot_node(accession),
                    _attr("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
                    _attr("rcsb_entry_info.nonpolymer_entity_count", "greater", 0),
                ],
            },
            limit=limit,
            timeout=timeout,
            sort_by_resolution=True,
        )
        if found:
            return {
                "tier": "B",
                "tier_note": SEARCH_TIERS["B"],
                "pdb_ids": found,
                "query_description": f"UniProt {accession} + рентген + связанный лиганд",
            }

        # Уровень C: хоть какая-то структура этого белка.
        found = _search(
            _uniprot_node(accession), limit=limit, timeout=timeout, sort_by_resolution=True
        )
        if found:
            return {
                "tier": "C",
                "tier_note": SEARCH_TIERS["C"],
                "pdb_ids": found,
                "query_description": f"UniProt {accession}",
            }

    # Уровень D: мишень не разрешилась — остаётся искать словами.
    text = " ".join(str(text_query or "").split())
    if text:
        found = _search(
            {"type": "terminal", "service": "full_text", "parameters": {"value": text}},
            limit=limit,
            timeout=timeout,
            sort_by_resolution=False,
        )
        if found:
            return {
                "tier": "D",
                "tier_note": SEARCH_TIERS["D"],
                "pdb_ids": found,
                "query_description": f"полнотекстовый поиск «{text}»",
            }

    return {"tier": None, "tier_note": None, "pdb_ids": [], "query_description": None}


def count_complexes(
    uniprot_accession: str,
    ligand_comp_ids: list[str],
    *,
    timeout: float = 12.0,
) -> int | None:
    """Сколько структур PDB содержат этот белок в комплексе с этим лигандом.

    Это независимое от BindingDB экспериментальное свидетельство: если по
    лиганду есть десятки кристаллов с одним белком и ни одного с другим, то
    первый — настоящая мишень, даже когда в выдаче BindingDB он представлен
    лишь мутантом. None означает «не удалось узнать», а не «ноль».
    """
    accession = str(uniprot_accession or "").strip().upper()
    comp_ids = [c for c in ligand_comp_ids or [] if c]
    if not accession or not comp_ids:
        return None
    ligand_node: dict[str, Any] = (
        _attr("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "exact_match", comp_ids[0])
        if len(comp_ids) == 1
        else {
            "type": "group",
            "logical_operator": "or",
            "nodes": [
                _attr("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "exact_match", c)
                for c in comp_ids
            ],
        }
    )
    payload = request_json(
        _RCSB_SEARCH,
        body={
            "query": {"type": "group", "logical_operator": "and", "nodes": [_uniprot_node(accession), ligand_node]},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1}},
        },
        timeout=timeout,
        max_retries=2,
        transport=_TRANSPORT,
    )
    if not payload:
        return None
    total = payload.get("total_count")
    return int(total) if isinstance(total, int) else 0


def get_entity_details(pdb_id: str, *, timeout: float = 12.0) -> dict[str, Any]:
    """Цепи, их UniProt и состав гетерогрупп структуры.

    Кейс требует назвать цепь и значимые кофакторы (например, Zn²⁺ для
    карбоангидразы). Это не выводится из метаданных записи — нужны отдельные
    запросы к полимерным и неполимерным сущностям.
    """
    clean_id = str(pdb_id).strip().upper()[:4]
    entry = _get_json_with_retry(f"{_RCSB_REST}/{clean_id}", timeout=timeout)
    if not entry:
        return {"pdb_id": clean_id, "chains": [], "heteroatoms": [], "available": False}

    ids = entry.get("rcsb_entry_container_identifiers") or {}
    chains: list[dict[str, Any]] = []
    for entity_id in ids.get("polymer_entity_ids") or []:
        data = _get_json_with_retry(
            f"{_RCSB_CORE}/polymer_entity/{clean_id}/{entity_id}", timeout=timeout
        )
        if not data:
            continue
        container = data.get("rcsb_polymer_entity_container_identifiers") or {}
        accessions = [
            ref.get("database_accession")
            for ref in container.get("reference_sequence_identifiers") or []
            if (ref.get("database_name") or "").lower().startswith("uniprot")
        ]
        chains.append(
            {
                "entity_id": entity_id,
                "auth_chains": container.get("auth_asym_ids") or [],
                "name": ((data.get("rcsb_polymer_entity") or {}).get("pdbx_description")),
                "uniprot_accessions": [acc for acc in accessions if acc],
            }
        )

    heteroatoms: list[dict[str, Any]] = []
    for entity_id in ids.get("non_polymer_entity_ids") or []:
        data = _get_json_with_retry(
            f"{_RCSB_CORE}/nonpolymer_entity/{clean_id}/{entity_id}", timeout=timeout
        )
        if not data:
            continue
        nonpoly = data.get("pdbx_entity_nonpoly") or {}
        comp_id = str(nonpoly.get("comp_id") or "").upper()
        container = data.get("rcsb_nonpolymer_entity_container_identifiers") or {}
        heteroatoms.append(
            {
                "comp_id": comp_id,
                "name": nonpoly.get("name"),
                "auth_chains": container.get("auth_asym_ids") or [],
                "role": (
                    "metal_ion"
                    if comp_id in _METAL_IONS
                    else "crystallization_artifact"
                    if comp_id in _CRYSTALLIZATION_ARTIFACTS
                    else "bound_ligand"
                ),
            }
        )

    return {
        "pdb_id": clean_id,
        "chains": chains,
        "heteroatoms": heteroatoms,
        "available": True,
        "fetched_at": _now_iso(),
    }


def significant_cofactors(details: dict[str, Any]) -> list[dict[str, Any]]:
    """Гетерогруппы, которые имеют отношение к функции, без буферов и глицерина."""
    return [
        het
        for het in details.get("heteroatoms") or []
        if het.get("role") in ("metal_ion", "bound_ligand")
    ]
