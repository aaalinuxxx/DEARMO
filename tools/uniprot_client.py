"""UniProt: имя мишени из BindingDB → проверенный accession.

Почему не спросить у модели. Qwen знает, что у карбоангидразы II accession
P00918, и ответит правильно — это и опасно: на незнакомом белке она ответит так
же уверенно и так же правдоподобно, а проверить вывод будет нечем. Поэтому
accession приходит только из REST UniProt, и только если имя из записи
действительно совпало с запрошенным. Не совпало — возвращаем None, и агент
честно пишет мишень без идентификатора.
"""
from __future__ import annotations

import re
import urllib.parse
from typing import Any

from tools.http_json import Transport, request_json

_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_FIELDS = "accession,id,protein_name,gene_primary,organism_name,cc_function"

_TRANSPORT = Transport()

# Слова, которые UniProt и BindingDB добавляют к одному и тому же белку
# по-разному; они не несут различающего смысла и мешают сравнению имён.
_NOISE = re.compile(
    r"\b(domain[- ]containing|protein|isoform|precursor|putative|probable|fragment)\b"
)


def transport_status() -> dict[str, Any]:
    return _TRANSPORT.as_dict()


def reset_transport_status() -> None:
    _TRANSPORT.reset()


def _normalize(text: str) -> str:
    lowered = str(text or "").lower()
    lowered = _NOISE.sub(" ", lowered)
    return " ".join(re.sub(r"[^a-z0-9]+", " ", lowered).split())


def _names_of(entry: dict[str, Any]) -> list[str]:
    """Рекомендованное и альтернативные имена записи — всё, с чем сравниваем."""
    desc = entry.get("proteinDescription") or {}
    names: list[str] = []
    recommended = (desc.get("recommendedName") or {}).get("fullName", {}).get("value")
    if recommended:
        names.append(recommended)
    for alt in desc.get("alternativeNames") or []:
        value = (alt.get("fullName") or {}).get("value")
        if value:
            names.append(value)
    for sub in desc.get("submissionNames") or []:
        value = (sub.get("fullName") or {}).get("value")
        if value:
            names.append(value)
    return names


def _match_quality(query: str, entry: dict[str, Any]) -> str | None:
    """`exact`, `partial` или None. None означает «не подтвердили»."""
    wanted = _normalize(query)
    if not wanted:
        return None
    variants = [_normalize(name) for name in _names_of(entry)]
    gene = ((entry.get("genes") or [{}])[0].get("geneName") or {}).get("value")
    if gene:
        variants.append(_normalize(gene))
    if wanted in variants:
        return "exact"
    for variant in variants:
        if not variant:
            continue
        if wanted in variant or variant in wanted:
            return "partial"
    return None


def _entry_to_record(entry: dict[str, Any], quality: str) -> dict[str, Any]:
    desc = entry.get("proteinDescription") or {}
    genes = entry.get("genes") or []
    # У записей TrEMBL рекомендованного имени нет, только submissionName.
    all_names = _names_of(entry)
    recommended = (desc.get("recommendedName") or {}).get("fullName", {}).get("value")
    return {
        "accession": entry.get("primaryAccession"),
        "entry_name": entry.get("uniProtkbId"),
        "protein_name": recommended or (all_names[0] if all_names else None),
        "gene": ((genes[0].get("geneName") or {}).get("value") if genes else None),
        "organism": (entry.get("organism") or {}).get("scientificName"),
        "reviewed": entry.get("entryType", "").lower().startswith("uniprotkb reviewed"),
        "name_match": quality,
        "url": f"https://www.uniprot.org/uniprotkb/{entry.get('primaryAccession')}",
        "source": "UniProt",
    }


def resolve_target(
    target_name: str,
    organism: str | None = None,
    *,
    timeout: float = 12.0,
) -> dict[str, Any] | None:
    """Найти accession по названию белка и организму.

    Сначала ищем среди проверенных вручную записей (Swiss-Prot); если там ничего
    нет, повторяем без этого ограничения — для нечеловеческих мишеней запись
    часто только в TrEMBL. Имя записи обязано совпасть с запрошенным, иначе
    поисковая выдача UniProt подставит соседний белок, а мы этого не заметим.
    """
    name = str(target_name or "").strip()
    if not name:
        return None

    records = resolve_target_parts(name, organism, timeout=timeout)
    return records[0] if records else None


def resolve_target_parts(
    target_name: str,
    organism: str | None = None,
    *,
    timeout: float = 12.0,
) -> list[dict[str, Any]]:
    """Все записи UniProt, которым соответствует название мишени.

    Составные мишени BindingDB называет слитно: «Breakpoint cluster region
    protein/Tyrosine-protein kinase ABL1» — это слитый белок BCR-ABL, и такой
    записи в UniProt нет, а у половин есть свои. Какая из половин является
    мишенью лекарства, по названию не видно: для иматиниба это киназа ABL1, но
    первой в строке стоит BCR. Поэтому выбор здесь не делается — возвращаются
    обе записи, а разрешает неоднозначность поиск структур: комплекс с этим
    лигандом есть только у настоящей мишени.
    """
    name = str(target_name or "").strip()
    if not name:
        return []

    whole = _resolve_one(name, organism, timeout=timeout)
    if whole is not None:
        return [whole]

    found: list[dict[str, Any]] = []
    for part in re.split(r"\s*/\s*", name):
        part = part.strip()
        if not part or part == name:
            continue
        record = _resolve_one(part, organism, timeout=timeout)
        if record is not None and all(r["accession"] != record["accession"] for r in found):
            record["resolved_from_part"] = part
            record["full_target_name"] = name
            found.append(record)
    return found


def _resolve_one(
    name: str, organism: str | None, *, timeout: float
) -> dict[str, Any] | None:
    parts = [f'"{name}"']
    if organism and organism.strip():
        parts.append(f'organism_name:"{organism.strip()}"')

    for reviewed_only in (True, False):
        query = " AND ".join(parts + (["reviewed:true"] if reviewed_only else []))
        url = (
            f"{_SEARCH}?query={urllib.parse.quote(query)}"
            f"&fields={_FIELDS}&format=json&size=10"
        )
        data = request_json(url, timeout=timeout, transport=_TRANSPORT)
        if not data:
            continue
        results = data.get("results") or []
        # Точное совпадение имени важнее позиции в выдаче.
        for wanted_quality in ("exact", "partial"):
            for entry in results:
                if _match_quality(name, entry) == wanted_quality:
                    record = _entry_to_record(entry, wanted_quality)
                    record["query"] = query
                    return record
    return None
