"""Простой RAG: заранее размеченные чанки в data/rag/."""
from __future__ import annotations

import json
from pathlib import Path

_RAG_ROOT = Path(__file__).resolve().parent.parent / "data" / "rag"


def retrieve(query_terms: list[str], *, limit: int = 3, min_score: float = 0.0) -> list[dict]:
    if not _RAG_ROOT.is_dir():
        return []
    terms = [t.lower() for t in query_terms if t]
    scored: list[tuple[float, dict]] = []
    for path in _RAG_ROOT.glob("*.json"):
        with path.open(encoding="utf-8") as f:
            chunk = json.load(f)
        text = (chunk.get("text") or "").lower()
        tags = " ".join(chunk.get("tags") or []).lower()
        hits = sum(1 for t in terms if t in text or t in tags)
        if not terms or hits == 0:
            continue
        score = hits / len(terms)
        if score < min_score:
            continue
        scored.append((score, {**chunk, "chunk_id": path.stem, "score": round(score, 3)}))
    scored.sort(key=lambda x: -x[0])
    return [c for _, c in scored[:limit]]
