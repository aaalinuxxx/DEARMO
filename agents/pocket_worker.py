"""Один карман — один воркер. Не пишет в общий state."""
from __future__ import annotations

import os
from pathlib import Path

from models.docking import PoseRow
from models.result import Evidence
from orchestration.work_queue import WorkItem
from tools.literature_rag import retrieve


def _min_score() -> float:
    return float(os.environ.get("ORC_RAG_MIN_SCORE", "0.5"))


def run_pocket(
    *,
    pocket_id: str,
    rows: list[PoseRow],
    ligand_id: str | None,
    query_terms: list[str],
) -> WorkItem:
    agent = f"pocket_worker_{pocket_id}"
    item = WorkItem(pocket_id=pocket_id, status="running", agent=agent)
    pocket_rows = [row for row in rows if row.pocket_id == pocket_id]
    if ligand_id:
        preferred = [row for row in pocket_rows if row.ligand_id == ligand_id]
        if preferred:
            pocket_rows = preferred
    if not pocket_rows:
        item.status = "failed"
        item.errors.append(f"Нет поз для {pocket_id}")
        return item
    best = min(pocket_rows, key=lambda row: row.vina_score)
    item.evidence.append(
        Evidence(
            kind="calculation",
            value=best.vina_score,
            unit="kcal/mol",
            method=f"vina_score_{pocket_id}",
            source="AutoDock Vina",
            context={
                "ligand_id": best.ligand_id,
                "pocket_id": pocket_id,
                "pose_file": best.pose_file,
                "note": "Оценка Vina — расчёт, не экспериментальная аффинность.",
            },
        )
    )
    if best.p2rank_score is not None:
        item.evidence.append(
            Evidence(
                kind="calculation",
                value=best.p2rank_score,
                method=f"p2rank_score_{pocket_id}",
                source="P2Rank",
                context={"pocket_id": pocket_id},
            )
        )
    chunks = retrieve(query_terms, min_score=_min_score(), limit=2)
    item.rag_hits = len(chunks)
    for chunk in chunks:
        item.evidence.append(
            Evidence(
                kind="database",
                value=chunk.get("title") or chunk.get("chunk_id"),
                method=f"rag_{pocket_id}",
                source=chunk.get("source", "rag"),
                context={"score": chunk.get("score"), "pocket_id": pocket_id},
            )
        )
    item.status = "done"
    return item
