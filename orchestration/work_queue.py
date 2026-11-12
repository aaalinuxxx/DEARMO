"""Очередь частичных результатов воркеров карманов."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from models.result import Evidence

WorkStatus = Literal["pending", "running", "done", "failed"]


@dataclass
class WorkItem:
    pocket_id: str
    status: WorkStatus = "pending"
    evidence: list[Evidence] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    agent: str = ""
    rag_hits: int = 0
