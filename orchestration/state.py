"""Независимое состояние обработки одного исследовательского запроса."""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

from models.request import ResearchRequest
from models.result import Evidence

MAX_WORKFLOW_STEPS = 40
MAX_REVIEW_RETRIES = 1


@dataclass
class WorkflowState:
    request: ResearchRequest
    scenario_id: str | None = None
    plan: list[str] = field(default_factory=list)
    evidence: list[Evidence] = field(default_factory=list)
    observations: dict[str, Any] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    agent_trace: list[dict[str, Any]] = field(default_factory=list)
    artifacts: dict[str, str] = field(default_factory=dict)
    rag_context: list[dict[str, Any]] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    review_passed: bool = False
    review_errors: list[str] = field(default_factory=list)
    step_count: int = 0
    explanation: str | None = None
    planner_source: str | None = None

    def add_evidence(self, evidence: Evidence) -> None:
        self.evidence.append(evidence)

    def add_error(self, message: str) -> None:
        if not message.strip():
            raise ValueError("Сообщение об ошибке не должно быть пустым")
        self.errors.append(message)

    def trace(self, agent: str, action: str, *, detail: str | None = None, duration_s: float | None = None) -> None:
        entry: dict[str, Any] = {
            "agent": agent,
            "action": action,
            "step": self.step_count,
        }
        if detail:
            entry["detail"] = detail
        if duration_s is not None:
            entry["duration_s"] = round(duration_s, 3)
        self.agent_trace.append(entry)

    def bump_step(self) -> None:
        self.step_count += 1
        if self.step_count > MAX_WORKFLOW_STEPS:
            raise RuntimeError("Превышен лимит шагов workflow")

    class _AgentStep:
        def __init__(self, state: "WorkflowState", agent: str, action: str):
            self._state = state
            self._agent = agent
            self._action = action
            self._t0 = time.perf_counter()

        def __enter__(self) -> "WorkflowState":
            self._state.bump_step()
            return self._state

        def __exit__(self, exc_type, exc, tb) -> None:
            detail = str(exc) if exc else None
            self._state.trace(
                self._agent,
                self._action,
                detail=detail,
                duration_s=time.perf_counter() - self._t0,
            )

    def agent_step(self, agent: str, action: str) -> _AgentStep:
        return WorkflowState._AgentStep(self, agent, action)
