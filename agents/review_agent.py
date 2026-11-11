"""Агент проверки: программная валидация Evidence и формулировок."""
from __future__ import annotations

import re

from models.result import EvidenceKind
from orchestration.state import WorkflowState

_ALLOWED_KINDS: set[EvidenceKind] = {"database", "calculation", "prediction", "hypothesis"}
_NUMBER_IN_TEXT = re.compile(r"-?\d+(?:[.,]\d+)?")


def _evidence_numbers(evidence: list) -> set[str]:
    nums: set[str] = set()
    for item in evidence:
        val = item.value
        if isinstance(val, (int, float)):
            nums.add(str(val))
            nums.add(f"{val:.2f}".rstrip("0").rstrip("."))
    return nums


def review(state: WorkflowState) -> WorkflowState:
    with state.agent_step("review_agent", "validate"):
        issues: list[str] = []
        for ev in state.evidence:
            if ev.kind not in _ALLOWED_KINDS:
                issues.append(f"Недопустимый kind: {ev.kind}")
            if ev.kind == "calculation" and not ev.method and not ev.source:
                issues.append("calculation без method/source")
        if state.explanation:
            known = _evidence_numbers(state.evidence)
            for match in _NUMBER_IN_TEXT.finditer(state.explanation):
                token = match.group(0).replace(",", ".")
                if token not in known and not any(token in k for k in known):
                    # Допускаем целые без точного совпадения только для малых шумов — иначе флаг
                    if len(token) > 1:
                        issues.append(f"Число в explanation без Evidence: {token}")
        if state.scenario_id and not state.evidence and not state.errors:
            issues.append("Нет evidence при успешном сценарии")
        state.review_errors = issues
        state.review_passed = len(issues) == 0
        if issues:
            state.add_error("Review: " + "; ".join(issues))
        state.trace("review_agent", "done", detail=f"passed={state.review_passed}")
    return state
