"""Агент проверки: программная валидация Evidence и численных утверждений с допуском на округление."""
from __future__ import annotations

import re

from models.result import EvidenceKind
from orchestration.state import WorkflowState

_ALLOWED_KINDS: set[EvidenceKind] = {"database", "calculation", "prediction", "hypothesis"}
_NUMBER_IN_TEXT = re.compile(r"-?\d+(?:[.,]\d+)?")


def _evidence_numeric_values(evidence: list) -> tuple[set[str], list[float]]:
    str_nums: set[str] = set()
    float_nums: list[float] = []
    for item in evidence:
        val = item.value
        if isinstance(val, (int, float)) and not isinstance(val, bool):
            fval = float(val)
            float_nums.append(fval)
            str_nums.add(str(val))
            str_nums.add(f"{fval:.2f}".rstrip("0").rstrip("."))
            str_nums.add(f"{fval:.1f}".rstrip("0").rstrip("."))
            str_nums.add(str(round(fval)))
    return str_nums, float_nums


def _is_supported_number(token: str, known_strs: set[str], known_floats: list[float]) -> bool:
    if token in known_strs or any(token in k for k in known_strs):
        return True
    try:
        val = float(token)
    except ValueError:
        return False
    # Числовой допуск на округление LLM (до 2% или 0.05 абсолютных единиц)
    for ref in known_floats:
        tol = max(0.05, abs(ref) * 0.02)
        if abs(val - ref) <= tol:
            return True
    return False


def review(state: WorkflowState) -> WorkflowState:
    with state.agent_step("review_agent", "validate"):
        issues: list[str] = []
        for ev in state.evidence:
            if ev.kind not in _ALLOWED_KINDS:
                issues.append(f"Недопустимый kind: {ev.kind}")
            if ev.kind == "calculation" and not ev.method and not ev.source:
                issues.append("calculation без method/source")
        if state.explanation:
            known_strs, known_floats = _evidence_numeric_values(state.evidence)
            for match in _NUMBER_IN_TEXT.finditer(state.explanation):
                token = match.group(0).replace(",", ".")
                if not _is_supported_number(token, known_strs, known_floats):
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
