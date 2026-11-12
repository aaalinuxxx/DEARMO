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


# Термин в тексте -> методы Evidence, которым это число имеет право принадлежать.
# Проверки «число вообще встречается среди фактов» недостаточно: пересказ может
# взять верное число и приписать его не тому измерению (IC50 вместо оценки Vina).
_TERM_TO_METHODS: tuple[tuple[re.Pattern[str], tuple[str, ...]], ...] = (
    (re.compile(r"vina[^.;\n]{0,40}?(-?\d+[.,]?\d*)", re.I), ("vina_score",)),
    (re.compile(r"(-?\d+[.,]?\d*)[^.;\n]{0,30}?ккал/моль", re.I), ("vina_score",)),
    (re.compile(r"qed[^.;\n]{0,40}?(-?\d+[.,]?\d*)", re.I), ("qed",)),
    (re.compile(r"logp[^.;\n]{0,40}?(-?\d+[.,]?\d*)", re.I), ("logp",)),
    (re.compile(r"молекулярн\w+ масс\w+[^.;\n]{0,30}?(-?\d+[.,]?\d*)", re.I), ("molecular_weight",)),
    (re.compile(r"разрешени\w+[^.;\n]{0,30}?(-?\d+[.,]?\d*)", re.I), ("resolution",)),
    (re.compile(r"\bIC50\b[^.;\n]{0,30}?(-?\d+[.,]?\d*)", re.I), ("IC50", "chembl_ic50")),
    (re.compile(r"\bKi\b[^.;\n]{0,30}?(-?\d+[.,]?\d*)", re.I), ("Ki", "chembl_ki")),
    (re.compile(r"\bKd\b[^.;\n]{0,30}?(-?\d+[.,]?\d*)", re.I), ("Kd",)),
    (re.compile(r"pubchem[^.;\n]{0,40}?(\d+)", re.I), ("pubchem_cid",)),
)


def _values_for_methods(evidence: list, methods: tuple[str, ...]) -> list[float]:
    wanted = {m.lower() for m in methods}
    values: list[float] = []
    for item in evidence:
        method = str(item.method or "").lower()
        # vina_score_site1 и best_vina_* — те же измерения под уточнёнными именами.
        if method in wanted or any(method.startswith(f"{m}_") or m in method for m in wanted):
            try:
                values.append(float(item.value))
            except (TypeError, ValueError):
                continue
    return values


def _check_number_method_binding(explanation: str, evidence: list) -> list[str]:
    """Число рядом с названием величины обязано принадлежать именно этой величине."""
    issues: list[str] = []
    for pattern, methods in _TERM_TO_METHODS:
        for match in pattern.finditer(explanation):
            try:
                stated = float(match.group(1).replace(",", "."))
            except (TypeError, ValueError):
                continue
            allowed = _values_for_methods(evidence, methods)
            if not allowed:
                continue
            if not any(abs(stated - ref) <= max(0.05, abs(ref) * 0.02) for ref in allowed):
                issues.append(
                    f"Число {match.group(1)} приписано величине «{methods[0]}», "
                    f"но среди фактов у неё значения {allowed}"
                )
    return issues


def review(state: WorkflowState) -> WorkflowState:
    with state.agent_step("review_agent", "validate"):
        issues: list[str] = []
        for ev in state.evidence:
            if ev.kind not in _ALLOWED_KINDS:
                issues.append(f"Недопустимый kind: {ev.kind}")
            if ev.kind == "calculation" and not ev.method and not ev.source:
                issues.append("calculation без method/source")
        # Детерминированный пересказ печатается прямо из Evidence, поэтому проверяется
        # только текст, сочинённый языковой моделью.
        if state.explanation and state.observations.get("explanation_source") != "deterministic":
            known_strs, known_floats = _evidence_numeric_values(state.evidence)
            for match in _NUMBER_IN_TEXT.finditer(state.explanation):
                token = match.group(0).replace(",", ".")
                if not _is_supported_number(token, known_strs, known_floats):
                    if len(token) > 1:
                        issues.append(f"Число в explanation без Evidence: {token}")
            issues.extend(_check_number_method_binding(state.explanation, state.evidence))
        if state.scenario_id and not state.evidence and not state.errors:
            issues.append("Нет evidence при успешном сценарии")
        state.review_errors = issues
        state.review_passed = len(issues) == 0
        if issues:
            state.add_error("Review: " + "; ".join(issues))
        state.trace("review_agent", "done", detail=f"passed={state.review_passed}")
        return state
