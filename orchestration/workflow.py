"""Маршрутизация сценариев, лимиты шагов, сборка ResearchResult."""
from __future__ import annotations

from agents import planner, review_agent
from models.request import ResearchRequest
from models.result import ResearchResult
from orchestration.executor import execute_plan
from orchestration.gigachat_client import narrate_with_gigachat
from orchestration.state import WorkflowState


def _build_explanation(state: WorkflowState) -> str:
    lines = []
    for ev in state.evidence[:12]:
        src = ev.source or ev.method or ev.kind
        lines.append(f"- [{ev.kind}] {ev.method or 'value'}={ev.value} ({src})")
    summary = "\n".join(lines) if lines else "Данных для объяснения нет."
    narrated = narrate_with_gigachat(summary, state.limitations)
    if narrated:
        return narrated
    if state.scenario_id:
        return f"Сценарий {state.scenario_id} выполнен. См. таблицу доказательств."
    return "Не удалось классифицировать запрос."


def _result_status(state: WorkflowState) -> str:
    if state.errors and not state.evidence:
        return "error"
    if any("уточн" in e.lower() for e in state.errors):
        return "needs_clarification"
    if not state.evidence:
        return "no_data"
    if not state.review_passed:
        return "error"
    return "success"


def run_workflow(request: ResearchRequest) -> tuple[ResearchResult, WorkflowState]:
    state = WorkflowState(request=request)
    try:
        state = planner.plan(state)
    except RuntimeError as exc:
        state.add_error(str(exc))
        return _finalize(state)

    if state.scenario_id is None:
        state.add_error("Не удалось определить сценарий. Уточните: мишень, докинг, сравнение или валидация позы.")
        return _finalize(state)

    if state.scenario_id == "target_for_ligand" and "ligand_query" not in state.observations:
        state.add_error("Укажите название лиганда или SMILES.")
        return _finalize(state)

    try:
        state = execute_plan(state)
    except Exception as exc:  # noqa: BLE001 — граница workflow
        state.add_error(f"Ошибка сценария: {exc}")

    state = review_agent.review(state)
    state.explanation = _build_explanation(state)
    return _finalize(state)


def _finalize(state: WorkflowState) -> tuple[ResearchResult, WorkflowState]:
    status = _result_status(state)
    if status == "needs_clarification" and not state.errors:
        state.add_error("Недостаточно входных данных для сценария.")
    result = ResearchResult(
        status=status,
        evidence=list(state.evidence),
        explanation=state.explanation,
        limitations=list(state.limitations),
        errors=list(state.errors),
    )
    return result, state
