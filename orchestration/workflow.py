"""Маршрутизация сценариев, лимиты шагов, ранний выход на уточнение и сборка ResearchResult."""
from __future__ import annotations

from agents import planner, review_agent
from models.request import ResearchRequest
from models.result import ResearchResult
from orchestration.executor import HANDLERS, execute_plan
from orchestration.llm_client import critique, narrate_with_llm
from orchestration.scenarios import missing_required_keys
from orchestration.state import WorkflowState


def _deterministic_explanation(state: WorkflowState) -> str:
    """Текст без LLM: перечисление собранного. Используется, если пересказ не прошёл проверку."""
    if not state.evidence:
        return "Не удалось классифицировать запрос." if not state.scenario_id else (
            f"Сценарий {state.scenario_id} выполнен, но измерений не собрано."
        )
    lines = [f"Выполнено: {', '.join(state.observations.get('completed_steps') or ['—'])}."]
    for ev in state.evidence[:12]:
        lines.append(f"- [{ev.kind}] {ev.method or 'значение'} = {ev.value} ({ev.source or '—'})")
    return "\n".join(lines)


def _build_explanation(state: WorkflowState) -> str:
    lines = []
    for ev in state.evidence[:12]:
        src = ev.source or ev.method or ev.kind
        lines.append(f"- [{ev.kind}] {ev.method or 'value'}={ev.value} ({src})")
    summary = "\n".join(lines) if lines else "Данных для объяснения нет."
    narrated = narrate_with_llm(summary, state.limitations)
    if narrated:
        return narrated
    return _deterministic_explanation(state)


def _evidence_digest(state: WorkflowState) -> str:
    """Сводка для критика: по одному представителю каждого измерения.

    Простая обрезка первых N фактов прятала от критика поздние шаги (оценку Vina,
    проверку позы), и он объявлял неполным ответ, который на деле был полным.
    """
    seen: set[str] = set()
    lines: list[str] = []
    for ev in state.evidence:
        method = str(ev.method or "value")
        if method in seen:
            continue
        seen.add(method)
        lines.append(f"- [{ev.kind}] {method}={ev.value} ({ev.source or '—'})")
    done = ", ".join(state.observations.get("completed_steps") or []) or "—"
    header = f"Выполненные шаги: {done}."
    return header + "\n" + ("\n".join(lines) if lines else "Фактов нет.")


def _critic_replan(state: WorkflowState) -> None:
    """LLM-критик смотрит на собранное и может потребовать ещё один шаг.

    Это обратная связь, а не конвейер: критик не придумывает числа, он называет
    недостающий сценарий. Допускается ровно один доп. шаг, чтобы не зациклиться.
    """
    if not state.evidence or state.observations.get("critic_ran"):
        return
    verdict = critique(state.request.text, _evidence_digest(state))
    if verdict is None:
        return
    state.observations["critic_ran"] = True
    state.observations["critic"] = verdict
    state.trace(
        "critic",
        "sufficient" if verdict["sufficient"] else "needs_more",
        detail=f"{verdict['source']}: {(verdict.get('missing') or 'всё есть')[:120]}",
    )
    if verdict["sufficient"]:
        return

    missing = verdict.get("missing") or "не названо"
    next_scenario = verdict.get("next_scenario")
    done = set(state.observations.get("completed_steps") or [])

    # Критик может назвать шаг, входных данных для которого нет (например, проверку
    # позы без докинга). Такой шаг не запускается: иначе в успешный ответ попадёт
    # ошибка о недостающем файле, которую пользователь не просил.
    blocked = missing_required_keys(state.observations, next_scenario) if next_scenario in HANDLERS else []
    if blocked:
        state.trace("critic", f"skip:{next_scenario}", detail=f"нет входных данных: {', '.join(blocked)}")
        state.limitations.append(
            f"Критик предложил шаг {next_scenario}, но для него нет данных ({', '.join(blocked)}). Шаг не выполнялся."
        )
        return

    if next_scenario in HANDLERS and next_scenario not in done:
        state.observations["scenario_ids"] = list(
            state.observations.get("scenario_ids") or []
        ) + [next_scenario]
        state.trace("critic", f"replan:{next_scenario}", detail=missing[:120])
        try:
            with state.agent_step("orchestrator", f"replan:{next_scenario}"):
                HANDLERS[next_scenario](state)
            done.add(next_scenario)
            state.observations["completed_steps"] = sorted(done)
        except Exception as exc:  # noqa: BLE001 — доп. шаг не должен ронять ответ
            state.limitations.append(f"Доп. шаг {next_scenario} по требованию критика не выполнен: {exc}")
    else:
        state.limitations.append(
            f"Критик отметил неполноту ответа: {missing}. Автоматически это не достраивается."
        )


def _needs_clarification(state: WorkflowState) -> bool:
    return any("уточн" in e.lower() for e in state.errors) or bool(
        state.observations.get("clarification_candidates")
    )


def _result_status(state: WorkflowState) -> str:
    # Проверка на неоднозначность/уточнение ДО проверки отсутствия evidence,
    # иначе неоднозначный запрос на этапе planner ошибочно получает status="error".
    if _needs_clarification(state):
        return "needs_clarification"
    if state.errors and not state.evidence:
        return "error"
    if not state.evidence:
        return "no_data"
    if not state.review_passed:
        return "error"
    # Что-то собрать удалось, но не всё, о чём просили. Называть это успехом нельзя:
    # пользователь должен видеть, что часть запроса осталась невыполненной.
    if state.errors:
        return "partial"
    planned = list(state.observations.get("scenario_ids") or [])
    completed = set(state.observations.get("completed_steps") or [])
    if planned and not completed.issuperset(planned):
        return "partial"
    return "success"


def run_workflow(request: ResearchRequest) -> tuple[ResearchResult, WorkflowState]:
    state = WorkflowState(request=request)
    try:
        state = planner.plan(state)
    except RuntimeError as exc:
        state.add_error(str(exc))
        return _finalize(state)

    if _needs_clarification(state):
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

    if _needs_clarification(state):
        return _finalize(state)

    _critic_replan(state)

    # Пересказ создаётся ДО проверки: иначе review_agent проверяет пустой текст
    # и любая ошибка LLM в числах уходит пользователю незамеченной.
    state.explanation = _build_explanation(state)
    state = review_agent.review(state)

    if not state.review_passed and state.observations.get("explanation_source") != "deterministic":
        # Факты собраны верно, подвёл только пересказ: заменяем его на перечисление
        # измерений и проверяем заново, вместо того чтобы хоронить весь ответ.
        state.limitations.append(
            "Пересказ LLM не прошёл проверку Evidence и заменён на перечисление измерений: "
            + "; ".join(state.review_errors)
        )
        state.errors = [err for err in state.errors if not err.startswith("Review:")]
        state.observations["explanation_source"] = "deterministic"
        state.explanation = _deterministic_explanation(state)
        state = review_agent.review(state)

    return _finalize(state)


def _finalize(state: WorkflowState) -> tuple[ResearchResult, WorkflowState]:
    status = _result_status(state)
    if status == "needs_clarification" and not state.errors:
        state.add_error("Недостаточно входных данных для сценария. Уточните запрос.")
    result = ResearchResult(
        status=status,
        evidence=list(state.evidence),
        explanation=state.explanation,
        limitations=list(state.limitations),
        errors=list(state.errors),
    )
    return result, state