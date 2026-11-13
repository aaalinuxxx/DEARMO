"""Маршрутизация сценариев, лимиты шагов, ранний выход на уточнение и сборка ResearchResult."""
from __future__ import annotations

from agents import planner, review_agent
from agents.target_agent import split_variant
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
    # Тот же порядок, что и у пересказа модели: ответ по существу (мишень, accession,
    # PDB, цепь) раньше дескрипторов. Простая обрезка первых 12 фактов теряла его.
    lines.extend(_summary_lines(state))
    for extra in (_caveat_line(state), _sources_line(state)):
        if extra:
            lines.append(extra)
    return "\n".join(lines)


# Факты, без которых ответ не отвечает на вопрос: они идут в пересказ первыми.
# Простая обрезка первых 12 фактов оставляла за бортом мишень, accession и PDB,
# потому что дескрипторы RDKit считаются раньше и занимали весь лимит.
_ANSWER_FIRST = (
    "pubchem_cid",
    "uniprot_accession",
    # Измерения — экспериментальное подтверждение выбора, ядро ответа. Раньше они
    # шли после дескрипторов и обрезались лимитом: в тексте было «константа не указана».
    "Ki",
    "Kd",
    "IC50",
    "EC50",
    "target_choice_rationale",
    "pdb_structure",
    "structure_chain",
    "resolution",
    "metal_ion",
    "bound_ligand",
    "structure_search_tier",
    "vina_score",
    "vina_delta_to_best_pose",
    "contact_gained_in_pose",
    "contact_lost_in_pose",
    "key_distance_by_pose",
    "pose_to_inspect",
    "geometry_verdict",
    "checks_passed_percent",
)

# Дескрипторы относятся к молекуле, а не к белку, и это надо говорить вслух.
_LIGAND_PROPERTIES = frozenset(
    {
        "molecular_weight",
        "logp",
        "tpsa",
        "qed",
        "hbd",
        "hba",
        "rotatable_bonds",
        "heavy_atoms",
        "sa_score",
        "fraction_csp3",
        "rings",
    }
)


def _summary_lines(state: WorkflowState, limit: int = 16) -> list[str]:
    """По одному представителю каждого вида факта, ответ по существу — первым."""
    by_method: dict[str, list] = {}
    for ev in state.evidence:
        by_method.setdefault(str(ev.method or "value"), []).append(ev)

    order = [m for m in _ANSWER_FIRST if m in by_method]
    order += [m for m in by_method if m not in order]

    lines: list[str] = []
    for method in order:
        for ev in by_method[method][:2 if method in _ANSWER_FIRST else 1]:
            context = ev.context or {}
            if method == "structure_search_tier":
                # Код уровня сам по себе ничего не говорит: в пересказ идёт смысл.
                lines.append(f"- [{ev.kind}] найденная структура — {context.get('meaning')} ({ev.source or '—'})")
                continue
            if method in _LIGAND_PROPERTIES:
                # Без явной привязки Qwen писала «молекулярная масса белка
                # составляет 222.251 г/моль»: это масса лиганда, а не белка.
                ligand = state.observations.get("ligand_query") or "лиганд"
                suffix = f" [свойство лиганда {ligand}]"
            else:
                # Название соединения в сводку не идёт: в «[6,7-BIS(2-METHOXY…]» ревью
                # находило числа 6 и 7, которых в Evidence нет. Для кода компонента
                # достаточно comp_id, для иона — его имя («ZINC ION»).
                if method == "bound_ligand":
                    target = None
                elif method == "metal_ion":
                    target = context.get("name")
                else:
                    target = context.get("target_name") or context.get("name")
                if target:
                    # «[A65S,N67Q,…]» в сводке модель теряла или пересказывала с искажениями.
                    base, variant = split_variant(str(target))
                    target = f"{base}, мутантный вариант" if variant else base
                suffix = f" [{target}]" if target and method not in ("pdb_structure",) else ""
            unit = f" {ev.unit}" if ev.unit else ""
            lines.append(f"- [{ev.kind}] {method}={ev.value}{unit}{suffix} ({ev.source or '—'})")
            if len(lines) >= limit:
                return lines
    return lines


def _caveat_line(state: WorkflowState) -> str:
    """Оговорка, которую нельзя оставлять только в списке ограничений.

    Если измерение относится к мутанту, а не к природному белку, это меняет смысл
    числа рядом с ним. Модель такую оговорку в пересказ не переносит, поэтому её
    добавляет код.
    """
    for lim in state.limitations:
        if "мутантным вариантам" in lim:
            return "Важно: " + lim
    return ""


def _sources_line(state: WorkflowState) -> str:
    """Ссылки на записи BindingDB. Собираются кодом: модели URL не доверяем."""
    links: list[str] = []
    for ev in state.evidence:
        link = (ev.context or {}).get("link")
        if ev.source == "BindingDB" and link and link not in links:
            links.append(str(link))
    if not links:
        return ""
    return "Записи измерений BindingDB: " + "; ".join(links[:3])


def _build_explanation(state: WorkflowState) -> str:
    lines = _summary_lines(state)
    summary = "\n".join(lines) if lines else "Данных для объяснения нет."
    narrated = narrate_with_llm(summary, state.limitations)
    if narrated:
        state.observations["explanation_source"] = "llm"
        extras = [part for part in (_caveat_line(state), _sources_line(state)) if part]
        return "\n\n".join([narrated, *extras])
    # Текст собран программой из фактов дословно. Без этой пометки ревью проверяло
    # его как сочинение модели и писало ложное «пересказ LLM не прошёл проверку»,
    # хотя LLM в этом прогоне вообще не участвовала.
    state.observations["explanation_source"] = "deterministic"
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
        # Замечание небольшой модели без названного следующего шага не проверено кодом:
        # на практике Qwen называла неполным ответ, где докинг в каждом кармане был
        # выполнен. Поэтому оно подаётся как мнение критика, а не как установленный факт.
        state.limitations.append(
            f"Замечание LLM-критика, не подтверждённое проверкой: {missing}. "
            "Шаг для устранения критик не назвал, автоматически это не достраивается."
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