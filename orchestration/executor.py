"""Исполнение плана: предусловия, повтор шага, журнал."""
from __future__ import annotations

from orchestration.handlers import (
    handle_compare_docking,
    handle_dock_pockets,
    handle_target_for_ligand,
    handle_validate_pose,
)
from orchestration.scenarios import missing_required_keys, producer_for
from orchestration.state import WorkflowState
from orchestration.task_graph import PlanStep, steps_from_ids

HANDLERS = {
    "target_for_ligand": handle_target_for_ligand,
    "dock_pockets": handle_dock_pockets,
    "compare_docking": handle_compare_docking,
    "validate_pose": handle_validate_pose,
}


def _repair(state: WorkflowState, scenario_id: str) -> None:
    """Достраивает недостающие входы из уже собранных наблюдений, без выдуманных чисел."""
    rec = state.observations.get("recommended_target") or {}
    if scenario_id in ("dock_pockets", "validate_pose") and not state.observations.get("pdb_id"):
        pdb_id = rec.get("pdb_id") or state.request.pdb_id
        if pdb_id:
            state.observations["pdb_id"] = pdb_id
            state.request.pdb_id = pdb_id
    ligand = state.observations.get("ligand_query") or state.observations.get("ligand_id")
    if ligand and not state.observations.get("ligand_id"):
        state.observations["ligand_id"] = ligand
    if scenario_id == "validate_pose" and state.observations.get("best_pose"):
        state.observations["pose_file"] = state.observations["best_pose"]


_MAX_DEPENDENCY_DEPTH = 2


def _ensure_inputs(state: WorkflowState, scenario_id: str, completed: list[str], depth: int = 0) -> None:
    """Достроить недостающие входы сценария, выполнив шаги, которые их производят.

    Запрос «проверь позу X» не содержит ни структуры, ни позы. Вместо отказа
    система сама planирует докинг, а при необходимости и подбор мишени перед ним.
    Глубина ограничена, чтобы цепочка зависимостей не могла зациклиться.
    """
    if depth >= _MAX_DEPENDENCY_DEPTH:
        return
    _repair(state, scenario_id)
    for key in missing_required_keys(state.observations, scenario_id):
        producer = producer_for(key)
        if not producer or producer == scenario_id or producer in completed:
            continue
        state.trace(
            "orchestrator",
            f"resolve:{key}",
            detail=f"для {scenario_id} не хватает {key} — добавлен шаг {producer}",
        )
        state.observations.setdefault("auto_steps", []).append(
            {"for": scenario_id, "missing": key, "added": producer}
        )
        _ensure_inputs(state, producer, completed, depth + 1)
        if producer not in completed:
            _run_step(state, PlanStep(scenario_id=producer, reason=f"dependency:{key}"), completed)


def execute_plan(state: WorkflowState) -> WorkflowState:
    ids = state.observations.get("scenario_ids") or ([state.scenario_id] if state.scenario_id else [])
    steps = steps_from_ids(ids, source=state.planner_source or "plan")
    completed: list[str] = []
    for step in steps:
        if step.scenario_id in completed:
            continue
        _ensure_inputs(state, step.scenario_id, completed)
        if step.scenario_id in completed:
            continue
        _run_step(state, step, completed)
    state.observations["completed_steps"] = completed
    # Порядок в plan — фактический порядок выполнения, включая шаги,
    # которые система добавила сама: иначе в выводе докинг стоял бы после
    # проверки позы, ради которой он и был запущен.
    planned = [s.scenario_id for s in steps]
    state.plan = [f"ok:{sid}" for sid in completed]
    state.plan += [f"fail:{sid}" for sid in planned if sid not in completed]
    return state


def _run_step(state: WorkflowState, step: PlanStep, completed: list[str]) -> None:
    handler = HANDLERS.get(step.scenario_id)
    if handler is None:
        state.add_error(f"Нет обработчика для {step.scenario_id}")
        return
    last_error: Exception | None = None
    for attempt in range(2):
        step.attempts = attempt + 1
        _repair(state, step.scenario_id)
        before = len(state.errors)
        try:
            with state.agent_step("orchestrator", f"run:{step.scenario_id}#{attempt + 1}"):
                handler(state)
            new_errors = state.errors[before:]
            if new_errors and attempt == 0:
                state.limitations.append(f"Повтор шага {step.scenario_id}: {new_errors[-1]}")
                continue
            completed.append(step.scenario_id)
            return
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            state.limitations.append(f"Сбой {step.scenario_id}, попытка {attempt + 1}: {exc}")
    if last_error is not None:
        state.add_error(f"Шаг {step.scenario_id} не выполнен: {last_error}")
