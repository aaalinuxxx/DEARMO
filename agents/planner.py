"""Планировщик: определяет сценарий, сущности, ограничения и последовательность действий."""
from orchestration.state import WorkflowState

def plan(state: WorkflowState) -> WorkflowState:
    """TODO: подключить LLM и проверить структурированный план; уточнять неоднозначность."""
    raise NotImplementedError("Планировщик ещё не реализован")
