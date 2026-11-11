"""Агент проверки: сопоставляет выводы с фактами и выявляет неподтверждённые утверждения."""
from orchestration.state import WorkflowState

def review(state: WorkflowState) -> WorkflowState:
    """TODO: согласовать проверки, случаи возврата на исправление и лимит повторов.

    Проверка LLM дополняет программную проверку чисел и источников, а не заменяет её.
    """
    raise NotImplementedError("Агент проверки ещё не реализован")
