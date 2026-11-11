"""Агент данных: получает структуры и экспериментальные записи через научные инструменты."""
from orchestration.state import WorkflowState

def collect_data(state: WorkflowState) -> WorkflowState:
    """TODO: выбрать обращения к PubChem/ChEMBL, сохранить источники и ошибки."""
    raise NotImplementedError("Агент данных ещё не реализован")
