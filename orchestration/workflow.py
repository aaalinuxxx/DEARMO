"""Управление порядком агентов, состоянием, ограничениями шагов и обработкой ошибок."""
from models.request import ResearchRequest
from models.result import ResearchResult

def run_workflow(request: ResearchRequest) -> ResearchResult:
    """TODO: определить маршруты сценариев, лимиты, повторы и вызовы агентов.

    При неполных входных данных запрашивать уточнение; не подменять ошибки данными.
    """
    raise NotImplementedError("Оркестрация ещё не реализована")
