"""Проверка базовых контрактов запроса и независимости состояний."""
import pytest
from models.request import ResearchRequest
from models.result import Evidence
from orchestration.state import WorkflowState


def test_empty_request_rejected():
    with pytest.raises(ValueError):
        ResearchRequest(text="  ")


def test_states_do_not_share_mutable_fields():
    a = WorkflowState(request=ResearchRequest(text="первый"))
    b = WorkflowState(request=ResearchRequest(text="второй"))
    a.add_evidence(Evidence(kind="calculation", value=1))
    a.add_error("Ошибка")
    assert b.evidence == []
    assert b.errors == []
