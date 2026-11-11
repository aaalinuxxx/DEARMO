"""Независимое состояние обработки одного исследовательского запроса."""
from dataclasses import dataclass, field
from typing import Any
from models.request import ResearchRequest
from models.result import Evidence

@dataclass
class WorkflowState:
    request: ResearchRequest
    plan: list[str] = field(default_factory=list)
    evidence: list[Evidence] = field(default_factory=list)
    observations: dict[str, Any] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    step_count: int = 0

    def add_evidence(self, evidence: Evidence) -> None:
        self.evidence.append(evidence)

    def add_error(self, message: str) -> None:
        if not message.strip():
            raise ValueError("Сообщение об ошибке не должно быть пустым")
        self.errors.append(message)

    # TODO (IT 1 + IT 2): согласовать маршруты агентов, лимиты и хранение истории.
