"""Общие модели результата и происхождения исследовательских данных."""
from dataclasses import dataclass, field
from typing import Any, Literal

# TODO (IT 2): согласовать состояния частичного успеха и формат ошибок API.
ResultStatus = Literal["success", "error", "no_data", "needs_clarification"]
EvidenceKind = Literal["database", "calculation", "prediction", "hypothesis"]

@dataclass
class Evidence:
    """Одно значение с указанием его происхождения и контекста."""
    kind: EvidenceKind
    value: Any
    unit: str | None = None
    source: str | None = None
    method: str | None = None
    context: dict[str, Any] = field(default_factory=dict)

@dataclass
class ResearchResult:
    """Структурированный ответ, отдельно от ошибок и ограничений."""
    status: ResultStatus
    evidence: list[Evidence] = field(default_factory=list)
    explanation: str | None = None
    limitations: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    # TODO (IT 2 + биологи): утвердить поля assay и связи утверждений с Evidence.
