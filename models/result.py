"""Общие модели результата и происхождения исследовательских данных."""
from dataclasses import dataclass, field
from typing import Any, Literal

# success — всё запрошенное выполнено; partial — часть шагов не дала результата;
# needs_clarification — запрос неоднозначен; no_data — источники ничего не вернули;
# error — расчёт или проверка не состоялись.
ResultStatus = Literal["success", "partial", "error", "no_data", "needs_clarification"]
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
