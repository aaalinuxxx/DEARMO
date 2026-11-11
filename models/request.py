"""Запрос исследователя и предварительный формат извлечённых параметров."""
from dataclasses import dataclass, field
from typing import Any

@dataclass
class ResearchRequest:
    """Исходный текст обязателен; остальные поля заполняются по мере обработки."""
    text: str
    scenario: str | None = None
    molecules: list[str] = field(default_factory=list)
    target_id: str | None = None
    # TODO (IT 2 + биологи): согласовать типы ограничений, единицы и допустимые значения.
    constraints: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.text, str) or not self.text.strip():
            raise ValueError("Текст запроса не должен быть пустым")
        self.text = self.text.strip()
