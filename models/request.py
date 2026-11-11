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
    ligand_ids: list[str] = field(default_factory=list)
    pdb_id: str | None = None
    poses_csv_path: str | None = None
    pose_file: str | None = None
    constraints: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.text, str) or not self.text.strip():
            raise ValueError("Текст запроса не должен быть пустым")
        self.text = self.text.strip()
