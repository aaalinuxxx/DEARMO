"""Запрос исследователя и предварительный формат извлечённых параметров с валидацией входов."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")

_NUMERIC_CONSTRAINT_KEYS = {
    "mw_max",
    "mw_min",
    "logp_max",
    "logp_min",
    "tpsa_max",
    "tpsa_min",
    "hbd_max",
    "hba_max",
    "rotatable_bonds_max",
    "qed_min",
    "ki_max_nm",
    "ic50_max_nm",
}

_MIN_MAX_PAIRS = (
    ("mw_min", "mw_max"),
    ("logp_min", "logp_max"),
    ("tpsa_min", "tpsa_max"),
)


@dataclass
class ResearchRequest:
    """Исходный текст обязателен; остальные поля заполняются и проверяются при инициализации."""

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

        if self.pdb_id is not None:
            if not isinstance(self.pdb_id, str) or not _PDB_ID_RE.match(self.pdb_id.strip()):
                raise ValueError(
                    f"Некорректный формат pdb_id '{self.pdb_id}': требуется ровно 4 алфавитно-цифровых символа (например, 3HS4)"
                )
            self.pdb_id = self.pdb_id.strip().upper()

        if self.constraints is not None:
            if not isinstance(self.constraints, dict):
                raise ValueError("Поле constraints должно быть словарём")
            for key, val in self.constraints.items():
                if key in _NUMERIC_CONSTRAINT_KEYS and val is not None:
                    if isinstance(val, bool) or not isinstance(val, (int, float)):
                        raise ValueError(f"Ограничение '{key}' должно быть числом, получено: {val!r}")
            for min_key, max_key in _MIN_MAX_PAIRS:
                if (
                    min_key in self.constraints
                    and max_key in self.constraints
                    and self.constraints[min_key] is not None
                    and self.constraints[max_key] is not None
                ):
                    if float(self.constraints[min_key]) > float(self.constraints[max_key]):
                        raise ValueError(
                            f"Некорректный диапазон ограничений: {min_key} ({self.constraints[min_key]}) > {max_key} ({self.constraints[max_key]})"
                        )
