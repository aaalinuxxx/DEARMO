"""Разбор составного запроса в упорядоченный план шагов и повтор при сбое."""
from __future__ import annotations

from dataclasses import dataclass

from orchestration.scenarios import SCENARIO_IDS, ScenarioId


@dataclass
class PlanStep:
    scenario_id: ScenarioId
    reason: str
    attempts: int = 0


def order_scenarios(ids: list[str]) -> list[str]:
    known = [sid for sid in SCENARIO_IDS if sid in ids]
    extra = [sid for sid in ids if sid not in known]
    return known + extra


def steps_from_ids(ids: list[str], *, source: str) -> list[PlanStep]:
    ordered = order_scenarios(ids)
    return [PlanStep(scenario_id=sid, reason=f"{source}:{sid}") for sid in ordered]
