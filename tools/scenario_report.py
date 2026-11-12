"""Сборка описания сценариев из фактического прогона.

Кейс требует показать «запрос, ожидаемый ответ, действия системы и фактический
результат». Единственный источник правды здесь — сам прогон: и demo_requests.json,
и docs/scenarios.md генерируются из него, поэтому разойтись с кодом не могут.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from models.result import ResearchResult
from orchestration.scenarios import SCENARIOS
from orchestration.state import WorkflowState

_ROOT = Path(__file__).resolve().parent.parent


def _agent_chain(state: WorkflowState) -> str:
    """Последовательность агентов без повторов подряд — «действия системы»."""
    chain: list[str] = []
    for row in state.agent_trace:
        agent = str(row.get("agent"))
        if not chain or chain[-1] != agent:
            chain.append(agent)
    return " → ".join(chain)


def _key_evidence(result: ResearchResult, limit: int = 6) -> list[str]:
    seen: set[str] = set()
    lines: list[str] = []
    for ev in result.evidence:
        method = str(ev.method or "value")
        if method in seen:
            continue
        seen.add(method)
        value = ev.value
        if isinstance(value, float):
            value = round(value, 3)
        text = str(value)
        if len(text) > 42:
            text = text[:42] + "…"
        unit = f" {ev.unit}" if ev.unit else ""
        lines.append(f"`{method}` = {text}{unit} [{ev.kind}, {ev.source or '—'}]")
        if len(lines) >= limit:
            break
    return lines


def collect(sid: str, result: ResearchResult, state: WorkflowState) -> dict[str, Any]:
    spec = SCENARIOS[sid]
    auto = state.observations.get("auto_steps") or []
    return {
        "id": sid,
        "title": spec.title,
        "request": spec.demo_request,
        "expected": spec.expected,
        "shows": spec.shows,
        "status": result.status,
        "planner": state.planner_source,
        "steps": state.plan,
        "auto_steps": [f"{s['added']} (не хватало {s['missing']})" for s in auto],
        "agents": _agent_chain(state),
        "evidence": _key_evidence(result),
        "artifacts": sorted(state.artifacts),
        "limitations": list(result.limitations),
        "errors": list(result.errors),
    }


def write_demo_requests(runs: list[dict[str, Any]], path: Path | None = None) -> Path:
    path = path or _ROOT / "data" / "demo_requests.json"
    payload = {
        "description": "Демо-запросы кейса СИНТЕЗ №1. Файл генерируется командой `python cli.py demo`.",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scenarios": [
            {
                "id": run["id"],
                "request": run["request"],
                "expected": run["expected"],
                "actual_status": run["status"],
            }
            for run in runs
        ],
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_scenarios_md(runs: list[dict[str, Any]], path: Path | None = None) -> Path:
    path = path or _ROOT / "docs" / "scenarios.md"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    out: list[str] = [
        "# Сценарии DEARMO (кейс СИНТЕЗ №1)",
        "",
        f"Файл сгенерирован командой `python cli.py demo` {stamp}. "
        "Колонка «фактический результат» — вывод этого прогона, а не описание намерений.",
        "",
        "| Сценарий | Статус | Шаги | Чем интересен |",
        "| --- | --- | --- | --- |",
    ]
    for run in runs:
        out.append(
            f"| {run['title']} | `{run['status']}` | {', '.join(run['steps']) or '—'} | {run['shows'] or '—'} |"
        )
    out.append("")

    for run in runs:
        out.extend(
            [
                f"## {run['title']} (`{run['id']}`)",
                "",
                f"**Запрос:** {run['request']}",
                "",
                f"**Ожидаемый ответ:** {run['expected'] or '—'}",
                "",
                f"**Действия системы:** {run['agents'] or '—'}",
                "",
                f"**Планировщик:** `{run['planner'] or '—'}`",
                "",
            ]
        )
        if run["auto_steps"]:
            out.extend(
                [
                    "**Шаги, добавленные системой самостоятельно:** "
                    + "; ".join(run["auto_steps"]),
                    "",
                ]
            )
        out.extend([f"**Фактический результат:** `{run['status']}`", ""])
        if run["evidence"]:
            out.append("Ключевые значения:")
            out.append("")
            out.extend(f"- {line}" for line in run["evidence"])
            out.append("")
        if run["artifacts"]:
            out.extend(["Артефакты: " + ", ".join(f"`{name}`" for name in run["artifacts"]), ""])
        if run["errors"]:
            out.append("Ошибки:")
            out.append("")
            out.extend(f"- {err}" for err in run["errors"])
            out.append("")
        if run["limitations"]:
            out.append("Ограничения, которые система назвала сама:")
            out.append("")
            out.extend(f"- {lim}" for lim in run["limitations"])
            out.append("")

    path.write_text("\n".join(out), encoding="utf-8")
    return path
