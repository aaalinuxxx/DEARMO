"""Экспорт результатов исследования (таблиц Evidence, сводной таблицы кандидатов и полного отчёта) в CSV и JSON."""
from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from models.result import ResearchResult
from orchestration.state import WorkflowState


def _now_compact() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _extract_cache_meta(context: dict[str, Any]) -> tuple[bool, str]:
    is_cached = bool(
        context.get("is_cached")
        or context.get("_from_cache")
        or context.get("_from_download")
    )
    fetched_at = str(context.get("fetched_at") or ("local-cache" if is_cached else "live"))
    return is_cached, fetched_at


def export_to_csv(
    result: ResearchResult,
    path: Path | str,
    *,
    state: WorkflowState | None = None,
) -> Path:
    """Сохранить таблицу доказательств (Evidence) в CSV с указанием источника и метаданных кэша."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "kind",
        "method",
        "value",
        "unit",
        "source",
        "is_cached",
        "fetched_at",
        "context_json",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for ev in result.evidence:
            ctx = ev.context if isinstance(ev.context, dict) else {}
            is_cached, fetched_at = _extract_cache_meta(ctx)
            writer.writerow(
                {
                    "kind": ev.kind,
                    "method": ev.method or "",
                    "value": "" if ev.value is None else ev.value,
                    "unit": ev.unit or "",
                    "source": ev.source or "",
                    "is_cached": is_cached,
                    "fetched_at": fetched_at,
                    "context_json": json.dumps(ctx, ensure_ascii=False, default=str),
                }
            )
    return out_path


def export_candidates_csv(
    candidates_table: list[dict[str, Any]],
    path: Path | str,
) -> Path:
    """Сохранить сводную таблицу сравнения кандидатов (с причинами исключения) в CSV."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "rank",
        "molecule",
        "smiles",
        "mw",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "qed",
        "ki_nm",
        "vina_score",
        "filter_status",
        "failed_reasons",
        "comparability_warnings",
        "evidence_sources",
        "missing_data",
        "preliminary_choice",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in candidates_table:
            serialized = dict(row)
            for list_key in ("failed_reasons", "comparability_warnings", "evidence_sources", "missing_data"):
                val = serialized.get(list_key)
                if isinstance(val, list):
                    serialized[list_key] = "; ".join(str(x) for x in val)
            writer.writerow(serialized)
    return out_path


def export_to_json(
    result: ResearchResult,
    state: WorkflowState | None,
    path: Path | str,
) -> Path:
    """Сохранить полный структурированный отчёт ResearchResult и журнал выполнения agent_trace в JSON."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    evidence_payload = []
    for ev in result.evidence:
        ctx = ev.context if isinstance(ev.context, dict) else {}
        is_cached, fetched_at = _extract_cache_meta(ctx)
        evidence_payload.append(
            {
                "kind": ev.kind,
                "method": ev.method,
                "value": ev.value,
                "unit": ev.unit,
                "source": ev.source,
                "is_cached": is_cached,
                "fetched_at": fetched_at,
                "context": ctx,
            }
        )

    payload: dict[str, Any] = {
        "exported_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": result.status,
        "scenario_id": state.scenario_id if state else None,
        "planner_source": state.planner_source if state else None,
        "plan": list(state.plan) if state else [],
        "completed_steps": list(state.observations.get("completed_steps", [])) if state else [],
        "request": {
            "text": state.request.text,
            "molecules": state.request.molecules,
            "target_id": state.request.target_id,
            "pdb_id": state.request.pdb_id,
            "constraints": state.request.constraints,
        }
        if state
        else None,
        "explanation": result.explanation,
        "limitations": list(result.limitations),
        "errors": list(result.errors),
        "evidence": evidence_payload,
        "candidates_table": state.observations.get("candidates_table") if state else None,
        "artifacts": dict(state.artifacts) if state else {},
        "agent_trace": list(state.agent_trace) if state else [],
    }

    with out_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2, default=str)
    return out_path


def export_run_bundle(
    result: ResearchResult,
    state: WorkflowState,
    export_dir: Path | str = "outputs/exports",
    *,
    prefix: str | None = None,
) -> dict[str, str]:
    """Сохранить комплект файлов запуска (CSV доказательств, CSV кандидатов при наличии и полный JSON)."""
    out_dir = Path(export_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scen = re.sub(r"[^a-zA-Z0-9_-]", "_", state.scenario_id or "query")
    stem = prefix or f"{_now_compact()}_{scen}"

    csv_path = export_to_csv(result, out_dir / f"{stem}_evidence.csv", state=state)
    state.artifacts["evidence_csv"] = str(csv_path)

    candidates = state.observations.get("candidates_table")
    if isinstance(candidates, list) and candidates:
        cand_path = export_candidates_csv(candidates, out_dir / f"{stem}_candidates.csv")
        state.artifacts["candidates_csv"] = str(cand_path)

    json_path = export_to_json(result, state, out_dir / f"{stem}_report.json")
    state.artifacts["report_json"] = str(json_path)

    return dict(state.artifacts)
