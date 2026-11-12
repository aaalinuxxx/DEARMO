"""Writer собирает report.md только из уже посчитанных Evidence."""
from __future__ import annotations

import uuid
from pathlib import Path

from orchestration.state import WorkflowState

_ROOT = Path(__file__).resolve().parent.parent


def write_report(state: WorkflowState) -> WorkflowState:
    with state.agent_step("writer", "report_md"):
        run_id = state.observations.get("run_id") or uuid.uuid4().hex[:8]
        state.observations["run_id"] = run_id
        out_dir = _ROOT / "outputs" / "runs" / run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "report.md"
        lines = [
            "# ORC docking report",
            "",
            f"- run_id: `{run_id}`",
            f"- pdb: {state.observations.get('pdb_id') or '—'}",
            f"- ligand: {state.observations.get('ligand_id') or state.observations.get('ligand_query') or '—'}",
            f"- planner: {state.planner_source or '—'}",
            "",
            "## Очередь карманов",
            "",
            "| pocket | status | agent |",
            "| --- | --- | --- |",
        ]
        for item in state.work_queue:
            lines.append(f"| {item.pocket_id} | {item.status} | {item.agent} |")
        lines.extend(["", "## Evidence", ""])
        for ev in state.evidence:
            if "vina" in (ev.method or "") or ev.method == "pdb_structure":
                lines.append(f"- [{ev.kind}] {ev.method} = {ev.value} ({ev.source})")
        lines.extend(
            [
                "",
                "Оценка Vina — расчёт. Она не является экспериментальной Ki или IC50.",
                "",
            ]
        )
        path.write_text("\n".join(lines), encoding="utf-8")
        state.artifacts["report_md"] = str(path)
        state.observations["pipeline_phase"] = "finalize"
    return state
