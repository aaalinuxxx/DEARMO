"""Терминальный интерфейс DEARMO."""
from __future__ import annotations

import argparse
import sys

from models.request import ResearchRequest
from orchestration.workflow import run_workflow


def _print_result(text: str) -> None:
    result, state = run_workflow(ResearchRequest(text=text))
    print(f"status: {result.status}")
    print(f"planner: {state.planner_source}")
    print("plan: " + ", ".join(state.plan))
    print("completed: " + ", ".join(state.observations.get("completed_steps") or []))
    if result.explanation:
        print()
        print(result.explanation)
    if result.limitations:
        print()
        print("limitations:")
        for item in result.limitations:
            print(f"  - {item}")
    if result.errors:
        print()
        print("errors:")
        for item in result.errors:
            print(f"  - {item}")
    print()
    print(f"{'kind':<14} {'method':<28} {'value':<18} source")
    for ev in result.evidence:
        value = str(ev.value)
        if len(value) > 16:
            value = value[:15] + "…"
        print(f"{ev.kind:<14} {(ev.method or ''):<28} {value:<18} {ev.source or ''}")
    if state.artifacts:
        print()
        print("artifacts:")
        for key, path in state.artifacts.items():
            print(f"  {key}: {path}")
    print()
    print("trace:")
    for row in state.agent_trace:
        print(f"  {row.get('step'):>2} {row.get('agent'):<16} {row.get('action')}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dearmo", description="DEARMO — запросы из терминала")
    sub = parser.add_subparsers(dest="cmd", required=True)

    ask = sub.add_parser("ask", help="Выполнить исследовательский запрос")
    ask.add_argument("text", nargs="+", help="Текст запроса")

    sub.add_parser("fetch-data", help="Скачать PDB, BindingDB и PubChem для ацетазоламида")

    args = parser.parse_args(argv)
    if args.cmd == "fetch-data":
        from tools.dataset_fetch import fetch_all

        for key, path in fetch_all().items():
            print(f"{key}: {path}")
        return 0
    _print_result(" ".join(args.text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
