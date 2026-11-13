"""Терминальный интерфейс DRUGDIS с отображением статуса кэша, экспортом CSV/JSON, сравнением кандидатов и интерактивным уточнением."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Callable

from models.request import ResearchRequest
from models.result import ResearchResult
from orchestration.scenarios import SCENARIOS
from orchestration.state import WorkflowState
from orchestration.workflow import run_workflow
from tools.export import export_candidates_csv, export_run_bundle
from tools import scenario_report
from tools.ranking import build_candidates_comparison_table


def _format_cache_badge(context: dict) -> str:
    if not isinstance(context, dict):
        return "[LIVE]"
    is_cached = bool(
        context.get("is_cached")
        or context.get("_from_cache")
        or context.get("_from_download")
    )
    if not is_cached:
        return "[LIVE]"
    fetched_at = context.get("fetched_at") or "local-cache"
    return f"[CACHE: {fetched_at}]"


def _print_result_objects(result: ResearchResult, state: WorkflowState) -> None:
    print(f"status: {result.status}")
    print(f"planner: {state.planner_source}")
    print("plan: " + ", ".join(state.plan))
    print("completed: " + ", ".join(state.observations.get("completed_steps") or []))
    for item in state.observations.get("tool_plans") or []:
        print(f"tools {item.get('scenario')} ({item.get('source')}): " + ", ".join(item.get("tools") or []))
    if state.observations.get("pipeline_phase"):
        print(f"phase: {state.observations['pipeline_phase']}")
    if state.work_queue:
        print("queue: " + ", ".join(f"{item.pocket_id}:{item.status}" for item in state.work_queue))
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
    print(f"{'kind':<14} {'method':<28} {'value':<18} {'source':<14} cache_status")
    for ev in result.evidence:
        value = str(ev.value)
        if len(value) > 16:
            value = value[:15] + "…"
        ctx = ev.context if isinstance(ev.context, dict) else {}
        badge = _format_cache_badge(ctx)
        print(f"{ev.kind:<14} {(ev.method or ''):<28} {value:<18} {(ev.source or ''):<14} {badge}")
    if state.artifacts:
        print()
        print("artifacts:")
        for key, path in state.artifacts.items():
            print(f"  {key}: {path}")
    print()
    print("trace:")
    for row in state.agent_trace:
        dur = f" ({row['duration_s']}s)" if row.get("duration_s") is not None else ""
        print(f"  {row.get('step'):>2} {row.get('agent'):<16} {row.get('action')}{dur}")


def _print_result(
    text: str,
    *,
    export_dir: str | None = None,
    as_json: bool = False,
    interactive: bool = False,
    input_fn: Callable[[str], str] = input,
) -> tuple[ResearchResult, WorkflowState]:
    req = ResearchRequest(text=text)
    result, state = run_workflow(req)

    # Интерактивный цикл доуточнения при неоднозначности (needs_clarification)
    if interactive and result.status == "needs_clarification":
        print("\nТребуется уточнение запроса:")
        for err in result.errors:
            print(f"  * {err}")
        candidates = state.observations.get("clarification_candidates") or []
        if candidates:
            print("Доступные кандидаты:")
            for idx, cand in enumerate(candidates, start=1):
                print(
                    f"  {idx}) {cand.get('title')} (CID: {cand.get('cid')}, SMILES: {cand.get('connectivity_smiles')})"
                )
        try:
            clarification = input_fn("Введите точное название, SMILES или номер кандидата: ").strip()
        except EOFError:
            clarification = ""

        if clarification:
            chosen_mol = clarification
            if clarification.isdigit() and candidates:
                c_idx = int(clarification) - 1
                if 0 <= c_idx < len(candidates):
                    chosen_mol = (
                        candidates[c_idx].get("connectivity_smiles")
                        or candidates[c_idx].get("title")
                        or clarification
                    )
            updated_req = ResearchRequest(
                text=f"{text} (уточнение: {chosen_mol})",
                molecules=[chosen_mol],
                ligand_ids=[chosen_mol],
                pdb_id=state.request.pdb_id,
                constraints=dict(state.request.constraints),
            )
            result, state = run_workflow(updated_req)

    if export_dir:
        export_run_bundle(result, state, export_dir=export_dir)

    if as_json:
        payload = {
            "status": result.status,
            "scenario_id": state.scenario_id,
            "explanation": result.explanation,
            "limitations": result.limitations,
            "errors": result.errors,
            "artifacts": state.artifacts,
            "evidence_count": len(result.evidence),
        }
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        _print_result_objects(result, state)

    return result, state


def _run_compare(molecules: list[str], target_id: str, export_dir: str | None) -> int:
    table = build_candidates_comparison_table(molecules, target_id=target_id)
    print(f"Сравнение кандидатов для мишени {target_id}:")
    print(
        f"{'Rank':<5} {'Molecule':<18} {'MW':<8} {'logP':<7} {'TPSA':<7} {'QED':<7} {'Ki(nM)':<10} {'Status':<10} Failed Reasons"
    )
    for row in table:
        reasons = "; ".join(row.get("failed_reasons") or []) or "—"
        print(
            f"{row['rank']:<5} {str(row['molecule'])[:17]:<18} "
            f"{str(row['mw'] if row['mw'] is not None else '—'):<8} "
            f"{str(row['logp'] if row['logp'] is not None else '—'):<7} "
            f"{str(row['tpsa'] if row['tpsa'] is not None else '—'):<7} "
            f"{str(row['qed'] if row['qed'] is not None else '—'):<7} "
            f"{str(row['ki_nm'] if row['ki_nm'] is not None else 'None'):<10} "
            f"{row['filter_status']:<10} {reasons}"
        )
    if export_dir:
        out_path = export_candidates_csv(table, Path(export_dir) / "candidates_comparison.csv")
        print(f"\nТаблица сравнения сохранена в: {out_path}")
    return 0


def _rel(value):
    """Путь относительно репозитория: в документации не должно быть домашних каталогов."""
    if isinstance(value, (list, tuple)):
        return [_rel(item) for item in value]
    text = str(value)
    root = str(Path(__file__).resolve().parent) + "/"
    return text[len(root):] if text.startswith(root) else text


def _carry_context(previous: dict, state: WorkflowState) -> dict:
    """Что следующий сценарий цепочки получает от предыдущего.

    Передаются только факты, которые предыдущий шаг действительно установил:
    молекула, структура, таблица поз и поза №1. Ничего не выдумывается.
    """
    obs = state.observations
    carried = dict(previous)
    ligand = obs.get("ligand_query")
    if ligand and not carried.get("molecules"):
        carried["molecules"] = [str(ligand)]
    if obs.get("pdb_id"):
        carried["pdb_id"] = str(obs["pdb_id"])
    poses = obs.get("poses_csv_path") or state.artifacts.get("poses_csv")
    if poses:
        carried["poses_csv_path"] = str(poses)
    pose = obs.get("best_pose") or state.artifacts.get("pose_pdbqt")
    if pose:
        carried["pose_file"] = str(pose)
    return carried


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="drugdis", description="DRUGDIS — мультиагентная система из терминала")
    sub = parser.add_subparsers(dest="cmd", required=True)

    ask = sub.add_parser("ask", help="Выполнить исследовательский запрос")
    ask.add_argument("text", nargs="+", help="Текст запроса")
    ask.add_argument("--export-dir", default=None, help="Папка для сохранения CSV и JSON отчётов")
    ask.add_argument("--json", action="store_true", help="Вывести краткий результат в формате JSON")
    ask.add_argument(
        "--interactive",
        action="store_true",
        help="Запрашивать уточнение в терминале при неоднозначности (needs_clarification)",
    )

    demo = sub.add_parser("demo", help="Прогнать все 4 демонстрационных сценария и сохранить контрольные отчёты")
    demo.add_argument("--export-dir", default="outputs/demo", help="Папка для сохранения контрольных результатов")
    demo.add_argument("--no-docs", action="store_true", help="Не обновлять demo_requests.json и docs/scenarios.md")

    comp = sub.add_parser("compare", help="Построить сводную таблицу сравнения нескольких молекул-кандидатов")
    comp.add_argument("molecules", nargs="+", help="Названия или SMILES молекул для сравнения")
    comp.add_argument("--target", default="", help="Мишень ChEMBL (ID или название). По умолчанию — активности по всем мишеням молекулы")
    comp.add_argument("--export-dir", default="outputs/exports", help="Папка для сохранения CSV таблицы")

    sub.add_parser("fetch-data", help="Скачать опорные PDB, BindingDB и PubChem в локальный кэш")

    args = parser.parse_args(argv)
    if args.cmd == "fetch-data":
        from tools.dataset_fetch import fetch_all

        for key, path in fetch_all().items():
            print(f"{key}: {path}")
        return 0

    if args.cmd == "demo":
        runs = []
        # Сценарии кейса — одна цепочка на одной молекуле: мишень -> докинг -> сравнение
        # поз -> проверка позы. Поэтому результат шага передаётся следующему, как это
        # делал бы пользователь в одном диалоге. Передача видна в выводе и в отчёте.
        carried: dict = {}
        for sid, spec in SCENARIOS.items():
            if not spec.demo_request:
                continue
            print(f"\n=== Сценарий: {sid} ({spec.title}) ===")
            if carried:
                print("контекст из предыдущего шага: " + ", ".join(f"{k}={_rel(v)}" for k, v in carried.items()))
            res, st = run_workflow(ResearchRequest(text=spec.demo_request, **carried))
            export_run_bundle(res, st, export_dir=args.export_dir, prefix=f"demo_{sid}")
            _print_result_objects(res, st)
            run = scenario_report.collect(sid, res, st)
            run["carried_context"] = {k: _rel(v) for k, v in carried.items()}
            runs.append(run)
            carried = _carry_context(carried, st)
        if runs and not args.no_docs:
            # Описание сценариев пишется из фактического прогона, а не правится руками.
            print("\nОбновлены файлы описания сценариев:")
            print(" ", scenario_report.write_demo_requests(runs))
            print(" ", scenario_report.write_scenarios_md(runs))
        failed = [run["id"] for run in runs if run["status"] not in ("success", "partial")]
        if failed:
            print("\nСценарии без результата: " + ", ".join(failed))
            return 1
        return 0

    if args.cmd == "compare":
        return _run_compare(args.molecules, args.target, args.export_dir)

    _print_result(
        " ".join(args.text),
        export_dir=args.export_dir,
        as_json=args.json,
        interactive=args.interactive,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
