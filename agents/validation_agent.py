"""Сценарий 4: агент проверки геометрии позы №1.

Агент не решает, прошла ли поза. Вердикт считает `tools.pose_validation` по
измерениям, а агент раскладывает готовый отчёт по Evidence, чтобы каждое
число в ответе можно было проследить до расчёта. Языковой модели остаётся
пересказ: в этом проекте она ошибается не числами, а их привязкой, и
«проверка пройдена» при проваленной проверке не поймала бы ни одна сверка
значений.
"""
from __future__ import annotations

import json
from pathlib import Path

from models.result import Evidence
from orchestration.state import WorkflowState
from tools import pose_validation
from tools.pose_geometry import FAILED, GeometryReport
from tools.pose_validation import PoseInputs, load_run_inputs, render_markdown

_REPO_ROOT = Path(__file__).resolve().parent.parent

#: Сколько измерений одной проверки попадает в Evidence. Нарушения идут
#: первыми, поэтому обрезается хвост пройденных, а не то, что объясняет отказ.
_MAX_MEASUREMENTS_PER_CHECK = 5


def _as_path(value: object) -> Path | None:
    if not value:
        return None
    path = Path(str(value))
    if not path.is_absolute():
        path = _REPO_ROOT / path
    return path if path.is_file() else None


def _run_dir(state: WorkflowState) -> Path | None:
    """Каталог прогона сценария 2. Перебором каталогов не ищется.

    Подстановка чужого прогона — ровно то, что кейс запрещает: проверять
    надо позу, посчитанную по этому же запросу.
    """
    for candidate in (
        state.observations.get("poses_csv_path"),
        state.artifacts.get("poses_csv"),
    ):
        if candidate:
            parent = Path(str(candidate)).parent
            if (parent / "run_meta.json").is_file():
                return parent
    pose_file = _as_path(state.observations.get("pose_file") or state.request.pose_file)
    if pose_file is not None:
        for parent in (pose_file.parent, pose_file.parent.parent):
            if (parent / "run_meta.json").is_file():
                return parent
    return None


def _inputs_without_meta(state: WorkflowState) -> PoseInputs | None:
    """Входы, собранные из состояния, когда run_meta.json рядом нет."""
    pose_file = _as_path(state.observations.get("pose_file") or state.request.pose_file)
    if pose_file is None:
        return None

    receptor = _as_path(
        state.observations.get("receptor_pdb_path")
        or state.artifacts.get("receptor_pdb")
    )
    notes: list[str] = []
    pdb_id = str(state.observations.get("pdb_id") or state.request.pdb_id or "").upper()
    if receptor is None and pdb_id:
        for folder in ("data/cache/pdb", "data/external/rcsb"):
            candidate = _REPO_ROOT / folder / f"{pdb_id}.pdb"
            if candidate.is_file():
                receptor = candidate
                notes.append(
                    f"Очищенный рецептор прогона не найден, взят {folder}/{pdb_id}.pdb: "
                    "в нём остаются вода и кристаллический лиганд."
                )
                break

    smiles = state.observations.get("smiles") or state.observations.get("canonical_smiles")
    comp_ids = state.observations.get("ligand_comp_ids") or []
    return PoseInputs(
        pose_file=pose_file,
        receptor_file=receptor,
        smiles=str(smiles) if smiles else None,
        crystal_pdb=_as_path(_REPO_ROOT / "data" / "cache" / "pdb" / f"{pdb_id}.pdb") if pdb_id else None,
        crystal_comp_id=str(comp_ids[0]) if comp_ids else None,
        pdb_id=pdb_id or None,
        ligand_id=state.observations.get("ligand_id") or state.observations.get("ligand_query"),
        notes=notes,
    )


def validate_pose_geometry(state: WorkflowState) -> WorkflowState:
    """Проверить геометрию позы №1 и записать вердикт с измерениями."""
    with state.agent_step("validation_agent", "pose_geometry"):
        run_dir = _run_dir(state)
        if run_dir is not None:
            inputs = load_run_inputs(run_dir)
        else:
            inputs = _inputs_without_meta(state)
        if inputs is None:
            state.add_error(
                "Нет позы для проверки геометрии. Укажите pose_file или выполните докинг."
            )
            return state

        report = pose_validation.validate_pose(
            inputs.pose_file,
            receptor_file=inputs.receptor_file,
            smiles=inputs.smiles,
            crystal_pdb=inputs.crystal_pdb,
            crystal_comp_id=inputs.crystal_comp_id,
            extra_notes=inputs.notes,
        )

        payload = report.as_dict()
        state.observations["geometry_report"] = payload
        state.observations["validation_checks"] = [
            {"name": check.name, "passed": check.passed, "detail": check.summary}
            for check in report.checks
            if check.applicable
        ]
        state.observations["geometry_verdict"] = report.verdict

        context = {
            "pdb_id": inputs.pdb_id,
            "ligand_id": inputs.ligand_id,
            "pocket_id": inputs.pocket_id,
            "pose_id": inputs.pose_id,
            "pose_file": str(inputs.pose_file),
        }

        # Вердикт и счёт — первое, что должно попасть в пересказ.
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=report.verdict_text,
                method="geometry_verdict",
                source="pose_geometry",
                context={**context, "score": report.score_text},
            )
        )
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=len(report.passed_checks),
                unit="проверок",
                method="geometry_checks_passed",
                source="pose_geometry",
                context={**context, "applicable": len(report.applicable)},
            )
        )
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=len(report.applicable),
                unit="проверок",
                method="geometry_checks_applicable",
                source="pose_geometry",
                context=context,
            )
        )
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=report.percent_passed,
                unit="%",
                method="checks_passed_percent",
                source="pose_geometry",
                context={**context, "score": report.score_text},
            )
        )

        # Каждое измерение отдельным фактом: кейс требует объяснимого вердикта,
        # а объяснить его можно только числом и парой атомов, а не статусом.
        for check in report.checks:
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=check.status,
                    method=f"geometry_check_{check.name}",
                    source="pose_geometry",
                    context={**context, "title": check.title, "summary": check.summary},
                )
            )
            shown = check.failures + [item for item in check.measurements if item.passed]
            for item in shown[:_MAX_MEASUREMENTS_PER_CHECK]:
                if item.value is None or isinstance(item.value, str):
                    continue
                state.add_evidence(
                    Evidence(
                        kind="calculation",
                        value=item.value,
                        unit=item.unit or None,
                        method=f"geometry_measure_{check.name}",
                        source="pose_geometry",
                        context={
                            **context,
                            "label": item.label,
                            "passed": item.passed,
                            "allowed": item.range_text,
                        },
                    )
                )

        for violation in report.critical_violations:
            state.limitations.append(f"Критическое нарушение геометрии. {violation}")
        for note in report.notes:
            state.limitations.append(note)
        for check in report.checks:
            if check.status == FAILED and not check.critical:
                state.limitations.append(f"Не пройдена проверка «{check.title}»: {check.summary}")

        report_path = _write_report(report, inputs)
        if report_path is not None:
            state.artifacts["geometry_report_md"] = str(report_path)
            state.artifacts["geometry_report_json"] = str(
                report_path.with_suffix(".json")
            )
    return state


def _write_report(report: GeometryReport, inputs: PoseInputs) -> Path | None:
    """Отчёт рядом с позой: тем же прогоном, который её и посчитал."""
    target_dir = inputs.pose_file.parent
    if target_dir.name == "poses":
        target_dir = target_dir.parent
    try:
        target_dir.mkdir(parents=True, exist_ok=True)
        markdown = target_dir / "geometry_report.md"
        markdown.write_text(render_markdown(report, inputs), encoding="utf-8")
        (target_dir / "geometry_report.json").write_text(
            json.dumps(report.as_dict(), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return markdown
    except OSError:
        return None


#: Прежнее имя шага. Маршрут сценария зовёт `pose_checks` именно через него.
validate_pose = validate_pose_geometry
