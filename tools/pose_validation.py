"""Сценарий 4: проверка геометрии позы №1 и вердикт по ней.

Модуль собирает шесть проверок из `tools.pose_geometry` в один отчёт и
отвечает на вопрос кейса: прошла ли поза и какие именно критерии.

Два правила, из-за которых отчёт устроен именно так.

Первое: вердикт детерминирован. Его считает код по измерениям, а не языковая
модель по тексту. Модель в этом проекте ошибается не числами, а привязкой
чисел, и «проверка пройдена» при проваленной проверке не поймала бы ни одна
сверка значений. Модели остаётся пересказ уже готового вердикта.

Второе: критическое нарушение не перевешивается арифметикой. Результат «4 из
5» при столкновении лиганда с белком — это всё равно отклонённая поза, и
итоговый статус обязан это показывать.

Поза берётся та, что посчитана в этом же прогоне: сценарий 2 записывает её в
run_meta.json → best_pose как глобальный минимум Vina по всем карманам при
одинаковых настройках скоринга. Здесь она не выбирается заново.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from tools.contacts import Atom, read_structure
from tools.pose_geometry import (
    FAILED,
    NOT_APPLICABLE,
    PASSED,
    GeometryCheck,
    GeometryReport,
    build_pose_molecule,
    check_bond_angles,
    check_bond_lengths,
    check_internal_clashes,
    check_metal_coordination,
    check_protein_ligand_clashes,
    check_ring_planarity,
    check_stereochemistry,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent

#: Порядок проверок в отчёте — порядок таблицы сценария 4 в описании кейса.
CHECK_ORDER = (
    "protein_ligand_clash",
    "bond_lengths",
    "bond_angles",
    "stereochemistry",
    "ring_planarity",
    "internal_clash",
    "metal_coordination",
)


@dataclass
class PoseInputs:
    """Всё, что нужно для проверки позы, с указанием, откуда это взято."""

    pose_file: Path
    receptor_file: Path | None = None
    smiles: str | None = None
    crystal_pdb: Path | None = None
    crystal_comp_id: str | None = None
    pdb_id: str | None = None
    ligand_id: str | None = None
    pocket_id: str | None = None
    pose_id: int | None = None
    vina_score: float | None = None
    notes: list[str] = field(default_factory=list)


def _resolve(path: Any) -> Path | None:
    """Путь из run_meta.json: относительные считаются от корня репозитория."""
    if not path:
        return None
    candidate = Path(str(path))
    if not candidate.is_absolute():
        candidate = _REPO_ROOT / candidate
    return candidate if candidate.is_file() else None


def load_run_inputs(run_dir: Path | str) -> PoseInputs:
    """Входы сценария 4 из каталога прогона сценария 2.

    Поза и рецептор не ищутся перебором каталогов: оба адреса записаны в
    run_meta.json тем же прогоном, который их и сделал. Перебор мог бы
    подставить позу из чужого запуска, а это ровно то, чего кейс запрещает.
    """
    run_dir = Path(run_dir)
    meta_path = run_dir / "run_meta.json"
    if not meta_path.is_file():
        raise FileNotFoundError(f"Нет run_meta.json в {run_dir}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))

    notes: list[str] = []
    best = meta.get("best_pose") or {}
    pose_file = _resolve(best.get("pose_file"))
    if pose_file is None:
        raise FileNotFoundError(f"В {meta_path} не указан файл позы №1 (best_pose.pose_file)")

    receptor = _resolve(meta.get("receptor_pdb_path"))
    if receptor is None:
        # Прогоны, записанные до появления receptor_pdb_path, обходятся
        # соседним файлом или кешем PDB. Факт замены попадает в отчёт:
        # молча подставлять другой белок нельзя.
        for fallback, label in (
            (run_dir / "receptor.pdb", "receptor.pdb рядом с позой"),
            (_REPO_ROOT / "data" / "cache" / "pdb" / f"{str(meta.get('pdb_id') or '').upper()}.pdb", "кеш PDB"),
            (_REPO_ROOT / "data" / "external" / "rcsb" / f"{str(meta.get('pdb_id') or '').upper()}.pdb", "загруженный PDB"),
        ):
            if fallback.is_file():
                receptor = fallback
                notes.append(f"receptor_pdb_path в run_meta.json не указан, взят {label}.")
                break

    control = meta.get("control") or {}
    return PoseInputs(
        pose_file=pose_file,
        receptor_file=receptor,
        smiles=meta.get("smiles"),
        crystal_pdb=_resolve(meta.get("pdb_path")),
        crystal_comp_id=control.get("crystal_comp_id"),
        pdb_id=meta.get("pdb_id"),
        ligand_id=meta.get("ligand_id"),
        pocket_id=best.get("pocket_id"),
        pose_id=best.get("pose_id"),
        vina_score=best.get("vina_score"),
        notes=notes,
    )


def crystal_ligand_atoms(pdb_file: Path | str, comp_id: str) -> list[Atom]:
    """Атомы кристаллического лиганда по коду компонента PDB.

    Копии не фильтруются: в 3HS4 ацетазоламид присутствует трижды, и нужная
    копия — та, что лежит в активном центре у иона цинка. Выбор делает сама
    проверка координации, беря ближайшую к металлу пару, вместо того чтобы
    угадывать копию по номеру остатка.
    """
    wanted = str(comp_id).strip().upper()
    if not wanted:
        return []
    return [
        atom
        for atom in read_structure(pdb_file)
        if atom.res_name.strip().upper() == wanted
    ]


def validate_pose(
    pose_file: Path | str,
    *,
    receptor_file: Path | str | None = None,
    smiles: str | None = None,
    crystal_pdb: Path | str | None = None,
    crystal_comp_id: str | None = None,
    extra_notes: list[str] | None = None,
) -> GeometryReport:
    """Шесть геометрических проверок позы и общий вердикт по ним."""
    report = GeometryReport(notes=list(extra_notes or []))

    pose_path = Path(pose_file)
    if not pose_path.is_file():
        report.checks.append(
            GeometryCheck(
                name="pose_file",
                title="Файл позы",
                status=FAILED,
                summary=f"Файл позы не найден: {pose_path}",
                critical=True,
            )
        )
        return report

    ligand = read_structure(pose_path)
    if not ligand:
        report.checks.append(
            GeometryCheck(
                name="pose_file",
                title="Файл позы",
                status=FAILED,
                summary=f"В файле позы нет координат атомов: {pose_path}",
                critical=True,
            )
        )
        return report

    protein: list[Atom] = []
    if receptor_file and Path(receptor_file).is_file():
        protein = read_structure(receptor_file)
    else:
        report.notes.append(
            "Координаты белка не переданы: проверки столкновений с белком и "
            "координации металла не выполнялись."
        )

    molecule, build_note = build_pose_molecule(ligand, smiles)
    if build_note:
        report.notes.append(build_note)

    reference: list[Atom] = []
    if crystal_pdb and crystal_comp_id and Path(crystal_pdb).is_file():
        reference = crystal_ligand_atoms(crystal_pdb, crystal_comp_id)
        if not reference:
            report.notes.append(
                f"В {Path(crystal_pdb).name} нет компонента {crystal_comp_id}: "
                "эталон для координации металла не взят."
            )

    checks = [
        check_protein_ligand_clashes(ligand, protein),
        check_bond_lengths(molecule),
        check_bond_angles(molecule),
        check_stereochemistry(molecule, smiles),
        check_ring_planarity(molecule),
        check_internal_clashes(molecule),
        check_metal_coordination(ligand, protein, reference_ligand=reference or None),
    ]
    order = {name: index for index, name in enumerate(CHECK_ORDER)}
    report.checks = sorted(checks, key=lambda check: order.get(check.name, 99))
    return report


def validate_run(run_dir: Path | str) -> tuple[GeometryReport, PoseInputs]:
    """Проверить позу №1 прогона сценария 2 по его run_meta.json."""
    inputs = load_run_inputs(run_dir)
    report = validate_pose(
        inputs.pose_file,
        receptor_file=inputs.receptor_file,
        smiles=inputs.smiles,
        crystal_pdb=inputs.crystal_pdb,
        crystal_comp_id=inputs.crystal_comp_id,
        extra_notes=inputs.notes,
    )
    return report, inputs


# --------------------------------------------------------------------------
# Представление отчёта
# --------------------------------------------------------------------------

_STATUS_TEXT = {
    PASSED: "пройдена",
    FAILED: "не пройдена",
    NOT_APPLICABLE: "не применимо",
}


def status_text(status: str) -> str:
    return _STATUS_TEXT.get(status, status)


def render_markdown(report: GeometryReport, inputs: PoseInputs | None = None) -> str:
    """Отчёт сценария 4 таблицей: что проверяли, что померили, что вышло."""
    lines: list[str] = ["# Проверка геометрии позы №1", ""]
    if inputs is not None:
        head = [
            f"- Структура: {inputs.pdb_id or '—'}",
            f"- Лиганд: {inputs.ligand_id or '—'}",
            f"- Поза: карман {inputs.pocket_id or '—'}, поза {inputs.pose_id or '—'}",
            f"- Оценка Vina: {inputs.vina_score if inputs.vina_score is not None else '—'} ккал/моль",
            f"- Файл позы: `{inputs.pose_file}`",
        ]
        lines.extend(head)
        lines.append("")

    lines.append(f"**Вердикт: {report.verdict_text}.** "
                 f"Пройдено {report.score_text} применимых проверок ({report.percent_passed}%).")
    lines.append("")
    lines.append("| Проверка | Статус | Измерение | Допустимо |")
    lines.append("| --- | --- | --- | --- |")
    for check in report.checks:
        shown = check.failures[:1] or check.measurements[:1]
        if shown:
            item = shown[0]
            measured = f"{item.label} = {item.value_text}"
            allowed = item.range_text
        else:
            measured = check.summary
            allowed = "—"
        lines.append(
            f"| {check.title} | {status_text(check.status)} | {measured} | {allowed} |"
        )

    if report.critical_violations:
        lines.extend(["", "## Критические нарушения", ""])
        lines.extend(f"- {item}" for item in report.critical_violations)

    if report.notes:
        lines.extend(["", "## Примечания", ""])
        lines.extend(f"- {note}" for note in report.notes)
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Совместимость с прежним интерфейсом
# --------------------------------------------------------------------------


@dataclass
class CheckResult:
    """Плоский результат проверки: имя, флаг, пояснение."""

    name: str
    passed: bool
    detail: str


def as_check_results(report: GeometryReport) -> list[CheckResult]:
    """Отчёт в плоский список. Неприменимые проверки не попадают в счёт."""
    return [
        CheckResult(name=check.name, passed=check.passed, detail=check.summary)
        for check in report.checks
        if check.applicable
    ]


def validate_pose_file(
    pose_file: Path | str,
    protein_file: Path | str | None = None,
    *,
    smiles: str | None = None,
) -> list[CheckResult]:
    """Прежняя сигнатура поверх нового отчёта."""
    report = validate_pose(pose_file, receptor_file=protein_file, smiles=smiles)
    results = as_check_results(report)
    if not results:
        return [
            CheckResult(
                name="pose_file",
                passed=False,
                detail="; ".join(check.summary for check in report.checks) or "Проверок не выполнено",
            )
        ]
    return results


def pass_rate(results: list[CheckResult]) -> float:
    """Процент пройденных проверок."""
    if not results:
        return 0.0
    passed = sum(1 for item in results if item.passed)
    return round(100.0 * passed / len(results), 1)
