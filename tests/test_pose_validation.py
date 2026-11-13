"""Сценарий 4 целиком: вердикт, счёт применимых проверок и объяснимость.

План защиты требует проверить сценарий 4 на обычной позе и на её контрольной
копии с намеренно созданным столкновением, убедившись, что система называет
нужную пару атомов и показывает нарушенный порог. Контрольная копия строится
здесь же: сдвиг позы к ближайшему остатку — самый короткий способ получить
столкновение, не трогая ничего другого в геометрии молекулы.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tools import pose_geometry as pg
from tools.contacts import heavy_atoms, read_structure
from tools.pose_validation import (
    CHECK_ORDER,
    load_run_inputs,
    pass_rate,
    render_markdown,
    validate_pose,
    validate_pose_file,
    validate_run,
)

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "outputs" / "docking"

BENZAMIDINE_RUN = RUNS / "3PTB_benzamidine_redock"
ACETAZOLAMIDE_RUN = RUNS / "3HS4_acetazolamide_redock"
ROSCOVITINE_RUN = RUNS / "2A4L_r_roscovitine_redock"


def _require(run: Path) -> Path:
    if not (run / "run_meta.json").is_file():
        pytest.skip(f"нет прогона {run.name}: сначала выполните докинг")
    return run


def make_clashing_copy(run: Path, destination: Path, *, target_fraction: float = 0.6) -> Path:
    """Контрольная копия позы с намеренно созданным столкновением.

    Поза целиком сдвигается вдоль вектора к ближайшему тяжёлому атому белка,
    пока эта пара не окажется заметно ближе порога. Внутренняя геометрия
    молекулы при переносе не меняется — значит провалиться обязана ровно
    одна проверка, столкновение с белком, и по ней видно, что вердикт
    реагирует именно на неё.
    """
    inputs = load_run_inputs(run)
    ligand = heavy_atoms(read_structure(inputs.pose_file))
    protein = [a for a in heavy_atoms(read_structure(inputs.receptor_file)) if not a.is_metal]
    candidates = pg.nearby_protein_atoms(ligand, protein)

    lig_atom, prot_atom, distance = min(
        ((l, p, l.distance_to(p)) for l in ligand for p in candidates),
        key=lambda item: item[2],
    )
    threshold = pg.clash_threshold(lig_atom.element, prot_atom.element, pg.PROTEIN_LIGAND_CLASH_FRACTION)
    scale = (distance - target_fraction * threshold) / distance
    shift = (
        (prot_atom.x - lig_atom.x) * scale,
        (prot_atom.y - lig_atom.y) * scale,
        (prot_atom.z - lig_atom.z) * scale,
    )

    lines = []
    for line in Path(inputs.pose_file).read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38]) + shift[0]
            y = float(line[38:46]) + shift[1]
            z = float(line[46:54]) + shift[2]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        lines.append(line)
    destination.write_text("\n".join(lines), encoding="utf-8")
    return destination


# --------------------------------------------------------------------------
# Обычная поза
# --------------------------------------------------------------------------


def test_real_pose_passes_and_every_check_is_reported() -> None:
    run = _require(BENZAMIDINE_RUN)
    report, inputs = validate_run(run)

    assert report.verdict == pg.VERDICT_PASSED, report.critical_violations
    assert report.percent_passed == 100.0
    # Все семь строк таблицы сценария 4 присутствуют, в том числе неприменимые.
    assert {check.name for check in report.checks} == set(CHECK_ORDER)
    assert inputs.pdb_id == "3PTB"


def test_percent_counts_only_applicable_checks() -> None:
    """Молекула без хиральных центров не штрафуется за их отсутствие."""
    run = _require(BENZAMIDINE_RUN)
    report, _ = validate_run(run)

    not_applicable = [c for c in report.checks if not c.applicable]
    assert not_applicable, "у бензамидина нет ни стереоцентров, ни металла"
    assert len(report.applicable) == len(report.checks) - len(not_applicable)
    assert report.score_text == f"{len(report.passed_checks)} из {len(report.applicable)}"


def test_every_check_keeps_the_numbers_behind_its_status() -> None:
    """Статус без измерения нельзя объяснить, поэтому числа сохраняются всегда."""
    run = _require(ROSCOVITINE_RUN)
    report, _ = validate_run(run)
    for check in report.checks:
        if check.status == pg.PASSED:
            assert check.measurements, f"{check.name} прошла, но ничего не померила"
            assert any(item.value is not None for item in check.measurements)


# --------------------------------------------------------------------------
# Контрольная копия со столкновением
# --------------------------------------------------------------------------


def test_clashing_control_copy_is_rejected_and_names_the_atom_pair(tmp_path: Path) -> None:
    """Главная проверка плана защиты: нужная пара атомов и нарушенный порог."""
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    clashed = make_clashing_copy(run, tmp_path / "clashed.pdbqt")

    report = validate_pose(
        clashed,
        receptor_file=inputs.receptor_file,
        smiles=inputs.smiles,
    )

    clash = report.check("protein_ligand_clash")
    assert clash.status == pg.FAILED
    worst = clash.failures[0]
    assert worst.value < worst.low, "расстояние обязано быть меньше порога"
    # Остаток белка назван, а не заменён на номер атома.
    assert any(ch.isdigit() for ch in worst.label) and "—" in worst.label
    assert report.verdict == pg.VERDICT_REJECTED
    assert report.critical_violations


def test_clash_overrides_the_arithmetic_of_passed_checks(tmp_path: Path) -> None:
    """«4 из 5» при столкновении — всё равно отклонённая поза."""
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    clashed = make_clashing_copy(run, tmp_path / "clashed.pdbqt")
    report = validate_pose(clashed, receptor_file=inputs.receptor_file, smiles=inputs.smiles)

    assert report.percent_passed >= 50.0, "большинство проверок всё ещё проходит"
    assert report.verdict == pg.VERDICT_REJECTED
    assert "отклонена" in report.verdict_text


def test_moving_the_pose_breaks_only_the_clash_check(tmp_path: Path) -> None:
    """Перенос не меняет внутреннюю геометрию: связи и кольца обязаны уцелеть."""
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    clashed = make_clashing_copy(run, tmp_path / "clashed.pdbqt")
    report = validate_pose(clashed, receptor_file=inputs.receptor_file, smiles=inputs.smiles)

    for name in ("bond_lengths", "bond_angles", "ring_planarity", "internal_clash"):
        assert report.check(name).status != pg.FAILED, name


# --------------------------------------------------------------------------
# Карбоангидраза: контакт с цинком
# --------------------------------------------------------------------------


def test_acetazolamide_metal_contact_is_measured_against_the_crystal() -> None:
    """Для 3HS4 кейс требует расстояние N–Zn и отклонение от кристалла."""
    run = _require(ACETAZOLAMIDE_RUN)
    report, inputs = validate_run(run)

    assert inputs.crystal_comp_id == "AZM"
    metal = report.check("metal_coordination")
    assert metal.applicable, metal.summary
    distance = metal.measurements[0]
    assert "ZN" in distance.label.upper()
    assert 1.5 < distance.value < pg.METAL_COORDINATION_MAX_A
    # Эталон измерен на кристалле, а не взят числом из кода.
    deviation = next(item for item in metal.measurements if "кристалл" in item.label.lower())
    assert deviation.value is not None


def test_pose_far_from_the_crystal_contact_does_not_silently_pass() -> None:
    """Поза №1 для 3HS4 воспроизводит кристалл плохо, и отчёт это показывает.

    Vina не моделирует координацию сульфонамида с цинком, поэтому поза уходит
    от кристаллической. Сценарий 4 обязан показать это числом, а не выдать
    «пройдено» по умолчанию.
    """
    run = _require(ACETAZOLAMIDE_RUN)
    report, _ = validate_run(run)
    metal = report.check("metal_coordination")

    if metal.status == pg.FAILED:
        assert metal.failures[0].value is not None
        assert report.verdict != pg.VERDICT_PASSED
    else:
        # Если докинг когда-нибудь воспроизведёт контакт, отклонение всё равно
        # должно быть измерено и показано.
        assert any("кристалл" in item.label.lower() for item in metal.measurements)


# --------------------------------------------------------------------------
# Отчёт и прежний интерфейс
# --------------------------------------------------------------------------


def test_markdown_report_states_verdict_score_and_thresholds() -> None:
    run = _require(ACETAZOLAMIDE_RUN)
    report, inputs = validate_run(run)
    text = render_markdown(report, inputs)

    assert "Вердикт" in text and report.verdict_text in text
    assert report.score_text in text
    assert "3HS4" in text
    assert "| Проверка | Статус | Измерение | Допустимо |" in text


def test_missing_pose_file_is_a_failure_not_a_crash(tmp_path: Path) -> None:
    report = validate_pose(tmp_path / "нет-такого.pdbqt")
    assert report.verdict == pg.VERDICT_REJECTED
    assert report.percent_passed == 0.0


def test_pose_without_receptor_says_so_instead_of_inventing_a_verdict() -> None:
    """Без координат белка проверка столкновений не выполняется, а не «проходит»."""
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    report = validate_pose(inputs.pose_file, smiles=inputs.smiles)

    assert report.check("protein_ligand_clash").status == pg.NOT_APPLICABLE
    assert any("белк" in note.lower() for note in report.notes)


def test_legacy_interface_still_works() -> None:
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    results = validate_pose_file(inputs.pose_file, inputs.receptor_file, smiles=inputs.smiles)
    assert results
    assert pass_rate(results) == 100.0
