"""RMSD между позами: симметрия, отсутствие выравнивания, отказ вместо догадки."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tools.rmsd import (
    Atom,
    RmsdError,
    crystal_ligand_atoms,
    min_distance_to_atoms,
    parse_ligand_atoms,
    pose_rmsd,
    symmetry_rmsd,
)

ROOT = Path(__file__).resolve().parent.parent

# Бензамидин: кольцо плюс амидиновая группа с двумя эквивалентными азотами.
BENZAMIDINE = """\
HETATM    1  C1  BEN A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  C2  BEN A   1       1.390   0.000   0.000  1.00  0.00           C
HETATM    3  C3  BEN A   1       2.085   1.204   0.000  1.00  0.00           C
HETATM    4  C4  BEN A   1       1.390   2.408   0.000  1.00  0.00           C
HETATM    5  C5  BEN A   1       0.000   2.408   0.000  1.00  0.00           C
HETATM    6  C6  BEN A   1      -0.695   1.204   0.000  1.00  0.00           C
HETATM    7  C7  BEN A   1      -0.730  -1.300   0.000  1.00  0.00           C
HETATM    8  N1  BEN A   1      -0.080  -2.450   0.000  1.00  0.00           N
HETATM    9  N2  BEN A   1      -2.060  -1.300   0.000  1.00  0.00           N
END
"""


def _shift(text: str, dx: float) -> str:
    out = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38]) + dx
            line = f"{line[:30]}{x:8.3f}{line[38:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def test_identical_poses_have_zero_rmsd():
    atoms = parse_ligand_atoms(BENZAMIDINE)
    assert symmetry_rmsd(atoms, atoms) == pytest.approx(0.0, abs=1e-9)


def test_rigid_shift_gives_exactly_that_shift():
    """Структуры не выравниваются: сдвиг на 1,5 Å обязан дать RMSD 1,5 Å."""
    probe = parse_ligand_atoms(_shift(BENZAMIDINE, 1.5))
    reference = parse_ligand_atoms(BENZAMIDINE)
    assert symmetry_rmsd(probe, reference) == pytest.approx(1.5, abs=1e-6)


def test_swapping_equivalent_nitrogens_does_not_inflate_rmsd():
    """N1 и N2 амидиновой группы неразличимы: перестановка не должна наказываться.

    Без учёта симметрии правильная поза получила бы здесь RMSD около 0,4 Å
    просто из-за порядка атомов в файле.
    """
    lines = BENZAMIDINE.splitlines()
    swapped = lines[:7] + [lines[8], lines[7]] + lines[9:]
    probe = parse_ligand_atoms("\n".join(swapped) + "\n")
    reference = parse_ligand_atoms(BENZAMIDINE)
    assert len(probe) == len(reference) == 9
    assert symmetry_rmsd(probe, reference) == pytest.approx(0.0, abs=1e-6)


def test_different_atom_counts_raise_instead_of_returning_a_number():
    """Молчаливое число здесь было бы хуже ошибки: его приняли бы за измерение."""
    reference = parse_ligand_atoms(BENZAMIDINE)
    with pytest.raises(RmsdError, match="Разное число тяжёлых атомов"):
        symmetry_rmsd(reference[:5], reference)


def test_different_composition_raises():
    reference = parse_ligand_atoms(BENZAMIDINE)
    probe = list(reference)
    probe[0] = Atom("O", probe[0].x, probe[0].y, probe[0].z)
    with pytest.raises(RmsdError, match="Разный состав"):
        symmetry_rmsd(probe, reference)


def test_hydrogens_are_excluded_from_heavy_atom_rmsd():
    with_h = BENZAMIDINE.replace(
        "END",
        "HETATM   10  H1  BEN A   1       5.000   5.000   5.000  1.00  0.00           H\nEND",
    )
    assert len(parse_ligand_atoms(with_h)) == 9


def test_autodock_polar_hydrogen_type_is_not_mistaken_for_an_element():
    """Тип HD в PDBQT — водород, а не что-то тяжёлое."""
    line = (
        "ATOM      1  HD1 UNL     1       1.000   2.000   3.000  0.00  0.00    +0.000 HD\n"
    )
    assert parse_ligand_atoms(line) == []


def test_autodock_types_map_back_to_real_elements():
    """OA, NA и SA — кислород, азот и сера, а не натрий."""
    text = (
        "ATOM      1  O   UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 OA\n"
        "ATOM      2  N   UNL     1       1.000   0.000   0.000  0.00  0.00    +0.000 NA\n"
        "ATOM      3  S   UNL     1       2.000   0.000   0.000  0.00  0.00    +0.000 SA\n"
    )
    assert [a.element for a in parse_ligand_atoms(text)] == ["O", "N", "S"]


def test_aromatic_carbon_type_a_is_carbon_not_an_unknown_element():
    """Meeko пишет ароматический углерод типом A; принять его за элемент нельзя."""
    line = "ATOM      1  C1  UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000  A\n"
    assert [a.element for a in parse_ligand_atoms(line)] == ["C"]


def test_only_first_model_of_a_multi_model_file_is_read():
    text = (
        "MODEL 1\n"
        "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00     C\n"
        "ENDMDL\n"
        "MODEL 2\n"
        "ATOM      1  C   UNL     1       9.000   9.000   9.000  0.00  0.00     C\n"
        "ENDMDL\n"
    )
    atoms = parse_ligand_atoms(text)
    assert len(atoms) == 1 and atoms[0].x == 0.0


def test_min_distance_measures_to_nearest_atom_not_to_centroid():
    atoms = [Atom("C", 0.0, 0.0, 0.0), Atom("C", 10.0, 0.0, 0.0)]
    assert min_distance_to_atoms((1.0, 0.0, 0.0), atoms) == pytest.approx(1.0)


def test_crystal_ligand_is_selected_by_comp_id():
    text = BENZAMIDINE.replace("END", "HETATM   10  O   HOH A   9       9.0   9.0   9.0  1.00  0.00           O\nEND")
    assert len(crystal_ligand_atoms(text, "BEN")) == 9
    with pytest.raises(RmsdError, match="нет лиганда"):
        crystal_ligand_atoms(text, "AZM")


def test_pose_rmsd_reads_files(tmp_path: Path):
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    a.write_text(BENZAMIDINE, encoding="utf-8")
    b.write_text(_shift(BENZAMIDINE, 2.0), encoding="utf-8")
    assert pose_rmsd(a, b) == pytest.approx(2.0, abs=1e-6)


def test_placeholder_demo_pose_is_not_benzamidine():
    """Страховка от возврата заглушки, выдаваемой за контрольный прогон.

    В data/examples/3ptb_ben лежала поза из пяти атомов в начале координат с
    записанным рядом RMSD 0,356 Å. Ни того ни другого никто не вычислял. Если
    такой файл снова начнут выдавать за бензамидин, тест это покажет.
    """
    pose = ROOT / "data" / "examples" / "3ptb_ben" / "poses" / "BEN_pose1.pdbqt"
    if not pose.is_file():
        pytest.skip("демо-файла нет")
    atoms = parse_ligand_atoms(pose.read_text(encoding="utf-8"))
    crystal = crystal_ligand_atoms(
        (ROOT / "data" / "external" / "rcsb" / "3PTB.pdb").read_text(encoding="utf-8", errors="replace"),
        "BEN",
    )
    assert len(crystal) == 9, "у бензамидина девять тяжёлых атомов"
    if len(atoms) == len(crystal):
        pytest.skip("демо-поза заменена настоящей")
    with pytest.raises(RmsdError):
        symmetry_rmsd(atoms, crystal)


@pytest.mark.integration
def test_real_redock_of_3ptb_reproduces_the_crystal_pose():
    """Сквозная проверка: докинг вслепую обязан воспроизвести кристалл.

    Требует Vina и P2Rank. Это и есть контрольный результат для защиты —
    в отличие от прежнего, он считается, а не записан в файл руками.
    """
    from tools.p2rank_runner import p2rank_binary
    from tools.vina_runner import vina_binary

    if vina_binary() is None or p2rank_binary() is None:
        pytest.skip("нужны Vina и P2Rank")

    import tempfile

    from tools.live_docking import dock_pockets

    with tempfile.TemporaryDirectory() as tmp:
        run = dock_pockets(
            "3PTB",
            "NC(=N)c1ccccc1",
            "benzamidine",
            Path(tmp),
            chain="A",
            crystal_comp_id="BEN",
        )
    assert run.box_source == "p2rank", "карманы обязаны прийти от P2Rank, а не от центра белка"
    assert len(run.pockets) == 3
    assert len({p.pocket_id for p in run.poses}) >= 2, "докинг идёт в каждый карман, не в один"
    control = run.control
    assert control["pocket_hit_rank"] == 1, "истинный карман трипсина должен быть первым"
    assert control["best_pose_rmsd_to_crystal_a"] < 2.0
    assert control["demo_check_passed"] is True
    # Первая поза каждого кармана — точка отсчёта для RMSD внутри кармана.
    for pose in run.poses:
        if pose.pose_id == 1:
            assert pose.rmsd_lb == 0.0 and pose.rmsd_ub == 0.0
        else:
            assert pose.rmsd_lb is not None and pose.rmsd_lb > 0
    assert math.isfinite(run.best.vina_score) and run.best.vina_score < 0


def test_pdbqt_halogen_types_are_not_truncated_to_carbon():
    """Тип «Cl» в колонках 78–79 PDBQT читался как «C»: RMSD гефитиниба не считался."""
    from tools.rmsd import _element_from_line

    def pdbqt(atom_type: str) -> str:
        return f"ATOM      1  X   UNL     1       0.000   0.000   0.000  1.00  0.00    -0.081 {atom_type:<2}"

    assert _element_from_line(pdbqt("Cl")) == "CL"
    assert _element_from_line(pdbqt("Br")) == "BR"
    assert _element_from_line(pdbqt("OA")) == "O"
    assert _element_from_line(pdbqt("A")) == "C"
