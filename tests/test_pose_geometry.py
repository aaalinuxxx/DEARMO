"""Проверки геометрии позы: каждая на данных, а не на заглушке.

Поза для проверок берётся из настоящего редокинга в outputs/docking, а не из
data/examples: файлы там — заглушки по 5-6 атомов с координатами на единичной
сетке у начала координат, и проверка геометрии на них измеряет пустоту.
Искажения (столкновение, растянутая связь, смятое кольцо) строятся прямо в
тесте, чтобы каждая проверка была показана и сработавшей, и несработавшей.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tools import pose_geometry as pg
from tools.contacts import Atom, heavy_atoms, read_structure
from tools.pose_validation import load_run_inputs

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "outputs" / "docking"

#: Прогоны из плана защиты: трипсин, карбоангидраза, CDK2.
BENZAMIDINE_RUN = RUNS / "3PTB_benzamidine_redock"
ACETAZOLAMIDE_RUN = RUNS / "3HS4_acetazolamide_redock"
ROSCOVITINE_RUN = RUNS / "2A4L_r_roscovitine_redock"

rdkit_only = pytest.mark.skipif(
    not pg._RDKIT_AVAILABLE, reason="проверки по молекуле требуют RDKit"
)


def _require(run: Path) -> Path:
    if not (run / "run_meta.json").is_file():
        pytest.skip(f"нет прогона {run.name}: сначала выполните докинг")
    return run


def _atom(serial: int, element: str, x: float, y: float, z: float, **kwargs) -> Atom:
    return Atom(
        serial=serial,
        name=kwargs.get("name", f"{element}{serial}"),
        element=element,
        res_name=kwargs.get("res_name", "LIG"),
        res_seq=kwargs.get("res_seq", 1),
        chain=kwargs.get("chain", "A"),
        x=x,
        y=y,
        z=z,
        is_hetatm=kwargs.get("is_hetatm", False),
    )


# --------------------------------------------------------------------------
# Пороги и примитивы
# --------------------------------------------------------------------------


def test_vdw_radii_come_from_a_table_not_a_default() -> None:
    """Радиусы различают элементы: иначе порог столкновения был бы одинаков."""
    carbon = pg.vdw_radius("C")
    nitrogen = pg.vdw_radius("N")
    oxygen = pg.vdw_radius("O")
    assert carbon > nitrogen > oxygen
    assert 1.0 < oxygen < 2.5


def test_clash_threshold_is_the_documented_fraction_of_the_radius_sum() -> None:
    """Порог кейса — 0,75 суммы радиусов ван-дер-Ваальса пары атомов."""
    expected = 0.75 * (pg.vdw_radius("C") + pg.vdw_radius("O"))
    assert pg.clash_threshold("C", "O", pg.PROTEIN_LIGAND_CLASH_FRACTION) == pytest.approx(expected)
    assert pg.PROTEIN_LIGAND_CLASH_FRACTION == 0.75


def test_plane_deviation_is_zero_for_a_flat_ring_and_grows_when_bent() -> None:
    flat = [(math.cos(a), math.sin(a), 0.0) for a in (0, 1.05, 2.09, 3.14, 4.19, 5.24)]
    assert pg._plane_deviation(flat) == pytest.approx(0.0, abs=1e-6)
    bent = list(flat)
    bent[0] = (bent[0][0], bent[0][1], 0.8)
    assert pg._plane_deviation(bent) > 0.25


def test_nearby_protein_atoms_keeps_the_pocket_and_drops_the_rest() -> None:
    ligand = [_atom(1, "C", 0.0, 0.0, 0.0)]
    near = _atom(2, "O", 3.0, 0.0, 0.0, res_name="ASP")
    far = _atom(3, "O", 99.0, 0.0, 0.0, res_name="ASP")
    selected = pg.nearby_protein_atoms(ligand, [near, far])
    assert near in selected and far not in selected


# --------------------------------------------------------------------------
# Столкновения лиганда с белком
# --------------------------------------------------------------------------


def test_clash_check_names_the_atom_pair_residue_distance_and_threshold() -> None:
    """Кейс требует пару атомов, остаток, расстояние и порог, а не только статус."""
    ligand = [_atom(1, "N", 0.0, 0.0, 0.0)]
    # Кислород остова Gly219 вплотную: заведомо ближе 0,75 суммы радиусов.
    protein = [_atom(2, "O", 1.2, 0.0, 0.0, name="O", res_name="GLY", res_seq=219)]
    check = pg.check_protein_ligand_clashes(ligand, protein)

    assert check.status == pg.FAILED
    assert check.critical, "столкновение с белком обязано быть критическим"
    worst = check.failures[0]
    assert "Gly219" in worst.label and "O" in worst.label
    assert worst.value == pytest.approx(1.2, abs=1e-3)
    assert worst.low == pytest.approx(pg.clash_threshold("N", "O", 0.75), abs=1e-3)


def test_clash_check_ignores_metal_coordination() -> None:
    """Цинк в 2,0 Å от азота — координация, а не столкновение.

    Без этого исключения правильно скоординированный сульфонамид всегда
    объявлялся бы столкновением: координационная связь короче любого порога
    по ван-дер-Ваальсу.
    """
    ligand = [_atom(1, "N", 0.0, 0.0, 0.0)]
    zinc = _atom(2, "ZN", 2.0, 0.0, 0.0, name="ZN", res_name="ZN", is_hetatm=True)
    assert zinc.is_metal
    assert pg.check_protein_ligand_clashes(ligand, [zinc]).status == pg.NOT_APPLICABLE


def test_clash_check_reports_the_closest_contact_even_when_it_passes() -> None:
    """У пройденной проверки тоже есть число: иначе её не с чем сопоставить."""
    ligand = [_atom(1, "C", 0.0, 0.0, 0.0)]
    protein = [_atom(2, "O", 4.0, 0.0, 0.0, name="OD1", res_name="ASP", res_seq=189)]
    check = pg.check_protein_ligand_clashes(ligand, protein)
    assert check.status == pg.PASSED
    assert check.measurements and check.measurements[0].value == pytest.approx(4.0)
    assert "Asp189" in check.measurements[0].label


def test_real_docked_pose_has_no_clashes_with_its_receptor() -> None:
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    ligand = read_structure(inputs.pose_file)
    protein = read_structure(inputs.receptor_file)
    check = pg.check_protein_ligand_clashes(ligand, protein)
    assert check.status == pg.PASSED, check.summary


# --------------------------------------------------------------------------
# Молекула позы, связи, углы, кольца
# --------------------------------------------------------------------------


@rdkit_only
def test_pose_molecule_takes_bond_orders_from_the_run_smiles() -> None:
    """Из PDBQT кратности связей не восстановить — шаблон обязателен."""
    run = _require(ROSCOVITINE_RUN)
    inputs = load_run_inputs(run)
    atoms = read_structure(inputs.pose_file)
    mol, note = pg.build_pose_molecule(atoms, inputs.smiles)

    assert note is None, note
    assert mol.GetNumAtoms() == len(heavy_atoms(atoms))
    assert any(atom.GetIsAromatic() for atom in mol.GetAtoms())


@rdkit_only
def test_pose_molecule_without_smiles_says_so_instead_of_guessing() -> None:
    run = _require(BENZAMIDINE_RUN)
    atoms = read_structure(load_run_inputs(run).pose_file)
    mol, note = pg.build_pose_molecule(atoms, None)
    assert mol is not None
    assert note and "SMILES" in note


@rdkit_only
def test_bond_lengths_and_angles_pass_on_a_real_pose() -> None:
    run = _require(ACETAZOLAMIDE_RUN)
    inputs = load_run_inputs(run)
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)
    lengths = pg.check_bond_lengths(mol)
    angles = pg.check_bond_angles(mol)

    assert lengths.status == pg.PASSED, lengths.summary
    assert angles.status == pg.PASSED, angles.summary
    # Отчёт обязан показывать измеренное значение и допустимый диапазон.
    first = lengths.measurements[0]
    assert first.low is not None and first.high is not None
    assert first.low < first.value < first.high


@rdkit_only
def test_stretched_bond_is_caught_with_its_measured_value_and_range() -> None:
    """Связь, растянутую вдвое, проверка обязана назвать и показать диапазон."""
    from rdkit.Chem import rdMolTransforms

    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)
    assert pg.check_bond_lengths(mol).status == pg.PASSED

    conf = mol.GetConformer()
    bond = mol.GetBondWithIdx(0)
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    original = rdMolTransforms.GetBondLength(conf, i, j)
    position = conf.GetAtomPosition(j)
    position.x += original  # атом уезжает примерно на длину связи
    conf.SetAtomPosition(j, position)

    check = pg.check_bond_lengths(mol)
    assert check.status == pg.FAILED
    broken = check.failures[0]
    assert broken.value > original
    assert broken.high is not None and broken.value > broken.high


@rdkit_only
def test_aromatic_ring_planarity_is_measured_and_a_bent_ring_fails() -> None:
    run = _require(ROSCOVITINE_RUN)
    inputs = load_run_inputs(run)
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)

    check = pg.check_ring_planarity(mol)
    assert check.status == pg.PASSED, check.summary
    assert all(item.value <= pg.RING_PLANARITY_TOLERANCE_A for item in check.measurements)

    ring = next(
        r for r in mol.GetRingInfo().AtomRings()
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
    )
    conf = mol.GetConformer()
    position = conf.GetAtomPosition(ring[0])
    position.z += 1.0  # выгибаем один атом кольца из плоскости
    conf.SetAtomPosition(ring[0], position)

    bent = pg.check_ring_planarity(mol)
    assert bent.status == pg.FAILED
    assert bent.failures[0].value > pg.RING_PLANARITY_TOLERANCE_A


@rdkit_only
def test_internal_clash_only_looks_at_atoms_far_apart_in_the_bond_graph() -> None:
    """Соседи через связь сближены валентностью, и порогу по радиусам не подчиняются."""
    run = _require(ROSCOVITINE_RUN)
    inputs = load_run_inputs(run)
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)
    check = pg.check_internal_clashes(mol)
    assert check.status == pg.PASSED, check.summary
    assert pg.INTERNAL_CLASH_MIN_BOND_SEPARATION >= 4


# --------------------------------------------------------------------------
# Стереохимия
# --------------------------------------------------------------------------


@rdkit_only
def test_defined_stereocentre_is_checked_and_preserved() -> None:
    """У (R)-росковитина центр задан, значит проверка применима и проходит."""
    run = _require(ROSCOVITINE_RUN)
    inputs = load_run_inputs(run)
    assert "@" in (inputs.smiles or ""), "в SMILES прогона должна быть стереохимия"
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)

    check = pg.check_stereochemistry(mol, inputs.smiles)
    assert check.status == pg.PASSED, check.summary
    assert check.measurements[0].value in ("R", "S")


@rdkit_only
def test_molecule_without_stereochemistry_is_not_applicable_not_failed() -> None:
    """«Не применимо» не должно занижать процент пройденности."""
    run = _require(BENZAMIDINE_RUN)
    inputs = load_run_inputs(run)
    mol, _ = pg.build_pose_molecule(read_structure(inputs.pose_file), inputs.smiles)
    check = pg.check_stereochemistry(mol, inputs.smiles)
    assert check.status == pg.NOT_APPLICABLE
    assert not check.applicable


@rdkit_only
def test_inverted_stereocentre_is_reported_as_changed() -> None:
    """Зеркальная поза обязана провалить проверку заданной конфигурации."""
    run = _require(ROSCOVITINE_RUN)
    inputs = load_run_inputs(run)
    atoms = read_structure(inputs.pose_file)
    mirrored = [
        Atom(
            serial=a.serial, name=a.name, element=a.element, res_name=a.res_name,
            res_seq=a.res_seq, chain=a.chain, x=-a.x, y=a.y, z=a.z, is_hetatm=a.is_hetatm,
        )
        for a in atoms
    ]
    mol, _ = pg.build_pose_molecule(mirrored, inputs.smiles)
    check = pg.check_stereochemistry(mol, inputs.smiles)
    assert check.status == pg.FAILED
    assert check.critical
    assert check.failures[0].value in ("R", "S", "не определена")


# --------------------------------------------------------------------------
# Координация металла
# --------------------------------------------------------------------------


def test_metal_check_is_not_applicable_without_a_metal() -> None:
    ligand = [_atom(1, "N", 0.0, 0.0, 0.0)]
    protein = [_atom(2, "O", 2.0, 0.0, 0.0, res_name="ASP")]
    assert pg.check_metal_coordination(ligand, protein).status == pg.NOT_APPLICABLE


def test_metal_check_measures_the_distance_and_the_deviation_from_the_crystal() -> None:
    """Эталон измеряется на кристалле, а не берётся числом из кода."""
    zinc = _atom(9, "ZN", 0.0, 0.0, 0.0, name="ZN", res_name="ZN", is_hetatm=True)
    pose = [_atom(1, "N", 2.6, 0.0, 0.0)]
    crystal = [_atom(1, "N", 1.95, 0.0, 0.0, res_name="AZM", is_hetatm=True)]

    check = pg.check_metal_coordination(pose, [zinc], reference_ligand=crystal)
    assert check.status == pg.FAILED, check.summary
    distance, deviation = check.measurements[0], check.measurements[1]
    assert distance.value == pytest.approx(2.6, abs=1e-3)
    assert deviation.value == pytest.approx(0.65, abs=1e-3)


def test_metal_check_passes_when_the_contact_matches_the_crystal() -> None:
    zinc = _atom(9, "ZN", 0.0, 0.0, 0.0, name="ZN", res_name="ZN", is_hetatm=True)
    pose = [_atom(1, "N", 2.05, 0.0, 0.0)]
    crystal = [_atom(1, "N", 1.95, 0.0, 0.0, res_name="AZM", is_hetatm=True)]
    assert pg.check_metal_coordination(pose, [zinc], reference_ligand=crystal).status == pg.PASSED


def test_metal_check_picks_the_crystal_copy_nearest_the_metal() -> None:
    """В 3HS4 три копии ацетазоламида, эталон — та, что в активном центре."""
    zinc = _atom(9, "ZN", 0.0, 0.0, 0.0, name="ZN", res_name="ZN", is_hetatm=True)
    pose = [_atom(1, "N", 2.0, 0.0, 0.0)]
    crystal = [
        _atom(1, "N", 30.0, 0.0, 0.0, res_name="AZM", res_seq=702, is_hetatm=True),
        _atom(2, "N", 1.95, 0.0, 0.0, res_name="AZM", res_seq=701, is_hetatm=True),
    ]
    check = pg.check_metal_coordination(pose, [zinc], reference_ligand=crystal)
    assert check.status == pg.PASSED
    assert "1.95" in check.summary


def test_ligand_too_far_from_the_metal_is_not_applicable_but_keeps_the_number() -> None:
    """Отсутствие координации — не провал геометрии, но число обязано остаться."""
    zinc = _atom(9, "ZN", 0.0, 0.0, 0.0, name="ZN", res_name="ZN", is_hetatm=True)
    pose = [_atom(1, "N", 8.0, 0.0, 0.0)]
    check = pg.check_metal_coordination(pose, [zinc])
    assert check.status == pg.NOT_APPLICABLE
    assert check.measurements[0].value == pytest.approx(8.0)
