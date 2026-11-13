"""Сценарий 3: контакты, сравнение поз и аналитика рядов.

Геометрия проверяется на синтетических структурах с заранее известными
расстояниями: тест должен падать от ошибки в формуле, а не от того, что
в репозитории сменился прогон докинга. Отдельно отмечены тесты, которым
нужны живые результаты сценария 2.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from models.docking import PoseRow
from tools import contacts as contacts_mod
from tools.contacts import (
    Atom,
    count_models,
    detect_contacts,
    heavy_atoms,
    ligand_charged_atoms,
    parse_atoms,
    protein_charge,
    read_structure,
)
from tools.pose_metrics import (
    PoseAnalysis,
    PoseKey,
    contact_matrix,
    key_atom_pairs,
    distance_series,
    pocket_baselines,
    protein_state_rows,
    select_for_comparison,
    series_findings,
    vina_deltas,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DOCKING = _REPO_ROOT / "outputs" / "docking"


# --------------------------------------------------------------------------
# Вспомогательные построители
# --------------------------------------------------------------------------


def _pdb_line(serial, name, res_name, res_seq, x, y, z, element, record="ATOM"):
    return (
        f"{record:<6}{serial:>5} {name:^4} {res_name:>3} A{res_seq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
    )


def _pdbqt_line(serial, name, res_name, res_seq, x, y, z, ad_type):
    return (
        f"ATOM  {serial:>5} {name:^4} {res_name:>3} A{res_seq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +0.000 {ad_type:>2}"
    )


def _row(pocket_id: str, pose_id: int, score: float, pose_file: str = "x.pdbqt") -> PoseRow:
    return PoseRow(
        ligand_id="lig",
        pocket_id=pocket_id,
        pocket_rank=int(pocket_id[-1]) if pocket_id[-1].isdigit() else 1,
        p2rank_score=1.0,
        pose_id=pose_id,
        vina_score=score,
        pose_file=pose_file,
    )


def _analysis(pocket_id, pose_id, score, ligand_atoms=(), contacts=()) -> PoseAnalysis:
    row = _row(pocket_id, pose_id, score)
    return PoseAnalysis(
        key=PoseKey.from_row(row),
        row=row,
        ligand_atoms=list(ligand_atoms),
        contacts=list(contacts),
    )


# --------------------------------------------------------------------------
# Чтение координат
# --------------------------------------------------------------------------


def test_pdb_element_comes_from_element_column():
    atoms = parse_atoms(_pdb_line(1, "OD1", "ASP", 189, 1.0, 2.0, 3.0, "O"), pdbqt=False)
    assert len(atoms) == 1
    atom = atoms[0]
    assert atom.element == "O"
    assert atom.res_name == "ASP"
    assert atom.res_seq == 189
    assert atom.residue_label == "Asp189"
    assert (atom.x, atom.y, atom.z) == (1.0, 2.0, 3.0)


@pytest.mark.parametrize(
    "ad_type,expected",
    [("A", "C"), ("C", "C"), ("NA", "N"), ("OA", "O"), ("SA", "S"), ("HD", "H"), ("Zn", "ZN")],
)
def test_pdbqt_element_restored_from_autodock_type(ad_type, expected):
    """В PDBQT колонка элемента занята типом AutoDock, а не элементом.

    Ароматический углерод Meeko пишет как «A»: без карты типов он ушёл бы
    в неизвестные элементы, а цинк прочитался бы углеродом.
    """
    atoms = parse_atoms(_pdbqt_line(1, "C1", "LIG", 1, 0.0, 0.0, 0.0, ad_type), pdbqt=True)
    assert atoms[0].element == expected


def test_multimodel_pdbqt_selects_one_model():
    """Файл Vina многомодельный: без выбора модели конформеры склеились бы."""
    text = "\n".join(
        [
            "MODEL 1",
            _pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"),
            "ENDMDL",
            "MODEL 2",
            _pdbqt_line(1, "N1", "LIG", 1, 9.0, 0.0, 0.0, "NA"),
            "ENDMDL",
        ]
    )
    assert count_models(text) == 2
    assert parse_atoms(text, pdbqt=True, model=1)[0].x == 0.0
    assert parse_atoms(text, pdbqt=True, model=2)[0].x == 9.0
    assert len(parse_atoms(text, pdbqt=True)) == 2


def test_hydrogens_dropped_from_heavy_atoms():
    atoms = parse_atoms(
        "\n".join(
            [
                _pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"),
                _pdbqt_line(2, "H1", "LIG", 1, 1.0, 0.0, 0.0, "HD"),
            ]
        ),
        pdbqt=True,
    )
    assert len(atoms) == 2
    assert [a.name for a in heavy_atoms(atoms)] == ["N1"]


# --------------------------------------------------------------------------
# Типы взаимодействий на заданной геометрии
# --------------------------------------------------------------------------


def _ligand_amidine(offset: float) -> list[Atom]:
    """Амидин: углерод с двумя азотами. Заряд распознаётся по графу связей."""
    return parse_atoms(
        "\n".join(
            [
                _pdbqt_line(1, "N1", "LIG", 1, offset, 0.0, 0.0, "NA"),
                _pdbqt_line(2, "C1", "LIG", 1, offset + 1.33, 0.0, 0.0, "C"),
                _pdbqt_line(3, "N2", "LIG", 1, offset + 2.0, 1.15, 0.0, "NA"),
                _pdbqt_line(4, "C2", "LIG", 1, offset + 2.0, -1.25, 0.0, "C"),
            ]
        ),
        pdbqt=True,
    )


def test_amidine_nitrogens_detected_as_positive():
    charges = ligand_charged_atoms(_ligand_amidine(0.0))
    assert charges.get(1) == 1
    assert charges.get(3) == 1
    assert 2 not in charges


def test_protein_charges_follow_residue_and_atom_name():
    assert protein_charge(Atom(1, "OD1", "O", "ASP", 189, "A", 0, 0, 0)) == -1
    assert protein_charge(Atom(1, "NZ", "N", "LYS", 10, "A", 0, 0, 0)) == 1
    assert protein_charge(Atom(1, "CA", "C", "ASP", 189, "A", 0, 0, 0)) == 0
    assert protein_charge(Atom(1, "OG", "O", "SER", 190, "A", 0, 0, 0)) == 0


def test_salt_bridge_found_and_distance_matches_coordinates():
    """Амидин против карбоксилата Asp189 — пара из описания кейса."""
    ligand = _ligand_amidine(0.0)
    protein = parse_atoms(
        _pdb_line(100, "OD1", "ASP", 189, -2.8, 0.0, 0.0, "O"), pdbqt=False
    )
    found = detect_contacts(ligand, protein, include_rings=False)
    assert len(found) == 1
    contact = found[0]
    assert contact.interaction_type == "salt_bridge"
    assert contact.residue_label == "Asp189"
    assert contact.distance_a == pytest.approx(2.8, abs=1e-3)


def test_hydrogen_bond_when_residue_is_not_charged():
    ligand = _ligand_amidine(0.0)
    protein = parse_atoms(_pdb_line(100, "OG", "SER", 190, -3.0, 0.0, 0.0, "O"), pdbqt=False)
    found = detect_contacts(ligand, protein, include_rings=False)
    assert [c.interaction_type for c in found] == ["hydrogen_bond"]
    assert found[0].distance_a == pytest.approx(3.0, abs=1e-3)


def test_metal_coordination_detected_for_zinc():
    """Азот лиганда у Zn²⁺ — проверка, которую кейс требует для 3HS4."""
    ligand = parse_atoms(_pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"), pdbqt=True)
    protein = parse_atoms(
        _pdb_line(500, "ZN", "ZN", 301, 2.1, 0.0, 0.0, "ZN", record="HETATM"), pdbqt=False
    )
    found = detect_contacts(ligand, protein, include_rings=False)
    assert [c.interaction_type for c in found] == ["metal"]
    assert found[0].distance_a == pytest.approx(2.1, abs=1e-3)


def test_metal_not_reported_beyond_cutoff():
    ligand = parse_atoms(_pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"), pdbqt=True)
    protein = parse_atoms(
        _pdb_line(500, "ZN", "ZN", 301, 3.4, 0.0, 0.0, "ZN", record="HETATM"), pdbqt=False
    )
    assert detect_contacts(ligand, protein, include_rings=False) == []


def test_hydrophobic_contact_between_carbons():
    ligand = parse_atoms(_pdbqt_line(1, "C1", "LIG", 1, 0.0, 0.0, 0.0, "C"), pdbqt=True)
    protein = parse_atoms(_pdb_line(100, "CB", "VAL", 213, 3.8, 0.0, 0.0, "C"), pdbqt=False)
    found = detect_contacts(ligand, protein, include_rings=False)
    assert [c.interaction_type for c in found] == ["hydrophobic"]


def test_specific_interaction_wins_over_hydrophobic_for_same_pair():
    """Одна пара атомов даёт один контакт — самый специфичный из подошедших."""
    ligand = parse_atoms(_pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"), pdbqt=True)
    protein = parse_atoms(
        "\n".join(
            [
                _pdb_line(100, "OD1", "ASP", 189, 2.6, 0.0, 0.0, "O"),
                _pdb_line(101, "ZN", "ZN", 301, -2.2, 0.0, 0.0, "ZN", record="HETATM"),
            ]
        ),
        pdbqt=False,
    )
    found = detect_contacts(ligand, protein, include_rings=False)
    by_residue = {c.residue_label: c.interaction_type for c in found}
    assert by_residue["Zn301"] == "metal"
    # Азот лиганда не распознан как заряженный, поэтому с Asp189 — водородная связь.
    assert by_residue["Asp189"] == "hydrogen_bond"
    assert len(found) == len({c.pair_key for c in found})


def test_distant_atoms_give_no_contacts():
    ligand = parse_atoms(_pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"), pdbqt=True)
    protein = parse_atoms(_pdb_line(100, "OD1", "ASP", 189, 50.0, 0.0, 0.0, "O"), pdbqt=False)
    assert detect_contacts(ligand, protein, include_rings=False) == []


def test_one_pi_cation_per_residue():
    """His даёт π-катион через ND1 и NE2; в отчёт идёт один, ближайший."""
    ring = "\n".join(
        _pdbqt_line(
            index + 1,
            f"C{index + 1}",
            "LIG",
            1,
            1.4 * math.cos(index * math.pi / 3),
            1.4 * math.sin(index * math.pi / 3),
            0.0,
            "A",
        )
        for index in range(6)
    )
    ligand = parse_atoms(ring, pdbqt=True)
    protein = parse_atoms(
        "\n".join(
            [
                _pdb_line(100, "ND1", "HIS", 94, 0.0, 0.0, 4.0, "N"),
                _pdb_line(101, "NE2", "HIS", 94, 0.0, 0.0, 4.6, "N"),
            ]
        ),
        pdbqt=False,
    )
    found = [c for c in detect_contacts(ligand, protein) if c.interaction_type == "pi_cation"]
    assert len(found) == 1
    assert found[0].distance_a == pytest.approx(4.0, abs=1e-3)


# --------------------------------------------------------------------------
# Сравнение поз
# --------------------------------------------------------------------------


def test_vina_deltas_are_relative_to_reference():
    analyses = [
        _analysis("pocket1", 1, -9.0),
        _analysis("pocket1", 2, -8.4),
        _analysis("pocket2", 1, -7.0),
    ]
    reference = PoseKey("lig", "pocket1", 1)
    deltas = vina_deltas(analyses, reference)
    by_pose = {d.key.pose_id if d.key.pocket_id == "pocket1" else d.key.pocket_id: d for d in deltas}
    assert by_pose[1].is_reference and by_pose[1].delta_to_reference == 0.0
    assert by_pose[2].delta_to_reference == pytest.approx(0.6)
    assert by_pose["pocket2"].delta_to_reference == pytest.approx(2.0)


def test_contact_changes_use_the_baseline_of_the_same_pocket():
    """Поза чужого кармана стоит в другом участке белка.

    Если считать её контакты относительно позы №1 всего прогона, «исчезло»
    окажется всё подряд, и смена участка притворится изменением взаимодействий.
    """
    ligand = parse_atoms(_pdbqt_line(1, "N1", "LIG", 1, 0.0, 0.0, 0.0, "NA"), pdbqt=True)
    asp = parse_atoms(_pdb_line(100, "OD1", "ASP", 189, 2.9, 0.0, 0.0, "O"), pdbqt=False)
    ser = parse_atoms(_pdb_line(200, "OG", "SER", 190, 2.9, 0.0, 0.0, "O"), pdbqt=False)

    with_asp = detect_contacts(ligand, asp, include_rings=False)
    with_ser = detect_contacts(ligand, ser, include_rings=False)
    analyses = [
        _analysis("pocket1", 1, -9.0, ligand, with_asp),
        _analysis("pocket1", 2, -8.0, ligand, with_ser),
        _analysis("pocket2", 1, -7.0, ligand, with_ser),
    ]
    matrix = contact_matrix(analyses, PoseKey("lig", "pocket1", 1))

    assert pocket_baselines(analyses) == {
        "pocket1": PoseKey("lig", "pocket1", 1),
        "pocket2": PoseKey("lig", "pocket2", 1),
    }
    changed_poses = {change.key.pose_id for change in matrix.changes}
    assert changed_poses == {2}, "базовая поза кармана сама с собой не сравнивается"
    kinds = {(c.residue_label, c.change) for c in matrix.changes}
    assert ("Asp189", "lost") in kinds
    assert ("Ser190", "gained") in kinds


def test_same_atom_pair_measured_in_every_pose():
    """PDF требует измерять одни и те же пары атомов во всех позах."""
    protein = parse_atoms(_pdb_line(100, "OD1", "ASP", 189, 0.0, 0.0, 0.0, "O"), pdbqt=False)
    near = parse_atoms(_pdbqt_line(7, "N1", "LIG", 1, 2.9, 0.0, 0.0, "NA"), pdbqt=True)
    far = parse_atoms(_pdbqt_line(7, "N1", "LIG", 1, 11.5, 0.0, 0.0, "NA"), pdbqt=True)

    reference = _analysis("pocket1", 1, -9.0, near, detect_contacts(near, protein, include_rings=False))
    moved = _analysis("pocket1", 2, -8.0, far, detect_contacts(far, protein, include_rings=False))
    assert moved.contacts == [], "на 11,5 Å контакта быть не должно"

    pairs = key_atom_pairs(reference)
    series = distance_series([reference, moved], protein, pairs)
    assert len(series) == 1
    values = series[0].values
    assert values["pocket1/поза 1"] == pytest.approx(2.9, abs=1e-3)
    # Контакта нет, а расстояние всё равно измерено: ряд не должен рваться.
    assert values["pocket1/поза 2"] == pytest.approx(11.5, abs=1e-3)


def test_missing_atom_gives_none_not_zero():
    """Пропуск нельзя писать нулём: ноль означал бы совпадение атомов."""
    protein = parse_atoms(_pdb_line(100, "OD1", "ASP", 189, 0.0, 0.0, 0.0, "O"), pdbqt=False)
    near = parse_atoms(_pdbqt_line(7, "N1", "LIG", 1, 2.9, 0.0, 0.0, "NA"), pdbqt=True)
    other = parse_atoms(_pdbqt_line(99, "N1", "LIG", 1, 2.9, 0.0, 0.0, "NA"), pdbqt=True)

    reference = _analysis("pocket1", 1, -9.0, near, detect_contacts(near, protein, include_rings=False))
    without = _analysis("pocket1", 2, -8.0, other, [])
    series = distance_series([reference, without], protein, key_atom_pairs(reference))
    assert series[0].values["pocket1/поза 2"] is None


def test_protein_state_table_has_a_row_per_pose_and_blank_for_missing():
    protein = parse_atoms(_pdb_line(100, "OD1", "ASP", 189, 0.0, 0.0, 0.0, "O"), pdbqt=False)
    near = parse_atoms(_pdbqt_line(7, "N1", "LIG", 1, 2.9, 0.0, 0.0, "NA"), pdbqt=True)
    reference = _analysis("pocket1", 1, -9.0, near, detect_contacts(near, protein, include_rings=False))
    second = _analysis("pocket1", 2, -8.0, near, [])
    deltas = vina_deltas([reference, second], reference.key)
    series = distance_series([reference, second], protein, key_atom_pairs(reference))

    fields, rows = protein_state_rows([reference, second], deltas, series)
    assert len(rows) == 2
    assert "n_salt_bridge" in fields and "contact_count" in fields
    assert rows[0]["vina_delta_to_reference"] == 0.0
    assert rows[1]["vina_delta_to_reference"] == pytest.approx(1.0)
    assert rows[0]["rmsd_to_crystal"] == "", "нет данных — пустая ячейка, не ноль"


def test_selection_keeps_reference_pocket_and_best_of_others():
    analyses = [_analysis("pocket1", index, -9.0 + index * 0.1) for index in range(1, 10)]
    analyses.append(_analysis("pocket2", 1, -5.0))
    analyses.append(_analysis("pocket2", 2, -6.5))
    analyses.append(_analysis("pocket3", 1, -4.0))
    reference = PoseKey("lig", "pocket1", 1)

    selected = select_for_comparison(analyses, reference, max_poses=11)
    keys = {(item.key.pocket_id, item.key.pose_id) for item in selected}
    assert len(selected) == 11
    assert ("pocket2", 2) in keys, "от чужого кармана берётся лучшая поза"
    assert ("pocket2", 1) not in keys
    assert ("pocket3", 1) in keys
    assert sum(1 for item in selected if item.key.pocket_id == "pocket1") == 9


def test_short_series_gives_no_findings():
    assert series_findings("ряд", ["a", "b"], [1.0, 9.0]) == []


def test_monotonic_series_reported_as_trend():
    found = series_findings("ряд", ["a", "b", "c", "d"], [1.0, 2.0, 3.0, 4.0])
    assert any(item.finding == "trend_up" for item in found)


def test_spike_detected_against_robust_baseline():
    labels = [f"p{index}" for index in range(8)]
    values = [3.0, 3.1, 2.9, 3.0, 12.0, 3.1, 2.95, 3.05]
    found = series_findings("расстояние", labels, values)
    spikes = [item for item in found if item.finding == "spike"]
    assert spikes and spikes[0].pose_label == "p4"
    assert spikes[0].value == pytest.approx(12.0)


def test_flat_series_has_no_spike():
    found = series_findings("ряд", [f"p{i}" for i in range(6)], [3.0] * 6)
    assert [item for item in found if item.finding == "spike"] == []


# --------------------------------------------------------------------------
# Живые прогоны сценария 2
# --------------------------------------------------------------------------


def _live_runs() -> list[Path]:
    if not _DOCKING.is_dir():
        return []
    return sorted(
        path.parent
        for path in _DOCKING.glob("*/poses.csv")
        if (path.parent / "run_meta.json").is_file() and (path.parent / "poses").is_dir()
    )


@pytest.mark.integration
def test_compare_poses_on_live_run():
    """Сквозной прогон сценария 3 на результате сценария 2."""
    runs = _live_runs()
    if not runs:
        pytest.skip("нет прогонов сценария 2 в outputs/docking")
    from models.request import ResearchRequest
    from orchestration.state import WorkflowState
    from agents.pose_comparison_agent import compare_poses

    run_dir = runs[0]
    state = WorkflowState(request=ResearchRequest(text="Сравни позы и покажи изменения контактов"))
    compare_poses(state, run_dir / "poses.csv")

    assert state.errors == []
    summary = state.observations["pose_comparison"]
    assert summary["compared"] is True
    for name in ("pose_contacts_csv", "protein_state_csv", "vina_by_pose_png", "contact_map_png"):
        assert Path(state.artifacts[name]).is_file()

    methods = {item.method for item in state.evidence}
    assert {"vina_score_by_pose", "vina_delta_to_best_pose", "pose_contact_count"} <= methods
    assert "pose_to_inspect" in methods


@pytest.mark.integration
def test_every_reported_distance_is_reproducible_from_coordinates():
    """Кейс требует, чтобы расстояние с графика сверялось по координатам вручную.

    Тест делает ровно эту сверку: берёт строку из pose_contacts.csv, открывает
    названный в ней файл позы и рецептор и пересчитывает расстояние по номерам
    атомов. Если в отчёт попадёт число, которого нет в координатах, тест упадёт.
    """
    import csv
    import json

    runs = _live_runs()
    if not runs:
        pytest.skip("нет прогонов сценария 2 в outputs/docking")
    from models.request import ResearchRequest
    from orchestration.state import WorkflowState
    from agents.pose_comparison_agent import compare_poses

    run_dir = runs[0]
    state = WorkflowState(request=ResearchRequest(text="Сравни позы"))
    compare_poses(state, run_dir / "poses.csv")
    assert state.errors == []

    meta = json.loads((run_dir / "run_meta.json").read_text(encoding="utf-8"))
    protein = {atom.serial: atom for atom in read_structure(meta["receptor_pdb_path"])}

    with Path(state.artifacts["pose_contacts_csv"]).open(encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if not row["ligand_atom"].startswith("ring(")]
    assert rows, "в отчёте нет ни одного атом-атомного контакта"

    checked = 0
    for row in rows[:40]:
        pose_path = Path(row["pose_file"])
        if not pose_path.is_absolute():
            pose_path = _REPO_ROOT / pose_path
        ligand = {atom.serial: atom for atom in read_structure(pose_path)}
        ligand_atom = ligand[int(row["ligand_atom_serial"])]
        protein_atom = protein[int(row["protein_atom_serial"])]
        assert ligand_atom.distance_to(protein_atom) == pytest.approx(
            float(row["distance_a"]), abs=1e-3
        )
        assert protein_atom.res_name == row["residue_name"]
        assert protein_atom.res_seq == int(row["residue_seq"])
        checked += 1
    assert checked >= 1


@pytest.mark.integration
def test_single_pose_table_does_not_invent_a_comparison():
    """На одной позе сценарий не выдумывает сравнение, а говорит, чего не хватает."""
    runs = _live_runs()
    if not runs:
        pytest.skip("нет прогонов сценария 2 в outputs/docking")
    import shutil
    import tempfile

    from models.docking import load_poses_csv, write_poses_csv
    from models.request import ResearchRequest
    from orchestration.state import WorkflowState
    from agents.pose_comparison_agent import compare_poses

    source = runs[0]
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "run"
        shutil.copytree(source, target)
        rows = load_poses_csv(target / "poses.csv")
        write_poses_csv(target / "poses.csv", rows[:1])

        state = WorkflowState(request=ResearchRequest(text="Сравни позы"))
        compare_poses(state, target / "poses.csv")

        assert state.observations["pose_comparison"]["compared"] is False
        assert any("одна поза" in text for text in state.limitations)
        assert "contact_map_png" not in state.artifacts
