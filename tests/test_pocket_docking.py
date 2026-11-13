"""Сценарий 2: карманы, все позы, поза №1 и схема таблицы.

Тесты намеренно не требуют установленных Vina и P2Rank: разбор вывода и схема
таблицы проверяются на зафиксированных образцах настоящего вывода, а живой
расчёт вынесен в отдельный integration-тест.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from models.docking import (
    POSE_CSV_FIELDS,
    PocketRow,
    PoseRow,
    best_pose_overall,
    load_pockets_csv,
    load_poses_csv,
    write_pockets_csv,
    write_poses_csv,
)
from tools.p2rank_runner import parse_predictions
from tools.pdbqt_prep import ad_type, clean_receptor_pdb, receptor_pdbqt
from tools.vina_runner import parse_vina_modes, split_pdbqt_models

ROOT = Path(__file__).resolve().parent.parent

# Настоящий вывод P2Rank 2.5.1 на 3PTB. Имена колонок дополнены пробелами
# ровно так, как их печатает сам P2Rank.
P2RANK_CSV = """\
name     ,  rank,   score, probability, sas_points, surf_atoms,   center_x,   center_y,   center_z, residue_ids, surf_atom_ids
pocket1  ,     1,    2.36,       0.062,         23,         18,    -1.5173,    14.4714,    17.4739, A_189 A_190 A_195 A_57, 292 1252
pocket2  ,     2,    1.95,       0.041,         20,         17,     9.4208,     3.3601,    11.7224, A_129 A_130 A_162, 812 815
pocket3  ,     3,    1.31,       0.014,         27,         14,    -7.6220,    15.0873,    28.3686, A_141 A_151 A_39, 154 159
"""

# Настоящая таблица режимов AutoDock Vina.
VINA_STDOUT = """\
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -6.155          0          0
   2       -5.267      2.658      3.617
   3       -5.201      2.494      3.333
"""

PDB_WITH_ZN = """\
ATOM      1  N   HIS A  94      10.000  20.000  30.000  1.00  0.00           N
ATOM      2  CA  HIS A  94      11.000  20.000  30.000  1.00  0.00           C
ATOM      3  N   HIS B  94      40.000  20.000  30.000  1.00  0.00           N
HETATM    4 ZN    ZN A 262      12.000  21.000  31.000  1.00  0.00          ZN
HETATM    5  O   HOH A 400      50.000  50.000  50.000  1.00  0.00           O
HETATM    6  C1  AZM A 300      13.000  21.000  31.000  1.00  0.00           C
END
"""


# --- P2Rank -----------------------------------------------------------------


def test_parse_predictions_reads_rank_score_and_center():
    pockets = parse_predictions(_write(P2RANK_CSV, "predictions.csv"))
    assert [p.pocket_id for p in pockets] == ["pocket1", "pocket2", "pocket3"]
    assert [p.rank for p in pockets] == [1, 2, 3]
    assert pockets[0].p2rank_score == pytest.approx(2.36)
    assert pockets[0].center == pytest.approx((-1.5173, 14.4714, 17.4739))
    assert pockets[0].surf_atom_count == 18
    assert "A_189" in pockets[0].residue_ids


def test_pockets_are_sorted_by_rank_not_by_file_order():
    shuffled = "\n".join(
        [P2RANK_CSV.splitlines()[0], P2RANK_CSV.splitlines()[3], P2RANK_CSV.splitlines()[1]]
    )
    pockets = parse_predictions(_write(shuffled + "\n", "shuffled.csv"))
    assert [p.rank for p in pockets] == [1, 3]


# --- Vina -------------------------------------------------------------------


def test_parse_vina_modes_reads_every_pose_with_both_rmsd_bounds():
    modes = parse_vina_modes(VINA_STDOUT)
    assert [m.mode for m in modes] == [1, 2, 3]
    assert modes[0].affinity == pytest.approx(-6.155)
    # Первая поза — точка отсчёта, обе оценки RMSD нулевые по построению.
    assert modes[0].rmsd_lb == 0.0 and modes[0].rmsd_ub == 0.0
    assert modes[1].rmsd_lb == pytest.approx(2.658)
    assert modes[1].rmsd_ub == pytest.approx(3.617)


def test_parse_vina_modes_ignores_header_and_separator_rows():
    assert parse_vina_modes("mode | affinity\n-----+------\n") == []


def test_split_pdbqt_models_writes_one_file_per_pose(tmp_path: Path):
    multi = tmp_path / "out.pdbqt"
    multi.write_text(
        "MODEL 1\nATOM      1  C   UNL     1       1.000   0.000   0.000  0.00  0.00     C\nENDMDL\n"
        "MODEL 2\nATOM      1  C   UNL     1       2.000   0.000   0.000  0.00  0.00     C\nENDMDL\n",
        encoding="utf-8",
    )
    files = split_pdbqt_models(multi, tmp_path / "poses", "pocket1")
    assert [f.name for f in files] == ["pocket1_pose1.pdbqt", "pocket1_pose2.pdbqt"]
    assert "MODEL" not in files[0].read_text(encoding="utf-8")
    assert "2.000" in files[1].read_text(encoding="utf-8")


# --- Рецептор ---------------------------------------------------------------


def test_metal_survives_receptor_cleanup_but_water_and_ligand_do_not():
    """Zn²⁺ — часть активного центра, а не примесь кристаллизации."""
    cleaned = clean_receptor_pdb(PDB_WITH_ZN)
    assert "ZN" in cleaned
    assert "HOH" not in cleaned
    assert "AZM" not in cleaned


def test_chain_filter_keeps_only_requested_chain():
    cleaned = clean_receptor_pdb(PDB_WITH_ZN, chain="A")
    assert " A  94" in cleaned
    assert " B  94" not in cleaned


def test_metal_keeps_its_autodock_type_instead_of_becoming_carbon():
    assert ad_type("ZN") == "Zn"
    assert ad_type("FE") == "Fe"
    assert ad_type("C") == "C"
    assert "Zn" in receptor_pdbqt(PDB_WITH_ZN)


# --- Таблица поз ------------------------------------------------------------


def test_pose_csv_schema_is_the_agreed_contract():
    """Сценарии 3 и 4 читают эти колонки по именам. Переименование — поломка."""
    assert POSE_CSV_FIELDS == (
        "ligand_id",
        "pocket_id",
        "pocket_rank",
        "p2rank_score",
        "pocket_center_x",
        "pocket_center_y",
        "pocket_center_z",
        "pose_id",
        "vina_score",
        "pose_file",
        "rmsd_lb",
        "rmsd_ub",
        "rmsd_to_first_pose",
        "rmsd_to_crystal",
    )


def test_poses_round_trip_through_csv(tmp_path: Path):
    rows = [
        PoseRow(
            ligand_id="benzamidine",
            pocket_id="pocket1",
            pocket_rank=1,
            p2rank_score=2.36,
            pose_id=1,
            vina_score=-6.155,
            pose_file="poses/pocket1_pose1.pdbqt",
            rmsd_lb=0.0,
            rmsd_ub=0.0,
            rmsd_to_first_pose=0.0,
            rmsd_to_crystal=0.411,
            pocket_center=(-1.5173, 14.4714, 17.4739),
        )
    ]
    back = load_poses_csv(write_poses_csv(tmp_path / "poses.csv", rows))
    assert back[0].vina_score == pytest.approx(-6.155)
    assert back[0].pocket_center == pytest.approx((-1.5173, 14.4714, 17.4739))
    assert back[0].rmsd_ub == 0.0


def test_missing_measurement_reads_as_none_not_zero(tmp_path: Path):
    """Пустая ячейка — это «не измерено». Ноль означал бы идеальное совпадение."""
    path = tmp_path / "poses.csv"
    path.write_text(
        ",".join(POSE_CSV_FIELDS) + "\n" + "lig,pocket1,1,,,,,1,-6.1,p.pdbqt,,,,\n",
        encoding="utf-8",
    )
    row = load_poses_csv(path)[0]
    assert row.rmsd_to_crystal is None
    assert row.rmsd_ub is None
    assert row.pocket_center is None


def test_old_short_schema_still_loads():
    """Таблицы в data/examples записаны до появления новых колонок."""
    legacy = ROOT / "data" / "examples" / "3ptb_ben" / "poses.csv"
    if not legacy.is_file():
        pytest.skip("нет сохранённой таблицы")
    rows = load_poses_csv(legacy)
    assert rows[0].pose_id == 1
    # Прежнее имя колонки продолжает читаться и попадает в новое поле.
    assert rows[0].rmsd_lb == rows[0].rmsd_to_first_pose


def test_best_pose_is_global_minimum_across_pockets():
    """Поза №1 — лучшая по всему прогону, а не лучшая в первом кармане."""
    rows = [
        PoseRow("lig", "pocket1", 1, 2.4, 1, -6.1, "a.pdbqt"),
        PoseRow("lig", "pocket2", 2, 1.9, 1, -7.3, "b.pdbqt"),
        PoseRow("lig", "pocket2", 2, 1.9, 2, -5.0, "c.pdbqt"),
    ]
    best = best_pose_overall(rows)
    assert best is not None
    assert (best.pocket_id, best.pose_id, best.vina_score) == ("pocket2", 1, -7.3)


def test_best_pose_tie_breaks_on_pocket_rank():
    """При равных оценках выигрывает карман с лучшим рангом, а не порядок строк."""
    rows = [
        PoseRow("lig", "pocket3", 3, 1.3, 1, -6.5, "c.pdbqt"),
        PoseRow("lig", "pocket1", 1, 2.4, 1, -6.5, "a.pdbqt"),
    ]
    best = best_pose_overall(rows)
    assert best is not None and best.pocket_id == "pocket1"


def test_best_pose_of_empty_table_is_none():
    assert best_pose_overall([]) is None


def test_pockets_round_trip_through_csv(tmp_path: Path):
    rows = [PocketRow("pocket1", 1, 2.36, (-1.5, 14.5, 17.5), 18, ("A_189", "A_195"))]
    back = load_pockets_csv(write_pockets_csv(tmp_path / "pockets.csv", rows))
    assert back[0].residue_ids == ("A_189", "A_195")
    assert back[0].center == pytest.approx((-1.5, 14.5, 17.5))


_TMP: Path | None = None


def _write(text: str, name: str) -> Path:
    """Временный файл для разборщиков, принимающих путь."""
    import tempfile

    global _TMP
    if _TMP is None:
        _TMP = Path(tempfile.mkdtemp(prefix="drugdis_pocket_"))
    path = _TMP / name
    path.write_text(text, encoding="utf-8")
    return path
