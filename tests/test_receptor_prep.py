"""Подготовка рецептора и честный вердикт контрольного редокинга."""
from __future__ import annotations

from types import SimpleNamespace

from agents import docking_agent
from models.request import ResearchRequest
from orchestration.state import WorkflowState
from tools import pdbqt_prep


def _pqr_line(serial: int, name: str, resname: str, chain: str, resnum: str) -> str:
    # Строка в том виде, как её пишет pdb2pqr: фиксированные колонки.
    return f"ATOM  {serial:5d}  {name:<3} {resname:>3} {chain}{resnum:>4}     -18.878 200.796  11.286 -0.4157 1.8240"


def test_pqr_with_four_digit_residue_numbers_is_split_into_fields():
    """Номер от 1000 слипался с цепью («A1007»), и Meeko падал на разборе."""
    text, with_icode = pdbqt_prep._free_format_pqr(_pqr_line(4568, "N", "MET", "A", "1007") + "\n")
    tokens = text.split()
    assert tokens[:6] == ["ATOM", "4568", "N", "MET", "A", "1007"]
    assert with_icode == 0


def test_insertion_code_becomes_its_own_field_and_numbering_is_kept():
    """«184A» Meeko не читает слитно, но принимает вставочный код отдельным полем."""
    line = "ATOM   1196  N   GLY A 184A    -18.878 200.796  11.286 -0.4157 1.8240"
    text, with_icode = pdbqt_prep._free_format_pqr(line + "\n")
    tokens = text.split()
    assert tokens[4:7] == ["A", "184", "A"], "номер остатка не меняется, код идёт отдельно"
    assert with_icode == 1


def test_metal_ions_are_returned_to_the_receptor():
    """pdb2pqr отбрасывает ионы; без Zn докинг в карбоангидразу считал бы не ту задачу."""
    pdb = (
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
        "HETATM  500 ZN    ZN A 301      10.500  11.500  12.500  1.00  0.00          ZN\n"
    )
    lines = pdbqt_prep._metal_lines(pdb, start_serial=100)
    assert len(lines) == 1
    assert lines[0].startswith("HETATM")
    assert lines[0].rstrip().endswith("Zn"), "тип AutoDock металла — Zn, не углерод"
    assert "10.500" in lines[0]


def _result(rmsd: float) -> SimpleNamespace:
    return SimpleNamespace(
        best=SimpleNamespace(pose_file="pose.pdbqt"),
        control={
            "available": True,
            "crystal_comp_id": "AZM",
            "pocket_hit_rank": 1,
            "pocket_hit_distance_a": 0.7,
            "pocket_hit_radius_a": 4.0,
            "best_pose_rmsd_to_crystal_a": rmsd,
            "rmsd_threshold_a": 2.0,
        },
    )


def test_failed_redock_check_is_said_in_words():
    """RMSD 7.9 Å без вердикта выглядел как обычное число рядом с «success»."""
    state = WorkflowState(request=ResearchRequest(text="докинг"))
    docking_agent._record_control(state, _result(7.9))
    assert state.observations["redock_check"] == "failed"
    assert any("Контроль редокинга не пройден" in lim for lim in state.limitations)


def test_passed_redock_check_does_not_overclaim():
    state = WorkflowState(request=ResearchRequest(text="докинг"))
    docking_agent._record_control(state, _result(1.3))
    assert state.observations["redock_check"] == "passed"
    assert any("не гарантия" in lim for lim in state.limitations)


def test_label_layout_terminates_when_curves_converge(tmp_path):
    """Сходящиеся кривые зацикливали раскладку подписей: демо висело на 100% CPU.

    Три ряда с одинаковым значением в конце линии заставляют сдвигать подписи
    друг от друга. Прежний цикл «пока ближе gap» не завершался из-за плавающей
    запятой. Тест реально рисует график и падает по таймеру, а не зависает.
    """
    import signal

    from tools import plotting

    if plotting.plt is None:
        import pytest

        pytest.skip("matplotlib не установлен")

    class _Series:
        def __init__(self, label):
            self.label = self.residue_label = label
            # 2.0 Å — значение, на котором прежний цикл гарантированно не завершался:
            # (2.0 + 0.4) - 2.0 меньше 0.4 в двоичной арифметике.
            self.values = {"pose1": 2.0, "pose2": 2.0, "pose3": 2.0}

    def _timeout(*_):
        raise AssertionError("раскладка подписей зациклилась")

    previous = signal.signal(signal.SIGALRM, _timeout)
    signal.alarm(30)
    try:
        out = plotting.write_distance_series(
            [_Series("Thr199"), _Series("His94"), _Series("Zn")],
            ["pose1", "pose2", "pose3"],
            tmp_path / "series.png",
        )
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    assert out.is_file() and out.stat().st_size > 0


def test_p2rank_output_is_labelled_as_prediction_not_calculation():
    """Балл кармана выдаёт обученная модель: это предсказание. Vina при этом остаётся расчётом."""
    from types import SimpleNamespace

    from agents import pocket_worker
    from models.docking import PocketRow

    state = WorkflowState(request=ResearchRequest(text="докинг"))
    pocket = PocketRow(pocket_id="pocket1", rank=1, p2rank_score=10.45, center=(1.0, 2.0, 3.0),
                       surf_atom_count=30, residue_ids=("A_1",))
    docking_agent._record_pockets(state, SimpleNamespace(pockets=[pocket]))

    kinds = {ev.method: ev.kind for ev in state.evidence}
    assert kinds["p2rank_pocket"] == "prediction"
    assert "машинного обучения" in state.evidence[0].context["note"]
