"""Живой Vina вызывается для новой пары и не подменяет её сохранённым прогоном."""
from pathlib import Path

from models.request import ResearchRequest
from orchestration.state import WorkflowState
from tools.pdbqt_prep import crystal_ligand_center, receptor_pdbqt
from agents import docking_agent


_PDB = """\
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  20.000  30.000  1.00  0.00           C
HETATM    3  O   HOH A   9      50.000  50.000  50.000  1.00  0.00           O
HETATM    4  C1  LIG A   2      12.000  21.000  31.000  1.00  0.00           C
HETATM    5  C2  LIG A   2      13.000  21.000  31.000  1.00  0.00           C
HETATM    6  N1  LIG A   2      12.500  22.000  31.000  1.00  0.00           N
END
"""


def test_receptor_drops_crystal_ligand_and_water():
    text = receptor_pdbqt(_PDB)
    assert "ALA" in text
    assert "HOH" not in text
    assert "LIG" not in text
    center = crystal_ligand_center(_PDB)
    assert center is not None
    assert abs(center[0] - 12.5) < 0.01


def _fake_run(tmp_path: Path, control: dict | None = None):
    """Результат докинга с двумя карманами: поза №1 лежит не в первом."""
    from models.docking import PocketRow, PoseRow
    from tools.live_docking import DockingRun

    pockets = [
        PocketRow("pocket1", 1, 2.4, (1.0, 2.0, 3.0), 10, ("A_189",)),
        PocketRow("pocket2", 2, 1.9, (4.0, 5.0, 6.0), 8, ("A_57",)),
    ]
    poses = [
        PoseRow("aspirin", "pocket1", 1, 2.4, 1, -4.2, str(tmp_path / "p1.pdbqt"), pocket_center=pockets[0].center),
        PoseRow("aspirin", "pocket2", 2, 1.9, 1, -5.1, str(tmp_path / "p2.pdbqt"), pocket_center=pockets[1].center),
    ]
    return DockingRun(
        pdb_id="1OXR",
        ligand_id="aspirin",
        pdb_path=tmp_path / "1OXR.pdb",
        receptor_pdb=tmp_path / "receptor.pdb",
        receptor_pdbqt=tmp_path / "receptor.pdbqt",
        ligand_pdbqt=tmp_path / "ligand.pdbqt",
        pockets=pockets,
        poses=poses,
        best=poses[1],
        poses_csv=tmp_path / "poses.csv",
        pockets_csv=tmp_path / "pockets.csv",
        run_meta=tmp_path / "run_meta.json",
        box_source="p2rank",
        box_size=(22.0, 22.0, 22.0),
        ligand_prep="rdkit-rigid",
        exhaustiveness=8,
        num_modes=9,
        seed=42,
        control=control or {},
    )


def test_new_pair_calls_live_vina_not_acetazolamide_cache(monkeypatch, tmp_path: Path):
    captured: dict = {}

    def fake_dock(pdb_id, smiles, ligand_id, output_dir, **kwargs):
        captured.update(pdb_id=pdb_id, smiles=smiles, ligand_id=ligand_id, kwargs=kwargs)
        return _fake_run(tmp_path)

    monkeypatch.setattr(docking_agent, "dock_pair", fake_dock)
    state = WorkflowState(request=ResearchRequest(text="Выполни докинг аспирина в 1OXR"))
    state.observations.update(
        {
            "pdb_id": "1OXR",
            "pdb_chain": "A",
            "ligand_query": "aspirin",
            "ligand_id": "aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "ligand_comp_ids": ["ASP"],
        }
    )
    state = docking_agent.run_or_load(state, Path("data/examples/acetazolamide"))
    assert captured["pdb_id"] == "1OXR"
    assert captured["smiles"].startswith("CC(=O)O")
    # Цепь и код эталона доходят до докинга, но эталон нужен только для контроля.
    assert captured["kwargs"]["chain"] == "A"
    assert captured["kwargs"]["crystal_comp_id"] == "ASP"
    assert state.artifacts["poses_csv"] == str(tmp_path / "poses.csv")
    assert state.observations["poses_csv_path"] == str(tmp_path / "poses.csv")
    assert state.observations["receptor_pdb_path"] == str(tmp_path / "receptor.pdb")
    scores = [ev.value for ev in state.evidence if ev.method == "vina_score"]
    # Поза №1 — глобальный минимум по обоим карманам, а не первая строка.
    assert scores == [-5.1]
    assert -7.82 not in scores
    best_ctx = next(ev.context for ev in state.evidence if ev.method == "vina_score")
    assert best_ctx["pocket_id"] == "pocket2"


def test_every_pocket_is_reported_with_rank_score_and_center(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(docking_agent, "dock_pair", lambda *a, **k: _fake_run(tmp_path))
    state = WorkflowState(request=ResearchRequest(text="Найди карманы"))
    state.observations.update({"pdb_id": "1OXR", "ligand_id": "aspirin", "smiles": "CC(=O)O"})
    docking_agent.run_or_load(state)
    pockets = [ev for ev in state.evidence if ev.method == "p2rank_pocket"]
    assert [ev.context["rank"] for ev in pockets] == [1, 2]
    assert pockets[0].value == 2.4
    assert pockets[0].context["center"] == [1.0, 2.0, 3.0]
    assert pockets[0].source == "P2Rank"


def test_control_metrics_are_reported_only_when_a_reference_existed(monkeypatch, tmp_path: Path):
    """Без эталонного лиганда RMSD и ранг кармана не выдумываются."""
    monkeypatch.setattr(docking_agent, "dock_pair", lambda *a, **k: _fake_run(tmp_path))
    state = WorkflowState(request=ResearchRequest(text="x"))
    state.observations.update({"pdb_id": "1OXR", "ligand_id": "aspirin", "smiles": "CC(=O)O"})
    docking_agent.run_or_load(state)
    methods = {ev.method for ev in state.evidence}
    assert "rmsd_to_crystal" not in methods
    assert "pocket_hit_rank" not in methods

    control = {
        "available": True,
        "crystal_comp_id": "ASP",
        "pocket_hit_rank": 2,
        "pocket_hit_distance_a": 1.234,
        "pocket_hit_radius_a": 4.0,
        "best_pose_rmsd_to_crystal_a": 0.411,
        "rmsd_threshold_a": 2.0,
    }
    monkeypatch.setattr(docking_agent, "dock_pair", lambda *a, **k: _fake_run(tmp_path, control))
    state = WorkflowState(request=ResearchRequest(text="x"))
    state.observations.update({"pdb_id": "1OXR", "ligand_id": "aspirin", "smiles": "CC(=O)O"})
    docking_agent.run_or_load(state)
    by_method = {ev.method: ev for ev in state.evidence}
    assert by_method["pocket_hit_rank"].value == 2
    assert by_method["rmsd_to_crystal"].value == 0.411
    assert by_method["rmsd_to_crystal"].unit == "Å"
