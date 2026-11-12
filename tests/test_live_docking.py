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


def test_new_pair_calls_live_vina_not_acetazolamide_cache(monkeypatch, tmp_path: Path):
    captured: dict = {}

    class _Result:
        vina_score = -4.2
        pose_file = tmp_path / "pose_1.pdbqt"
        poses_csv = tmp_path / "poses.csv"
        receptor_pdbqt = tmp_path / "receptor.pdbqt"
        box_source = "crystal_ligand"
        ligand_prep = "rdkit-rigid"

    def fake_dock(pdb_id, smiles, ligand_id, output_dir, exhaustiveness=8):
        captured["pdb_id"] = pdb_id
        captured["smiles"] = smiles
        captured["ligand_id"] = ligand_id
        _Result.poses_csv.write_text(
            "ligand_id,pocket_id,pocket_rank,p2rank_score,pose_id,vina_score,pose_file,rmsd_to_first_pose,rmsd_to_crystal\n"
            f"aspirin,site1,1,,1,-4.2,{_Result.pose_file},0.0,\n",
            encoding="utf-8",
        )
        return _Result()

    monkeypatch.setattr(docking_agent, "dock_pair", fake_dock)
    state = WorkflowState(request=ResearchRequest(text="Выполни докинг аспирина в 1OXR"))
    state.observations.update(
        {
            "pdb_id": "1OXR",
            "ligand_query": "aspirin",
            "ligand_id": "aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        }
    )
    state = docking_agent.run_or_load(state, Path("data/examples/acetazolamide"))
    assert captured["pdb_id"] == "1OXR"
    assert captured["smiles"].startswith("CC(=O)O")
    assert state.artifacts["poses_csv"] == str(tmp_path / "poses.csv")
    scores = [ev.value for ev in state.evidence if ev.method == "vina_score"]
    assert scores == [-4.2]
    assert -7.82 not in scores
