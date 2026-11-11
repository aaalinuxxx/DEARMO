from pathlib import Path

from models.docking import best_vina_per_ligand_pocket, load_poses_csv

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "data" / "examples" / "acetazolamide" / "poses.csv"


def test_load_poses_and_best_per_pocket():
    rows = load_poses_csv(CSV)
    assert len(rows) >= 6
    best = best_vina_per_ligand_pocket(rows)
    assert best[("acetazolamide", "pocket1")].vina_score == -7.82


def test_ligand_order_changes_between_pockets():
    rows = load_poses_csv(CSV)
    best = best_vina_per_ligand_pocket(rows)
    ligands = sorted({k[0] for k in best})

    def rank(pocket: str) -> list[str]:
        return sorted(ligands, key=lambda lig: best[(lig, pocket)].vina_score)

    assert rank("pocket1") != rank("pocket2")
