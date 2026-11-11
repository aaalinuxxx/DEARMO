import pytest
from pathlib import Path

from models.docking import load_poses_csv

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "data" / "examples" / "3ptb_ben" / "poses.csv"


@pytest.mark.integration
def test_ben_redock_rmsd_golden():
    if not CSV.is_file():
        pytest.skip("no golden file")
    rows = load_poses_csv(CSV)
    assert rows[0].rmsd_to_crystal == pytest.approx(0.356, abs=0.01)
