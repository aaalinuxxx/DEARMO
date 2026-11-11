from pathlib import Path

from tools.pose_validation import pass_rate, validate_pose_file

ROOT = Path(__file__).resolve().parent.parent
POSE = ROOT / "data" / "examples" / "acetazolamide" / "poses" / "acetazolamide_pocket1_pose1.pdbqt"


def test_basic_pose_validation():
    checks = validate_pose_file(POSE)
    assert pass_rate(checks) >= 50.0
