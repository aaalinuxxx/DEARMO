"""Чтение/запись артефактов докинга."""
from __future__ import annotations

from pathlib import Path

from models.docking import POSE_CSV_FIELDS, PoseRow, load_poses_csv


def load_run_directory(example_dir: Path) -> dict[str, str]:
    example_dir = Path(example_dir)
    artifacts: dict[str, str] = {}
    report = example_dir / "report.md"
    poses = example_dir / "poses.csv"
    if report.is_file():
        artifacts["report_md"] = str(report)
    if poses.is_file():
        artifacts["poses_csv"] = str(poses)
    return artifacts


def read_report_summary(path: Path) -> str:
    path = Path(path)
    if not path.is_file():
        return ""
    text = path.read_text(encoding="utf-8")
    return text[:2000]


def poses_from_dir(example_dir: Path) -> list[PoseRow]:
    csv_path = example_dir / "poses.csv"
    return load_poses_csv(csv_path)
