"""Чтение/запись артефактов докинга."""
from __future__ import annotations

from pathlib import Path

from models.docking import POSE_CSV_FIELDS, PoseRow, load_poses_csv


# Что искать в каталоге прогона и под каким именем класть в artifacts.
# Рецептор попадает сюда отдельно: проверке геометрии (сценарий 4) нужны
# координаты белка, а раньше при чтении сохранённого прогона на них ничто
# не указывало и путь приходилось угадывать.
_RUN_FILES = {
    "report_md": "report.md",
    "poses_csv": "poses.csv",
    "pockets_csv": "pockets.csv",
    "receptor_pdbqt": "receptor.pdbqt",
    "receptor_pdb": "receptor.pdb",
    "run_meta": "run_meta.json",
}


def load_run_directory(example_dir: Path) -> dict[str, str]:
    example_dir = Path(example_dir)
    artifacts: dict[str, str] = {}
    for key, name in _RUN_FILES.items():
        path = example_dir / name
        if path.is_file():
            artifacts[key] = str(path)
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
