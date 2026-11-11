"""Модели строк таблицы докинга и загрузки CSV."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

POSE_CSV_FIELDS = (
    "ligand_id",
    "pocket_id",
    "pocket_rank",
    "p2rank_score",
    "pose_id",
    "vina_score",
    "pose_file",
    "rmsd_to_first_pose",
    "rmsd_to_crystal",
)


@dataclass
class PoseRow:
    ligand_id: str
    pocket_id: str
    pocket_rank: int
    p2rank_score: float | None
    pose_id: int
    vina_score: float
    pose_file: str
    rmsd_to_first_pose: float | None = None
    rmsd_to_crystal: float | None = None

    @classmethod
    def from_dict(cls, row: dict[str, str]) -> "PoseRow":
        def _float(key: str) -> float | None:
            raw = row.get(key, "").strip()
            if not raw:
                return None
            return float(raw)

        def _int(key: str) -> int:
            return int(float(row[key]))

        return cls(
            ligand_id=row["ligand_id"],
            pocket_id=row["pocket_id"],
            pocket_rank=_int("pocket_rank"),
            p2rank_score=_float("p2rank_score"),
            pose_id=_int("pose_id"),
            vina_score=float(row["vina_score"]),
            pose_file=row["pose_file"],
            rmsd_to_first_pose=_float("rmsd_to_first_pose"),
            rmsd_to_crystal=_float("rmsd_to_crystal"),
        )


def load_poses_csv(path: Path | str) -> list[PoseRow]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    rows: list[PoseRow] = []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(PoseRow.from_dict(row))
    return rows


def best_vina_per_ligand_pocket(rows: list[PoseRow]) -> dict[tuple[str, str], PoseRow]:
    best: dict[tuple[str, str], PoseRow] = {}
    for row in rows:
        key = (row.ligand_id, row.pocket_id)
        if key not in best or row.vina_score < best[key].vina_score:
            best[key] = row
    return best
