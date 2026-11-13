"""Модели строк таблицы докинга и загрузки CSV.

Схема poses.csv — контракт между сценарием 2 (докинг, производит таблицу) и
сценарием 3 (сравнение поз, читает её). Одна строка — одна поза одного кармана.
Пустая ячейка означает «не измерено» и читается как None: ноль вместо пропуска
превращал бы отсутствие расчёта в идеальное совпадение.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

POSE_CSV_FIELDS = (
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

# Схема pockets.csv — выход P2Rank до докинга. Нужна, чтобы ранг и балл кармана
# можно было показать даже для кармана, в котором докинг не дал ни одной позы.
POCKET_CSV_FIELDS = (
    "pocket_id",
    "rank",
    "p2rank_score",
    "center_x",
    "center_y",
    "center_z",
    "surf_atom_count",
    "residue_ids",
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
    # Нижняя и верхняя оценки RMSD относительно позы №1 того же кармана,
    # как их считает сам Vina (столбцы rmsd l.b. / rmsd u.b.).
    rmsd_lb: float | None = None
    rmsd_ub: float | None = None
    # Прежнее имя той же величины: равен rmsd_lb. Оставлен, чтобы уже записанные
    # таблицы и читающий их код не пришлось переписывать одним коммитом.
    rmsd_to_first_pose: float | None = None
    rmsd_to_crystal: float | None = None
    pocket_center: tuple[float, float, float] | None = None

    @classmethod
    def from_dict(cls, row: dict[str, str]) -> "PoseRow":
        """Прочитать строку, не требуя колонок, появившихся позже.

        Таблицы в data/examples записаны по короткой схеме без центра кармана и
        без оценок RMSD от Vina. Они должны читаться и дальше, иначе сохранённые
        прогоны пришлось бы перегенерировать ради одного нового столбца.
        """

        def _float(key: str) -> float | None:
            raw = (row.get(key) or "").strip()
            if not raw:
                return None
            try:
                return float(raw)
            except ValueError:
                return None

        def _int(key: str, default: int = 0) -> int:
            raw = (row.get(key) or "").strip()
            if not raw:
                return default
            try:
                return int(float(raw))
            except ValueError:
                return default

        center_vals = (_float("pocket_center_x"), _float("pocket_center_y"), _float("pocket_center_z"))
        center = center_vals if all(v is not None for v in center_vals) else None

        rmsd_lb = _float("rmsd_lb")
        legacy_first = _float("rmsd_to_first_pose")
        if rmsd_lb is None:
            rmsd_lb = legacy_first
        if legacy_first is None:
            legacy_first = rmsd_lb

        return cls(
            ligand_id=row["ligand_id"],
            pocket_id=row["pocket_id"],
            pocket_rank=_int("pocket_rank"),
            p2rank_score=_float("p2rank_score"),
            pose_id=_int("pose_id", default=1),
            vina_score=float(row["vina_score"]),
            pose_file=row["pose_file"],
            rmsd_lb=rmsd_lb,
            rmsd_ub=_float("rmsd_ub"),
            rmsd_to_first_pose=legacy_first,
            rmsd_to_crystal=_float("rmsd_to_crystal"),
            pocket_center=center,  # type: ignore[arg-type]
        )

    def to_dict(self) -> dict[str, object]:
        """Строка для csv.DictWriter по POSE_CSV_FIELDS. None становится пустой ячейкой."""
        cx, cy, cz = self.pocket_center if self.pocket_center else (None, None, None)
        values: dict[str, object] = {
            "ligand_id": self.ligand_id,
            "pocket_id": self.pocket_id,
            "pocket_rank": self.pocket_rank,
            "p2rank_score": self.p2rank_score,
            "pocket_center_x": cx,
            "pocket_center_y": cy,
            "pocket_center_z": cz,
            "pose_id": self.pose_id,
            "vina_score": self.vina_score,
            "pose_file": self.pose_file,
            "rmsd_lb": self.rmsd_lb,
            "rmsd_ub": self.rmsd_ub,
            "rmsd_to_first_pose": self.rmsd_to_first_pose,
            "rmsd_to_crystal": self.rmsd_to_crystal,
        }
        return {key: ("" if values[key] is None else values[key]) for key in POSE_CSV_FIELDS}


@dataclass
class PocketRow:
    """Карман-кандидат от P2Rank, до всякого докинга."""

    pocket_id: str
    rank: int
    p2rank_score: float | None
    center: tuple[float, float, float]
    surf_atom_count: int | None = None
    residue_ids: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "pocket_id": self.pocket_id,
            "rank": self.rank,
            "p2rank_score": "" if self.p2rank_score is None else self.p2rank_score,
            "center_x": self.center[0],
            "center_y": self.center[1],
            "center_z": self.center[2],
            "surf_atom_count": "" if self.surf_atom_count is None else self.surf_atom_count,
            "residue_ids": ";".join(self.residue_ids),
        }

    @classmethod
    def from_dict(cls, row: dict[str, str]) -> "PocketRow":
        def _float(key: str) -> float | None:
            raw = (row.get(key) or "").strip()
            if not raw:
                return None
            try:
                return float(raw)
            except ValueError:
                return None

        residues = (row.get("residue_ids") or "").strip()
        return cls(
            pocket_id=row["pocket_id"],
            rank=int(float(row["rank"])),
            p2rank_score=_float("p2rank_score"),
            center=(
                float(row["center_x"]),
                float(row["center_y"]),
                float(row["center_z"]),
            ),
            surf_atom_count=int(float(row["surf_atom_count"])) if (row.get("surf_atom_count") or "").strip() else None,
            residue_ids=tuple(part for part in residues.split(";") if part),
        )


def write_poses_csv(path: Path | str, rows: list[PoseRow]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.to_dict())
    return path


def write_pockets_csv(path: Path | str, rows: list[PocketRow]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POCKET_CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.to_dict())
    return path


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


def load_pockets_csv(path: Path | str) -> list[PocketRow]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as f:
        return [PocketRow.from_dict(row) for row in csv.DictReader(f)]


def best_vina_per_ligand_pocket(rows: list[PoseRow]) -> dict[tuple[str, str], PoseRow]:
    best: dict[tuple[str, str], PoseRow] = {}
    for row in rows:
        key = (row.ligand_id, row.pocket_id)
        if key not in best or row.vina_score < best[key].vina_score:
            best[key] = row
    return best


def best_pose_overall(rows: list[PoseRow]) -> PoseRow | None:
    """Поза №1 по Vina: самая отрицательная оценка среди всех карманов прогона.

    Сравнение корректно только потому, что все позы прогона посчитаны одной
    функцией скоринга с одними настройками. Между разными прогонами так
    сравнивать нельзя.
    """
    if not rows:
        return None
    return min(rows, key=lambda r: (r.vina_score, r.pocket_rank, r.pose_id))
