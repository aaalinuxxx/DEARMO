"""Обёртка P2Rank: предсказание карманов и разбор predictions.csv.

Карманы — единственный источник центров боксов для докинга «вслепую». Центр по
кристаллическому лиганду здесь намеренно не используется: иначе система знала бы
ответ заранее, и ранг найденного кармана ничего бы не доказывал.
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

from models.docking import PocketRow

_ROOT = Path(__file__).resolve().parent.parent
_BUNDLED_P2RANK = _ROOT / "tools" / "bin" / "p2rank" / "prank"
_BUNDLED_JDK = _ROOT / "tools" / "bin" / "jdk21"


def p2rank_binary() -> str | None:
    configured = os.environ.get("P2RANK_BIN")
    if configured and Path(configured).is_file():
        return configured
    if _BUNDLED_P2RANK.is_file():
        return str(_BUNDLED_P2RANK)
    return shutil.which("prank")


def bundled_java_home() -> str | None:
    """JDK рядом с проектом, если он положен.

    P2Rank 2.5 собран на Groovy 4, который не читает классы новых JDK и падает
    на «Unsupported class file major version». Поэтому системная java не годится
    по умолчанию, и рядом кладётся отдельный JDK 21.
    """
    if not _BUNDLED_JDK.is_dir():
        return None
    for candidate in (_BUNDLED_JDK / "Contents" / "Home", _BUNDLED_JDK):
        if (candidate / "bin" / "java").is_file():
            return str(candidate)
    return None


def _p2rank_env() -> dict[str, str]:
    env = dict(os.environ)
    java_home = bundled_java_home()
    # Явно заданный JAVA_HOME уважается: пользователь мог указать подходящий JDK сам.
    if java_home and not env.get("JAVA_HOME"):
        env["JAVA_HOME"] = java_home
    return env


def parse_predictions(csv_path: Path | str) -> list[PocketRow]:
    """Разобрать predictions.csv. P2Rank дополняет имена колонок пробелами."""
    rows: list[PocketRow] = []
    with Path(csv_path).open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = {(k or "").strip(): (v or "").strip() for k, v in raw.items() if k}
            name = row.get("name")
            if not name:
                continue

            def _float(key: str) -> float | None:
                value = row.get(key, "")
                try:
                    return float(value)
                except ValueError:
                    return None

            cx, cy, cz = _float("center_x"), _float("center_y"), _float("center_z")
            if cx is None or cy is None or cz is None:
                # Карман без центра бесполезен: бокс по нему не поставить.
                continue
            try:
                rank = int(float(row.get("rank", "")))
            except ValueError:
                rank = len(rows) + 1
            surf = row.get("surf_atoms", "")
            rows.append(
                PocketRow(
                    pocket_id=name,
                    rank=rank,
                    p2rank_score=_float("score"),
                    center=(cx, cy, cz),
                    surf_atom_count=int(float(surf)) if surf else None,
                    residue_ids=tuple(part for part in row.get("residue_ids", "").split() if part),
                )
            )
    rows.sort(key=lambda r: r.rank)
    return rows


def run_p2rank(protein_pdb: Path, output_dir: Path) -> Path:
    """Запустить предсказание и вернуть путь к predictions.csv."""
    binary = p2rank_binary()
    if not binary:
        raise RuntimeError(
            "P2Rank не найден. Задайте P2RANK_BIN или распакуйте релиз в tools/bin/p2rank."
        )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "predict", "-f", str(protein_pdb), "-o", str(output_dir)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=_p2rank_env())
    candidates = sorted(output_dir.glob("**/*predictions.csv"))
    if not candidates:
        tail = (proc.stdout or "")[-1500:] + (proc.stderr or "")[-1500:]
        raise RuntimeError(f"P2Rank не создал predictions.csv.\n{tail}")
    return candidates[0]


def predict_pockets(
    protein_pdb: Path,
    output_dir: Path,
    top_n: int = 3,
) -> list[PocketRow]:
    """Карманы-кандидаты по убыванию балла P2Rank, не более top_n.

    Сценарий 2 берёт до трёх карманов. Если P2Rank нашёл меньше, возвращается
    сколько есть: добивать список выдуманными центрами нельзя.
    """
    csv_path = run_p2rank(Path(protein_pdb), Path(output_dir))
    pockets = parse_predictions(csv_path)
    if not pockets:
        raise RuntimeError(f"P2Rank не нашёл ни одного кармана в {protein_pdb}")
    return pockets[: max(1, top_n)]
