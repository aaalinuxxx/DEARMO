"""AutoDock Vina через subprocess."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def vina_binary() -> str | None:
    return os.environ.get("VINA_BIN") or shutil.which("vina")


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    out_pose: Path | None = None,
    exhaustiveness: int = 8,
) -> tuple[float, Path]:
    binary = vina_binary()
    if not binary:
        raise RuntimeError("Vina не найден. Задайте VINA_BIN.")
    out_pose = out_pose or ligand_pdbqt.with_suffix(".out.pdbqt")
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        binary,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(cx),
        "--center_y",
        str(cy),
        "--center_z",
        str(cz),
        "--size_x",
        str(sx),
        "--size_y",
        str(sy),
        "--size_z",
        str(sz),
        "--out",
        str(out_pose),
        "--exhaustiveness",
        str(exhaustiveness),
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    score = None
    for line in proc.stdout.splitlines():
        if line.strip().startswith("1 "):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    score = float(parts[1])
                except ValueError:
                    pass
    if score is None:
        raise RuntimeError("Не удалось разобрать score Vina")
    return score, out_pose
