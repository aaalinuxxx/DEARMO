"""Обёртка P2Rank (subprocess)."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def p2rank_binary() -> str | None:
    return os.environ.get("P2RANK_BIN") or shutil.which("prank")


def run_p2rank(protein_pdb: Path, output_dir: Path) -> Path:
    binary = p2rank_binary()
    if not binary:
        raise RuntimeError("P2Rank не найден. Задайте P2RANK_BIN или установите prank.")
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "predict", str(protein_pdb), "-o", str(output_dir)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    csv_candidates = list(output_dir.glob("**/predictions.csv"))
    if not csv_candidates:
        raise RuntimeError("P2Rank не создал predictions.csv")
    return csv_candidates[0]
