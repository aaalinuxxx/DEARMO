"""AutoDock Vina через subprocess: запуск, разбор всех поз, нарезка на файлы.

Vina печатает таблицу режимов и пишет все позы в один многомодельный PDBQT.
Раньше отсюда забиралась только первая строка таблицы, поэтому в системе
физически не существовало ни остальных поз, ни оценок RMSD между ними —
а сценарий 2 обязан показать и то и другое.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VinaMode:
    """Одна строка таблицы режимов Vina.

    rmsd_lb и rmsd_ub — расстояние до ПЕРВОЙ позы того же запуска, как его
    считает сам Vina (столбцы «rmsd l.b.» и «rmsd u.b.»). Для первой позы оба
    равны нулю по построению, и это не измерение, а точка отсчёта.
    """

    mode: int
    affinity: float
    rmsd_lb: float
    rmsd_ub: float


# Строка таблицы: номер режима, аффинность, две оценки RMSD. Заголовок и
# разделители под шаблон не подходят, поэтому отдельно их отбрасывать не нужно.
_MODE_LINE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$"
)


def vina_binary() -> str | None:
    configured = os.environ.get("VINA_BIN")
    if configured and Path(configured).is_file():
        return configured
    bundled = Path(__file__).resolve().parent / "bin" / "vina"
    if bundled.is_file():
        return str(bundled)
    return shutil.which("vina")


def parse_vina_modes(stdout: str) -> list[VinaMode]:
    """Все режимы из таблицы Vina, в порядке печати (он же порядок MODEL в --out)."""
    modes: list[VinaMode] = []
    for line in stdout.splitlines():
        match = _MODE_LINE.match(line)
        if not match:
            continue
        mode, affinity, lb, ub = match.groups()
        modes.append(
            VinaMode(
                mode=int(mode),
                affinity=float(affinity),
                rmsd_lb=float(lb),
                rmsd_ub=float(ub),
            )
        )
    return modes


def split_pdbqt_models(multi_model: Path, out_dir: Path, prefix: str) -> list[Path]:
    """Разложить многомодельный PDBQT по одному файлу на позу.

    Сценарий 3 и проверка геометрии работают с одной позой за раз, и разбор
    мультимодели в каждом читателе — лишний источник расхождений. Нумерация
    файлов идёт по порядку MODEL, то есть совпадает с номером режима Vina.
    """
    text = Path(multi_model).read_text(encoding="utf-8", errors="replace")
    out_dir.mkdir(parents=True, exist_ok=True)

    blocks: list[list[str]] = []
    current: list[str] | None = None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = []
            continue
        if line.startswith("ENDMDL"):
            if current is not None:
                blocks.append(current)
            current = None
            continue
        if current is not None:
            current.append(line)

    # Файл без MODEL/ENDMDL — это одна поза, записанная без обёртки.
    if not blocks:
        stripped = [ln for ln in text.splitlines() if ln.strip()]
        blocks = [stripped] if stripped else []

    paths: list[Path] = []
    for index, block in enumerate(blocks, start=1):
        path = out_dir / f"{prefix}_pose{index}.pdbqt"
        path.write_text("\n".join(block).rstrip() + "\n", encoding="utf-8")
        paths.append(path)
    return paths


def run_vina_modes(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    out_pose: Path | None = None,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int | None = 42,
) -> tuple[list[VinaMode], Path]:
    """Запустить Vina и вернуть ВСЕ позы вместе с многомодельным файлом.

    seed фиксирован по умолчанию: без него два прогона одной пары дают разные
    числа, и воспроизвести показанный на защите результат нельзя.
    """
    binary = vina_binary()
    if not binary:
        raise RuntimeError("Vina не найден. Задайте VINA_BIN или положите бинарник в tools/bin/vina.")
    out_pose = Path(out_pose) if out_pose else ligand_pdbqt.with_suffix(".out.pdbqt")
    out_pose.parent.mkdir(parents=True, exist_ok=True)
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        binary,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--out", str(out_pose),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    modes = parse_vina_modes(proc.stdout)
    if not modes:
        raise RuntimeError(
            "Vina не вернул ни одной позы: таблица режимов пуста.\n" + proc.stdout[-2000:]
        )
    return modes, out_pose


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    out_pose: Path | None = None,
    exhaustiveness: int = 8,
) -> tuple[float, Path]:
    """Лучшая оценка и файл поз. Обёртка над run_vina_modes для старых вызовов."""
    modes, pose_path = run_vina_modes(
        receptor_pdbqt,
        ligand_pdbqt,
        center,
        size=size,
        out_pose=out_pose,
        exhaustiveness=exhaustiveness,
    )
    return modes[0].affinity, pose_path
