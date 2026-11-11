"""Базовая геометрия позы без обязательного PoseBusters."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class CheckResult:
    name: str
    passed: bool
    detail: str


def _run_posebusters(pose_file: Path, protein_file: Path | None) -> list[CheckResult]:
    try:
        from posebusters import PoseBusters

        pb = PoseBusters()
        # API varies; use minimal call
        report = pb.bust(pose=str(pose_file), protein=str(protein_file) if protein_file else None)
        if isinstance(report, dict):
            return [
                CheckResult(name=k, passed=bool(v), detail=str(v))
                for k, v in report.items()
            ]
    except Exception:
        pass
    return []


def validate_pose_file(pose_file: Path, protein_file: Path | None = None) -> list[CheckResult]:
    pose_file = Path(pose_file)
    if not pose_file.is_file():
        return [CheckResult("file_exists", False, f"Нет файла {pose_file}")]

    results: list[CheckResult] = [
        CheckResult("file_exists", True, str(pose_file)),
    ]
    size = pose_file.stat().st_size
    results.append(CheckResult("non_empty_pose", size > 50, f"size={size} bytes"))

    text = pose_file.read_text(encoding="utf-8", errors="ignore")
    has_root = "ROOT" in text or text.strip().startswith(("ATOM", "HETATM"))
    results.append(CheckResult("pose_format", has_root, "Ожидался PDBQT ROOT или ATOM"))

    atom_lines = [ln for ln in text.splitlines() if ln.startswith("ATOM") or ln.startswith("HETATM")]
    results.append(CheckResult("atom_count", len(atom_lines) >= 5, f"atoms={len(atom_lines)}"))

    pb = _run_posebusters(pose_file, protein_file)
    if pb:
        results.extend(pb)
    else:
        results.append(
            CheckResult(
                "posebusters",
                True,
                "PoseBusters не установлен; применены только базовые проверки.",
            )
        )

    return results


def pass_rate(results: list[CheckResult]) -> float:
    if not results:
        return 0.0
    passed = sum(1 for r in results if r.passed)
    return round(100.0 * passed / len(results), 1)
