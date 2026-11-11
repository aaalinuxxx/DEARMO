"""Тонкая обёртка PoseBusters."""
from __future__ import annotations

from pathlib import Path

from tools.pose_validation import CheckResult, pass_rate, validate_pose_file

__all__ = ["validate_pose_file", "CheckResult", "pass_rate"]
