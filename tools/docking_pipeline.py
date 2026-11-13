"""Отдельный CLI докинга одной пары: PDB + SMILES -> AutoDock Vina."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="DEARMO docking pipeline (P2Rank + Vina)")
    parser.add_argument("--pdb", required=True, help="PDB ID или путь к PDB")
    parser.add_argument("--smiles", required=True, help="SMILES лиганда")
    parser.add_argument("--output", type=Path, default=Path("outputs/docking_run"))
    parser.add_argument("--ligand-id", default="ligand", help="Имя лиганда для poses.csv")
    args = parser.parse_args(argv)

    from tools.live_docking import dock_pair

    try:
        result = dock_pair(args.pdb, args.smiles, args.ligand_id, args.output)
    except Exception as exc:
        print(f"Живой докинг не выполнен: {exc}", file=sys.stderr)
        return 1
    print(f"vina_score: {result.vina_score}")
    print(f"pose: {result.pose_file}")
    print(f"box: {result.box_source} center={result.center}")
    print(f"ligand_prep: {result.ligand_prep}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
