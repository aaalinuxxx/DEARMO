"""CLI докинга: опциональный live P2Rank+Vina; иначе справка."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="DEARMO docking pipeline (P2Rank + Vina)")
    parser.add_argument("--pdb", required=True, help="PDB ID или путь к PDB")
    parser.add_argument("--smiles", required=True, help="SMILES лиганда")
    parser.add_argument("--output", type=Path, default=Path("outputs/docking_run"))
    parser.add_argument("--live", action="store_true", help="Запустить P2Rank/Vina (нужны бинарники)")
    args = parser.parse_args(argv)

    if not args.live:
        print(
            "Режим без --live: сохраните прогон в data/examples/ вручную или используйте кеш.",
            file=sys.stderr,
        )
        return 2

    from tools.p2rank_runner import run_p2rank
    from tools.vina_runner import run_vina

    pdb_path = Path(args.pdb)
    if len(args.pdb) == 4:
        cache = Path(__file__).resolve().parent.parent / "data" / "cache" / "pdb" / f"{args.pdb.upper()}.pdb"
        if cache.is_file():
            pdb_path = cache
    run_p2rank(pdb_path, args.output)
    print("P2Rank завершён; подготовьте receptor/ligand PDBQT и вызовите run_vina.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
