"""Живой докинг одной пары: PDB + SMILES → PDBQT → AutoDock Vina."""
from __future__ import annotations

import csv
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from models.docking import POSE_CSV_FIELDS
from tools.pdbqt_prep import crystal_ligand_center, ligand_pdbqt_from_smiles, protein_center, receptor_pdbqt, write_text
from tools.vina_runner import run_vina, vina_binary

_ROOT = Path(__file__).resolve().parent.parent


@dataclass
class LiveDockResult:
    vina_score: float
    pose_file: Path
    poses_csv: Path
    receptor_pdbqt: Path
    ligand_pdbqt: Path
    center: tuple[float, float, float]
    box_source: str
    ligand_prep: str
    pdb_path: Path


def locate_pdb(pdb_id: str) -> Path:
    clean = pdb_id.strip().upper()
    if len(clean) != 4:
        path = Path(pdb_id)
        if path.is_file():
            return path
        raise FileNotFoundError(f"Нет PDB-файла: {pdb_id}")
    for folder in (_ROOT / "data" / "external" / "rcsb", _ROOT / "data" / "cache" / "pdb"):
        candidate = folder / f"{clean}.pdb"
        if candidate.is_file():
            return candidate
    dest = _ROOT / "data" / "cache" / "pdb" / f"{clean}.pdb"
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{clean}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "DEARMO/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        dest.write_bytes(resp.read())
    return dest


def dock_pair(
    pdb_id: str,
    smiles: str,
    ligand_id: str,
    output_dir: Path,
    *,
    exhaustiveness: int = 8,
) -> LiveDockResult:
    if vina_binary() is None:
        raise RuntimeError("AutoDock Vina не найден. Положите бинарник в tools/bin/vina или задайте VINA_BIN.")
    pdb_path = locate_pdb(pdb_id)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    output_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = write_text(output_dir / "receptor.pdbqt", receptor_pdbqt(pdb_text))
    ligand_text, ligand_prep = ligand_pdbqt_from_smiles(smiles)
    ligand_path = write_text(output_dir / "ligand.pdbqt", ligand_text)
    center = crystal_ligand_center(pdb_text)
    if center is None:
        center = protein_center(pdb_text)
        box_source = "protein_center"
        size = (22.0, 22.0, 22.0)
    else:
        box_source = "crystal_ligand"
        size = (20.0, 20.0, 20.0)
    pose_path = output_dir / "pose_1.pdbqt"
    score, pose_path = run_vina(
        receptor_path,
        ligand_path,
        center,
        size=size,
        out_pose=pose_path,
        exhaustiveness=exhaustiveness,
    )
    csv_path = output_dir / "poses.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_CSV_FIELDS)
        writer.writeheader()
        writer.writerow(
            {
                "ligand_id": ligand_id,
                "pocket_id": "site1",
                "pocket_rank": 1,
                "p2rank_score": "",
                "pose_id": 1,
                "vina_score": score,
                "pose_file": str(pose_path),
                "rmsd_to_first_pose": 0.0,
                "rmsd_to_crystal": "",
            }
        )
    # Паспорт прогона: по нему и только по нему прогон потом сопоставляется с парой.
    (output_dir / "run_meta.json").write_text(
        json.dumps(
            {
                "pdb_id": str(pdb_id).upper(),
                "ligand_id": ligand_id,
                "smiles": smiles,
                "box_source": box_source,
                "ligand_prep": ligand_prep,
                "engine": "AutoDock Vina",
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return LiveDockResult(
        vina_score=score,
        pose_file=pose_path,
        poses_csv=csv_path,
        receptor_pdbqt=receptor_path,
        ligand_pdbqt=ligand_path,
        center=center,
        box_source=box_source,
        ligand_prep=ligand_prep,
        pdb_path=pdb_path,
    )
