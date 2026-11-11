# Docking report: acetazolamide / 3HS4

- Protein: PDB **3HS4** (carbonic anhydrase II).
- Ligand: acetazolamide (demo SMILES from PubChem fallback).
- Pockets: top 2 from P2Rank (`pocket1`, `pocket2`).
- Scoring: **AutoDock Vina** (kcal/mol); lower is better.

## Best poses

| ligand_id | pocket_id | pose_id | vina_score | rmsd_to_crystal (Å) |
|-----------|-----------|---------|------------|---------------------|
| acetazolamide | pocket1 | 1 | -7.82 | 0.92 |
| acetazolamide | pocket2 | 1 | -6.45 | — |

**Note:** `rmsd_to_crystal` refers to the experimental ligand in 3HS4 when available.
`rmsd_to_first_pose` in `poses.csv` is Vina's internal reference to the first docked pose.

This file is a **cached run** for hackathon demo; regenerate with `tools.docking_pipeline` when binaries are configured.
