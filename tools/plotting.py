"""Графики сравнения докинга по карманам."""
from __future__ import annotations

from pathlib import Path

from models.docking import PoseRow, best_vina_per_ligand_pocket

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None


def pocket_ligand_matrix(rows: list[PoseRow]) -> tuple[list[str], list[str], list[list[float]]]:
    best = best_vina_per_ligand_pocket(rows)
    ligands = sorted({k[0] for k in best})
    pockets = sorted({k[1] for k in best}, key=lambda p: (int(p.replace("pocket", "") or 0) if p.replace("pocket", "").isdigit() else 99, p))
    matrix = []
    for pocket in pockets:
        row_vals = []
        for lig in ligands:
            cell = best.get((lig, pocket))
            row_vals.append(cell.vina_score if cell else float("nan"))
        matrix.append(row_vals)
    return ligands, pockets, matrix


def write_heatmap(rows: list[PoseRow], output_path: Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None:
        output_path.write_text("matplotlib not installed", encoding="utf-8")
        return output_path
    ligands, pockets, matrix = pocket_ligand_matrix(rows)
    fig, ax = plt.subplots(figsize=(max(4, len(ligands)), max(3, len(pockets))))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(ligands)), labels=ligands, rotation=45, ha="right")
    ax.set_yticks(range(len(pockets)), labels=pockets)
    ax.set_title("Лучший Vina score по лиганд × карман")
    fig.colorbar(im, ax=ax, label="Vina score (kcal/mol)")
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if val == val:  # not nan
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=120)
    plt.close(fig)
    return output_path
