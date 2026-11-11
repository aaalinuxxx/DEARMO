"""Агент анализа: сравнение карманов и графики."""
from __future__ import annotations

from pathlib import Path

from models.docking import best_vina_per_ligand_pocket, load_poses_csv
from models.result import Evidence
from orchestration.state import WorkflowState
from tools.plotting import pocket_ligand_matrix, write_heatmap

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUT = _REPO_ROOT / "outputs" / "analysis"


def compare_pockets(state: WorkflowState, csv_path: Path) -> WorkflowState:
    with state.agent_step("analysis_agent", "compare_pockets"):
        rows = load_poses_csv(csv_path)
        best = best_vina_per_ligand_pocket(rows)
        ligands, pockets, matrix = pocket_ligand_matrix(rows)

        for (lig, pocket), row in sorted(best.items()):
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=row.vina_score,
                    unit="kcal/mol",
                    method=f"best_vina_{lig}_{pocket}",
                    source="poses.csv",
                    context={"ligand_id": lig, "pocket_id": pocket, "pose_file": row.pose_file},
                )
            )

        if len(ligands) >= 2 and len(pockets) >= 2:
            ranks_by_pocket: dict[str, list[str]] = {}
            for pocket in pockets:
                cols = [(lig, best.get((lig, pocket))) for lig in ligands]
                cols = [(lig, r) for lig, r in cols if r is not None]
                cols.sort(key=lambda x: x[1].vina_score)
                ranks_by_pocket[pocket] = [lig for lig, _ in cols]
            pocket_ids = list(ranks_by_pocket.keys())
            if len(pocket_ids) >= 2 and ranks_by_pocket[pocket_ids[0]] != ranks_by_pocket[pocket_ids[1]]:
                state.add_evidence(
                    Evidence(
                        kind="calculation",
                        value="order_changes",
                        method="ligand_ranking_by_pocket",
                        source="analysis",
                        context={"ranks": ranks_by_pocket},
                    )
                )

        _OUT.mkdir(parents=True, exist_ok=True)
        heatmap_path = write_heatmap(rows, _OUT / "vina_heatmap.png")
        state.artifacts["heatmap_png"] = str(heatmap_path)
        state.limitations.append("Тепловая карта строится по лучшему Vina score на пару лиганд×карман.")
    return state
