"""ORC: карманы обрабатываются параллельными воркерами, чужой прогон не подставляется."""
import pytest

from agents import docking_agent
from models.request import ResearchRequest
from orchestration.workflow import run_workflow


def test_saved_run_is_used_only_for_the_same_pair(monkeypatch):
    """Живой докинг выключен — берётся сохранённый прогон ровно этой пары, с двумя карманами."""
    monkeypatch.setenv("DOCKING_LIVE", "0")
    result, state = run_workflow(
        ResearchRequest(
            text="Выполни докинг acetazolamide в структуру 3HS4",
            molecules=["acetazolamide"],
            ligand_ids=["acetazolamide"],
            pdb_id="3HS4",
        )
    )
    assert result.status == "success"
    assert state.observations.get("docking_mode") == "saved"
    assert [item.pocket_id for item in state.work_queue] == ["pocket1", "pocket2"]
    assert all(item.status == "done" for item in state.work_queue)
    assert state.artifacts.get("report_md")
    assert any(row["agent"] == "orc" for row in state.agent_trace)


def test_unknown_pair_does_not_borrow_another_molecules_poses(monkeypatch, tmp_path):
    """Для пары без своего прогона докинг честно не выполняется, а не берёт чужие позы."""
    monkeypatch.setenv("DOCKING_LIVE", "0")
    # Каталоги сохранённых прогонов пустые: в них заведомо нет ничего подходящего.
    monkeypatch.setattr(docking_agent, "_EXAMPLES", tmp_path / "examples")
    monkeypatch.setattr(docking_agent, "_ROOT", tmp_path)
    result, state = run_workflow(
        ResearchRequest(
            text="Выполни докинг erlotinib в структуру 1M17",
            molecules=["erlotinib"],
            ligand_ids=["erlotinib"],
            pdb_id="1M17",
        )
    )
    assert result.status != "success"
    vina = [e for e in result.evidence if e.method == "vina_score"]
    assert not vina, "оценка Vina не должна появляться без расчёта по этой паре"


@pytest.mark.parametrize("ligand,pdb", [("erlotinib", "1M17")])
def test_live_docking_scores_a_fresh_pair(ligand, pdb, monkeypatch):
    """Живой Vina считает пару, для которой в репозитории нет сохранённого прогона."""
    monkeypatch.setenv("DOCKING_LIVE", "1")
    result, state = run_workflow(
        ResearchRequest(
            text=f"Выполни докинг {ligand} в структуру {pdb}",
            molecules=[ligand],
            ligand_ids=[ligand],
            pdb_id=pdb,
        )
    )
    assert state.observations.get("docking_mode") == "live"
    scores = [float(e.value) for e in result.evidence if e.method == "vina_score"]
    assert scores and all(-20.0 < s < 0.0 for s in scores)
