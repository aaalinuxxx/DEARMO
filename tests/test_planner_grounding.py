"""Извлечение планировщика не перебивает то, что известно, и не принимает заглушки."""
from __future__ import annotations

from agents import planner
from models.request import ResearchRequest
from orchestration.llm_client import PlannerExtraction
from orchestration.state import WorkflowState


def _ext(**kw) -> PlannerExtraction:
    base = dict(
        scenario_id="dock_pockets", ligands=[], smiles=None, target_hint=None, pdb_id=None,
        poses_csv_path=None, pose_file=None, ligand_id=None, explanation_hint="", source="ollama:test",
    )
    base.update(kw)
    return PlannerExtraction(**base)


def _apply(request: ResearchRequest, ext: PlannerExtraction) -> WorkflowState:
    return planner._apply_extraction(WorkflowState(request=request), ext)


def test_placeholder_ligand_id_from_prompt_example_is_not_used():
    """Qwen копировала «X» из примера: результаты докинга ложились в каталог 3HS4_x."""
    state = _apply(
        ResearchRequest(text="Найди три кармана и проведи докинг", molecules=["acetazolamide"], pdb_id="3HS4"),
        _ext(ligand_id="X", ligands=["X"]),
    )
    assert state.observations["ligand_id"] == "acetazolamide"


def test_ligand_id_follows_the_known_molecule_not_the_models_guess():
    state = _apply(
        ResearchRequest(text="Проведи докинг", molecules=["erlotinib"]),
        _ext(ligand_id="imatinib"),
    )
    assert state.observations["ligand_id"] == "erlotinib"


def test_ungrounded_ligand_id_without_any_known_molecule_is_dropped():
    state = _apply(ResearchRequest(text="Проведи докинг"), _ext(ligand_id="X"))
    assert "ligand_id" not in state.observations
    assert any("отброшен" in lim for lim in state.limitations)


def test_ligand_id_named_in_the_text_is_accepted():
    state = _apply(ResearchRequest(text="докинг imatinib в 2OIQ"), _ext(ligand_id="imatinib"))
    assert state.observations["ligand_id"] == "imatinib"


def test_invented_path_does_not_overwrite_the_chain_context(tmp_path):
    """«path/to/universal_properties.csv» затирал верный путь таблицы поз."""
    real = tmp_path / "poses.csv"
    real.write_text("ligand_id\n", encoding="utf-8")
    state = _apply(
        ResearchRequest(text="Сравни позы", poses_csv_path=str(real)),
        _ext(scenario_id="compare_docking", poses_csv_path="path/to/universal_properties.csv"),
    )
    assert state.request.poses_csv_path == str(real)
    assert state.observations["poses_csv_path"] == str(real)


def test_invented_path_is_dropped_and_reported():
    state = _apply(
        ResearchRequest(text="Сравни позы"),
        _ext(scenario_id="compare_docking", poses_csv_path="path/to/universal_properties.csv"),
    )
    assert not state.request.poses_csv_path
    assert "poses_csv_path" not in state.observations
    assert any("такого файла нет" in lim for lim in state.limitations)


def test_existing_path_from_the_model_is_accepted_when_user_gave_none(tmp_path):
    real = tmp_path / "poses.csv"
    real.write_text("ligand_id\n", encoding="utf-8")
    state = _apply(
        ResearchRequest(text="Сравни позы"),
        _ext(scenario_id="compare_docking", poses_csv_path=str(real)),
    )
    assert state.request.poses_csv_path == str(real)


def test_null_string_from_the_model_means_not_given():
    """Qwen пишет пустое значение строкой «null»: это не путь и не повод для предупреждения."""
    state = _apply(
        ResearchRequest(text="Сравни позы"),
        _ext(scenario_id="compare_docking", poses_csv_path="null", pose_file="None"),
    )
    assert not state.request.poses_csv_path
    assert not state.limitations, state.limitations


def test_path_from_the_chain_is_visible_to_the_executor(tmp_path):
    """Готовая таблица поз должна попасть в observations, даже если модель ничего не предложила.

    Иначе исполнитель считал её отсутствующей и запускал докинг заново, а отчёт
    писал, что система «сама добавила шаг».
    """
    real = tmp_path / "poses.csv"
    real.write_text("ligand_id\n", encoding="utf-8")
    state = _apply(
        ResearchRequest(text="Сравни позы", poses_csv_path=str(real)),
        _ext(scenario_id="compare_docking", poses_csv_path=None),
    )
    assert state.observations["poses_csv_path"] == str(real)
