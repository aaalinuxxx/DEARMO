"""Маршрутизация workflow на ключевых словах, без вызова Qwen.

Молекулы здесь намеренно не те, что лежат в предзагруженных данных: маршрут
не должен зависеть от конкретного соединения.
"""
from models.request import ResearchRequest
from orchestration.scenarios import classify_by_keywords
from orchestration.workflow import run_workflow


def test_keyword_classifier_targets():
    assert classify_by_keywords("Подбери белок для докинга этого соединения") == "target_for_ligand"


def test_scenario1_target_for_unseen_molecule():
    """Мишень и структура подбираются для молекулы, которой нет ни в одном словаре кода."""
    result, state = run_workflow(
        ResearchRequest(
            text="Подбери белок для докинга erlotinib. Обоснуй выбор и предложи структуру из PDB.",
            molecules=["erlotinib"],
        )
    )
    assert state.scenario_id == "target_for_ligand"
    assert result.status in ("success", "no_data", "partial")
    if result.status == "success":
        cids = [e.value for e in result.evidence if e.method == "pubchem_cid"]
        # CID приходит из живого PubChem, а не из подстановки ацетазоламида (1986).
        assert cids and 1986 not in cids


def test_compound_request_runs_target_then_dock_then_pose():
    """Составной запрос проходит три шага подряд с живым докингом."""
    result, state = run_workflow(
        ResearchRequest(
            text="Подбери белок для erlotinib, затем выполни докинг в 1M17 и проверь правдоподобность позы.",
            molecules=["erlotinib"],
            pdb_id="1M17",
        )
    )
    assert result.status == "success"
    assert state.observations["completed_steps"] == [
        "target_for_ligand",
        "dock_pockets",
        "validate_pose",
    ]
    scores = [e.value for e in result.evidence if e.method == "vina_score"]
    assert scores, "докинг обязан дать оценку Vina"
    assert all(-20.0 < float(s) < 0.0 for s in scores), "оценка Vina должна быть физичной"
    assert state.observations.get("docking_mode") == "live"


def test_validate_pose_builds_its_own_pose_instead_of_borrowing():
    """Запрос без позы и без структуры: система достраивает докинг сама, чужую позу не берёт."""
    result, state = run_workflow(
        ResearchRequest(
            text="Проверь правдоподобность позы соединения erlotinib.",
            molecules=["erlotinib"],
        )
    )
    checks = [e for e in result.evidence if e.method == "checks_passed_percent"]
    if not checks:
        # Допустимо не справиться (нет сети), но тогда и результата проверки нет.
        assert result.status != "success"
        return

    # Проверка позы состоялась — значит, поза была посчитана здесь же.
    assert any(e.method == "vina_score" for e in result.evidence)
    assert state.observations.get("docking_mode") == "live"
    added = {step["added"] for step in state.observations.get("auto_steps") or []}
    assert "dock_pockets" in added, "докинг должен быть добавлен автоматически"


def test_dependency_resolution_orders_steps_by_execution():
    """Недостающий вход достраивается производящим шагом, и он идёт в плане первым."""
    _result, state = run_workflow(
        ResearchRequest(
            text="Сравни лиганды по карманам, меняется ли порядок по Vina.",
            molecules=["erlotinib"],
            pdb_id="1M17",
        )
    )
    assert state.plan[:2] == ["ok:dock_pockets", "ok:compare_docking"]
    missing = {step["missing"] for step in state.observations.get("auto_steps") or []}
    assert "poses_csv_path" in missing
