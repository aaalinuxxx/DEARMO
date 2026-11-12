"""Тесты ранжирования кандидатов, сводной таблицы сравнения, экспорта CSV/JSON и CLI."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import cli
from models.request import ResearchRequest
from models.result import Evidence, ResearchResult
from orchestration.state import WorkflowState
from tools import bindingdb_client, pdb_client
from tools.export import export_candidates_csv, export_run_bundle, export_to_csv, export_to_json
from tools.ranking import build_candidates_comparison_table, rank_candidates


def test_rank_candidates_orders_passed_first_and_records_failed_reasons():
    candidates = [
        {
            "molecule": "Heavy_Rejected",
            "smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
            "ki_nm": 1.0,
            "unit": "nM",
        },
        {
            "molecule": "Potent_Passed",
            "smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
            "ki_nm": 8.0,
            "unit": "nM",
        },
        {
            "molecule": "Weak_Passed",
            "smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
            "ki_nm": 45.0,
            "unit": "nM",
        },
    ]
    # Первому кандидату искусственно зададим несоответствие порогу ki_max_nm=50, а проверим отсев по прямому дескриптору
    candidates[0] = {
        "molecule": "Heavy_Rejected",
        "mw": 620.0,
        "logp": 6.1,
        "tpsa": 155.0,
        "hbd": 6,
        "hba": 12,
        "rotatable_bonds": 12,
        "qed": 0.12,
        "ki_nm": 1.0,
        "unit": "nM",
    }
    ranked = rank_candidates(candidates)
    assert ranked[0]["molecule"] == "Potent_Passed"
    assert ranked[0]["rank"] == 1
    assert ranked[0]["filter_status"] == "passed"

    assert ranked[1]["molecule"] == "Weak_Passed"
    assert ranked[1]["rank"] == 2
    assert ranked[1]["filter_status"] == "passed"

    assert ranked[2]["molecule"] == "Heavy_Rejected"
    assert ranked[2]["rank"] == 3
    assert ranked[2]["filter_status"] == "rejected"
    assert len(ranked[2]["failed_reasons"]) >= 3


def test_rank_candidates_handles_missing_and_incomparable_activity():
    candidates = [
        {
            "molecule": "No_Activity_Data",
            "smiles": "CCO",
            "ki_nm": None,
        },
        {
            "molecule": "Incomparable_Percent",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "ki": 85.0,
            "unit": "%",
        },
        {
            "molecule": "Valid_Ki_uM",
            "smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
            "ki": 0.012,
            "unit": "uM",
        },
    ]
    ranked = rank_candidates(candidates)
    assert ranked[0]["molecule"] == "Valid_Ki_uM"
    assert ranked[0]["ki_nm"] == pytest.approx(12.0)
    assert ranked[0]["preliminary_choice"] is False

    by_name = {r["molecule"]: r for r in ranked}
    assert by_name["No_Activity_Data"]["ki_nm"] is None
    assert "experimental_Ki" in by_name["No_Activity_Data"]["missing_data"]
    assert by_name["No_Activity_Data"]["preliminary_choice"] is True

    assert by_name["Incomparable_Percent"]["ki_nm"] is None
    assert by_name["Incomparable_Percent"]["comparable"] is False
    assert by_name["Incomparable_Percent"]["comparability_warnings"]


def test_build_candidates_comparison_table_end_to_end():
    table = build_candidates_comparison_table(
        ["acetazolamide", "methazolamide", "ethoxzolamide"],
        target_id="CHEMBL205",
    )
    assert len(table) == 3
    assert all(r["filter_status"] == "passed" for r in table)
    # Ethoxzolamide имеет Ki = 8.0 nM, Acetazolamide = 12.0 nM, Methazolamide = 14.0 nM
    ki_vals = [r["ki_nm"] for r in table if r["ki_nm"] is not None]
    assert len(ki_vals) >= 1 and ki_vals == sorted(ki_vals)


def test_export_csv_and_json_with_cache_metadata(tmp_path: Path):
    req = ResearchRequest(text="Подбери белок для ацетазоламида")
    state = WorkflowState(request=req, scenario_id="target_for_ligand")
    ev = Evidence(
        kind="database",
        value=12.0,
        unit="nM",
        source="ChEMBL",
        method="chembl_ki",
        context={"is_cached": True, "fetched_at": "2026-11-12T00:00:00Z"},
    )
    state.add_evidence(ev)
    state.observations["candidates_table"] = build_candidates_comparison_table(["acetazolamide", "methazolamide"])
    result = ResearchResult(
        status="success",
        evidence=[ev],
        explanation="Ki ацетазоламида составляет 12.0 нМ.",
        limitations=["Тестовое ограничение"],
    )

    artifacts = export_run_bundle(result, state, export_dir=tmp_path, prefix="test_run")
    assert Path(artifacts["evidence_csv"]).is_file()
    assert Path(artifacts["candidates_csv"]).is_file()
    assert Path(artifacts["report_json"]).is_file()

    with Path(artifacts["evidence_csv"]).open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["is_cached"] == "True"
    assert rows[0]["fetched_at"] == "2026-11-12T00:00:00Z"

    with Path(artifacts["report_json"]).open(encoding="utf-8") as f:
        report = json.load(f)
    assert report["status"] == "success"
    assert report["evidence"][0]["is_cached"] is True
    assert report["evidence"][0]["fetched_at"] == "2026-11-12T00:00:00Z"
    assert len(report["candidates_table"]) == 2


def test_bindingdb_and_pdb_clients_expose_cache_metadata():
    bdb = bindingdb_client.fetch_activities("ацетазоламид")
    assert "is_cached" in bdb
    assert "fetched_at" in bdb and bdb["fetched_at"]

    pdb_meta = pdb_client.get_structure_metadata("3HS4")
    assert pdb_meta["pdb_id"] == "3HS4"
    assert "is_cached" in pdb_meta
    assert "fetched_at" in pdb_meta and pdb_meta["fetched_at"]


def test_cli_ask_with_export_and_interactive_clarification(tmp_path: Path):
    # 1. Обычный запрос с экспортом
    code = cli.main(["ask", "Подбери белок для докинга ацетазоламида", "--export-dir", str(tmp_path)])
    assert code == 0
    assert any(tmp_path.glob("*_evidence.csv"))
    assert any(tmp_path.glob("*_report.json"))

    # 2. Интерактивное доуточнение при неоднозначном запросе
    res, _ = cli._print_result(
        "Подбери белок для докинга сульфонамида",
        interactive=True,
        input_fn=lambda _prompt: "1",  # выбираем 1-го кандидата (Acetazolamide)
    )
    assert res.status == "success"
