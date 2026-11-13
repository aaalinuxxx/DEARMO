"""Редокинг-бенчмарк: докинг вслепую в структуры с известным кристаллическим лигандом.

Единственная метрика, по которой можно честно сказать «докинг работает», — доля
систем, где поза №1 лежит не дальше 2 Å от кристаллической. Здесь каждая система
проходит полный путь сценария 2: центры боксов выбирает P2Rank по структуре без
лиганда, координаты кристаллического лиганда в расчёт не подаются и служат только
контролем. Провалы не скрываются и попадают в таблицу наравне с успехами.

    python -m tools.redock_benchmark              # все системы, пишет docs/redocking.md
    python -m tools.redock_benchmark 3PTB:benzamidine 3HS4:acetazolamide
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Бенчмарк проверяет расчёт, а не языковую модель: ответ LLM здесь не участвует.
os.environ.setdefault("DEARMO_LLM_PROVIDER", "keyword")

from models.request import ResearchRequest  # noqa: E402
from orchestration.workflow import run_workflow  # noqa: E402
from tools.live_docking import RMSD_SUCCESS_A  # noqa: E402

_ROOT = Path(__file__).resolve().parent.parent

# (PDB, название лиганда). Это перечень задач, а не химические знания: SMILES
# лиганда система получает из PubChem по названию, код кристаллического
# компонента — из RCSB по структуре молекулы.
DEFAULT_SYSTEMS: tuple[tuple[str, str], ...] = (
    ("3PTB", "benzamidine"),
    ("3HS4", "acetazolamide"),
    ("1M17", "erlotinib"),
    ("4WKQ", "gefitinib"),
    ("2OIQ", "imatinib"),
)

_REQUEST = (
    "Найди три наиболее перспективных кармана, проведи докинг в каждом "
    "и покажи лучшую позу по Vina."
)


def _slug(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name.strip())


def _min_rmsd(poses_csv: Path) -> float | None:
    if not poses_csv.is_file():
        return None
    values: list[float] = []
    with poses_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            raw = (row.get("rmsd_to_crystal") or "").strip()
            if raw:
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
    return min(values) if values else None


def run_system(pdb_id: str, ligand: str) -> dict:
    started = time.perf_counter()
    result, state = run_workflow(
        ResearchRequest(text=f"{_REQUEST} Лиганд: {ligand}, структура {pdb_id}.", molecules=[ligand], pdb_id=pdb_id)
    )
    row: dict = {
        "pdb_id": pdb_id,
        "ligand": ligand,
        "status": result.status,
        "seconds": round(time.perf_counter() - started),
        "errors": list(state.errors)[:2],
    }
    docking_dir = _ROOT / "outputs" / "docking"
    meta_file = None
    for candidate in sorted(docking_dir.glob(f"{pdb_id.upper()}_*/run_meta.json")):
        meta = json.loads(candidate.read_text(encoding="utf-8"))
        if str(meta.get("ligand_id", "")).lower().startswith(_slug(ligand).lower()[:6]):
            meta_file = candidate
            break
    if meta_file is None:
        return row

    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    control = meta.get("control") or {}
    best = meta.get("best_pose") or {}
    row.update(
        {
            "receptor_prep": meta.get("receptor_prep"),
            "pockets": meta.get("pocket_count"),
            "best_vina": best.get("vina_score"),
            "crystal_comp": control.get("crystal_comp_id"),
            "pocket_rank": control.get("pocket_hit_rank"),
            "pocket_distance": control.get("pocket_hit_distance_a"),
            "rmsd_pose1": control.get("best_pose_rmsd_to_crystal_a"),
            "rmsd_min": _min_rmsd(meta_file.parent / "poses.csv"),
            "control_available": bool(control.get("available")),
            "reason": control.get("reason"),
        }
    )
    rmsd = row.get("rmsd_pose1")
    row["success"] = (rmsd is not None) and rmsd <= RMSD_SUCCESS_A
    return row


def _fmt(value, spec: str = ".2f") -> str:
    return "—" if value is None else format(value, spec)


def write_report(rows: list[dict], path: Path) -> Path:
    measured = [r for r in rows if r.get("control_available")]
    passed = [r for r in measured if r.get("success")]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# Редокинг-бенчмарк",
        "",
        f"Сгенерировано командой `python -m tools.redock_benchmark` {stamp}.",
        "",
        "Докинг вслепую: центры боксов даёт P2Rank по структуре без лиганда, координаты "
        "кристаллического лиганда в расчёт не подаются и служат только контролем. "
        f"Успех — поза №1 в пределах {RMSD_SUCCESS_A:g} Å от кристаллической.",
        "",
        f"**Итог: {len(passed)} из {len(measured)} систем** "
        + (f"({100 * len(passed) / len(measured):.0f} %)." if measured else "(контроль не посчитан)."),
        "",
        "| Структура | Лиганд | Подготовка рецептора | Карман кристалла (ранг / расстояние, Å) "
        "| Vina поза №1, ккал/моль | RMSD поза №1, Å | Мин. RMSD среди поз, Å | ≤ порога |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        if not r.get("control_available"):
            why = r.get("reason") or "; ".join(r.get("errors") or []) or "нет данных"
            lines.append(f"| {r['pdb_id']} | {r['ligand']} | — | — | — | — | — | контроль не посчитан: {why} |")
            continue
        lines.append(
            f"| {r['pdb_id']} | {r['ligand']} | {r.get('receptor_prep') or '—'} | "
            f"{r.get('pocket_rank') or '—'} / {_fmt(r.get('pocket_distance'))} | "
            f"{_fmt(r.get('best_vina'), '.2f')} | {_fmt(r.get('rmsd_pose1'))} | "
            f"{_fmt(r.get('rmsd_min'))} | {'да' if r.get('success') else 'нет'} |"
        )
    lines += [
        "",
        "Оценка Vina — расчётная величина, а не экспериментальная константа. Сравнивать её "
        "между системами нельзя: она зависит от размера лиганда и белка.",
        "",
        "Провал по RMSD при верно найденном кармане означает, что скоринговая функция Vina "
        "ставит поза №1 не в нативную ориентацию. Это известное ограничение метода, а не "
        "ошибка подготовки: геометрия и типы атомов проверяются отдельно.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    systems = (
        tuple(tuple(item.split(":", 1)) for item in args if ":" in item)  # type: ignore[arg-type]
        or DEFAULT_SYSTEMS
    )
    rows: list[dict] = []
    for pdb_id, ligand in systems:
        print(f"→ {pdb_id} / {ligand} …", flush=True)
        row = run_system(pdb_id, ligand)
        rows.append(row)
        print(
            f"   {row['status']}, {row['seconds']} с, RMSD поза №1 = {_fmt(row.get('rmsd_pose1'))} Å, "
            f"успех: {row.get('success')}",
            flush=True,
        )
    path = write_report(rows, _ROOT / "docs" / "redocking.md")
    measured = [r for r in rows if r.get("control_available")]
    passed = [r for r in measured if r.get("success")]
    print(f"\nУспешно {len(passed)} из {len(measured)}. Отчёт: {path.relative_to(_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
