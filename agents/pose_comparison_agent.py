"""Сценарий 3: сравнение поз одного прогона докинга.

Что делает агент и чего он не делает.

Считает код: оценки Vina по позам и разницу с позой отсчёта, контакты белок–
лиганд по координатам, карту «поза × остаток», расстояния между одними и теми
же парами атомов во всех позах и закономерности в этих рядах.

Языковая модель не порождает ни одного числа. Ей достаётся готовый список
кандидатов, и вернуть она может только индекс в нём. Вывод «какую позу стоит
проверить подробнее» строится детерминированными правилами; модель вправе
выбрать другого кандидата из того же списка и переформулировать причину.
"""
from __future__ import annotations

import json
from pathlib import Path

from models.docking import PoseRow, load_poses_csv
from models.result import Evidence
from orchestration.llm_client import decide
from orchestration.state import WorkflowState
from tools import pose_metrics
from tools.contacts import Atom, read_structure
from tools.plotting import write_contact_map, write_distance_series, write_vina_by_pose
from tools.pose_metrics import (
    PoseAnalysis,
    PoseKey,
    collect_findings,
    contact_matrix,
    contact_rows,
    distance_series,
    key_atom_pairs,
    protein_state_rows,
    select_for_comparison,
    vina_deltas,
    write_csv,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Сколько пар атомов попадает на график расстояний. Больше — нечитаемо.
MAX_KEY_PAIRS = 6
# Сколько поз показывать на графиках. В таблицах остаются все.
MAX_PLOT_POSES = 12


def _resolve(path_value: str | None, run_dir: Path) -> Path | None:
    """Путь из таблицы: абсолютный, относительный каталогу прогона или корню репозитория."""
    if not path_value:
        return None
    candidate = Path(path_value)
    if candidate.is_file():
        return candidate
    for base in (run_dir, _REPO_ROOT):
        joined = base / path_value
        if joined.is_file():
            return joined
    name = candidate.name
    for folder in (run_dir, run_dir / "poses"):
        joined = folder / name
        if joined.is_file():
            return joined
    return None


def _load_run_meta(run_dir: Path) -> dict:
    meta_path = run_dir / "run_meta.json"
    if not meta_path.is_file():
        return {}
    try:
        data = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _load_receptor(run_dir: Path, meta: dict) -> tuple[list[Atom], str | None]:
    """Координаты белка. PDB предпочтительнее PDBQT: в PDBQT элемент — тип AutoDock."""
    candidates = [
        meta.get("receptor_pdb_path"),
        str(run_dir / "receptor.pdb"),
        meta.get("receptor_pdbqt"),
        str(run_dir / "receptor.pdbqt"),
    ]
    for value in candidates:
        path = _resolve(value, run_dir)
        if path is None:
            continue
        try:
            atoms = read_structure(path)
        except (OSError, ValueError):
            continue
        if atoms:
            return atoms, str(path)
    return [], None


def _reference_key(meta: dict, analyses: list[PoseAnalysis]) -> PoseKey | None:
    """Поза №1 по Vina. Берётся из паспорта прогона, а не вычисляется заново."""
    best = meta.get("best_pose")
    if isinstance(best, dict):
        pocket_id = best.get("pocket_id")
        pose_id = best.get("pose_id")
        if pocket_id is not None and pose_id is not None:
            for item in analyses:
                if item.key.pocket_id == str(pocket_id) and item.key.pose_id == int(pose_id):
                    return item.key
    return pose_metrics.best_pose_key(analyses)


def _build_analyses(
    rows: list[PoseRow],
    protein_atoms: list[Atom],
    run_dir: Path,
    state: WorkflowState,
) -> list[PoseAnalysis]:
    from tools.contacts import detect_contacts

    analyses: list[PoseAnalysis] = []
    missing: list[str] = []
    for row in rows:
        pose_path = _resolve(row.pose_file, run_dir)
        if pose_path is None:
            missing.append(row.pose_file)
            continue
        try:
            ligand_atoms = read_structure(pose_path)
        except (OSError, ValueError):
            missing.append(str(pose_path))
            continue
        if not ligand_atoms:
            missing.append(str(pose_path))
            continue
        analyses.append(
            PoseAnalysis(
                key=PoseKey.from_row(row),
                row=row,
                ligand_atoms=ligand_atoms,
                contacts=detect_contacts(ligand_atoms, protein_atoms),
            )
        )
    if missing:
        state.limitations.append(
            f"Файлы координат не найдены для {len(missing)} поз из {len(rows)}; "
            "эти позы в сравнение не вошли."
        )
    analyses.sort(key=lambda item: (item.key.pocket_id, item.key.pose_id))
    return analyses


# --------------------------------------------------------------------------
# Факты
# --------------------------------------------------------------------------


def _add_score_evidence(state: WorkflowState, deltas, analyses: list[PoseAnalysis]) -> None:
    by_key = {item.key: item for item in analyses}
    for delta in deltas:
        item = by_key.get(delta.key)
        context = {
            "ligand_id": delta.key.ligand_id,
            "pocket_id": delta.key.pocket_id,
            "pose_id": delta.key.pose_id,
            "pose_file": item.row.pose_file if item else None,
        }
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=delta.vina_score,
                unit="kcal/mol",
                method="vina_score_by_pose",
                source="poses.csv",
                context=context,
            )
        )
        if not delta.is_reference:
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=delta.delta_to_reference,
                    unit="kcal/mol",
                    method="vina_delta_to_best_pose",
                    source="poses.csv",
                    context={**context, "reference_pose": _reference_label(deltas)},
                )
            )


def _reference_label(deltas) -> str | None:
    for delta in deltas:
        if delta.is_reference:
            return delta.key.label
    return None


def _add_rmsd_evidence(state: WorkflowState, analyses: list[PoseAnalysis]) -> None:
    """RMSD считает сценарий 2; сценарий 3 только переносит их в факты.

    Нули позы 1 каждого кармана — точка отсчёта, а не измерение, и в факты
    не попадают, иначе разброс по карману окажется заниженным.
    """
    for item in analyses:
        context = {
            "ligand_id": item.key.ligand_id,
            "pocket_id": item.key.pocket_id,
            "pose_id": item.key.pose_id,
        }
        if item.key.pose_id != 1:
            for value, method in (
                (item.row.rmsd_lb, "pose_rmsd_lb_to_first"),
                (item.row.rmsd_ub, "pose_rmsd_ub_to_first"),
            ):
                if value is None:
                    continue
                state.add_evidence(
                    Evidence(
                        kind="calculation",
                        value=round(value, 3),
                        unit="Å",
                        method=method,
                        source="poses.csv",
                        context=context,
                    )
                )
        if item.row.rmsd_to_crystal is not None:
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=round(item.row.rmsd_to_crystal, 3),
                    unit="Å",
                    method="pose_rmsd_to_crystal",
                    source="poses.csv",
                    context={**context, "role": "control"},
                )
            )


def _add_contact_evidence(state: WorkflowState, analyses: list[PoseAnalysis], matrix) -> None:
    for item in analyses:
        context = {
            "ligand_id": item.key.ligand_id,
            "pocket_id": item.key.pocket_id,
            "pose_id": item.key.pose_id,
        }
        counts = {name: value for name, value in item.counts().items() if value}
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=len(item.contacts),
                method="pose_contact_count",
                source="pose_contacts.csv",
                context={**context, "by_type": counts},
            )
        )
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=len(item.residues()),
                method="pose_contact_residue_count",
                source="pose_contacts.csv",
                context=context,
            )
        )

    for change in matrix.changes:
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=change.residue_label,
                method=(
                    "contact_gained_in_pose" if change.change == "gained" else "contact_lost_in_pose"
                ),
                source="pose_contacts.csv",
                context={
                    "ligand_id": change.key.ligand_id,
                    "pocket_id": change.key.pocket_id,
                    "pose_id": change.key.pose_id,
                    "interaction_type": change.interaction_type,
                },
            )
        )


def _add_distance_evidence(state: WorkflowState, series, analyses: list[PoseAnalysis]) -> None:
    by_label = {item.key.label: item.key for item in analyses}
    for entry in series:
        for label, value in entry.values.items():
            key = by_label.get(label)
            if value is None or key is None:
                continue
            state.add_evidence(
                Evidence(
                    kind="calculation",
                    value=value,
                    unit="Å",
                    method="key_distance_by_pose",
                    source="protein_state.csv",
                    context={
                        "ligand_id": key.ligand_id,
                        "pocket_id": key.pocket_id,
                        "pose_id": key.pose_id,
                        "residue": entry.residue_label,
                        "ligand_atom": f"{entry.ligand_atom}#{entry.ligand_atom_serial}",
                        "protein_atom": f"{entry.protein_atom}#{entry.protein_atom_serial}",
                        "interaction_type": entry.interaction_type,
                    },
                )
            )


def _add_findings_evidence(state: WorkflowState, findings) -> None:
    for finding in findings:
        method = "pose_series_spike" if finding.finding in ("spike", "jump") else "pose_series_trend"
        value = finding.value if finding.finding in ("spike", "jump") else finding.finding
        state.add_evidence(
            Evidence(
                kind="calculation",
                value=value,
                method=method,
                source="protein_state.csv",
                context={
                    "series": finding.series,
                    "pose": finding.pose_label,
                    "baseline": finding.baseline,
                    "deviation": finding.deviation,
                    "detail": finding.detail,
                },
            )
        )


# --------------------------------------------------------------------------
# Вывод по сравнению
# --------------------------------------------------------------------------


def _candidates(
    analyses: list[PoseAnalysis], deltas, matrix, findings, reference: PoseKey
) -> list[dict]:
    """Позы, которые стоит проверить подробнее, с причинами. Считается кодом."""
    reasons: dict[str, list[str]] = {}
    reference_label = _reference_label(deltas)
    by_label = {item.key.label: item for item in analyses}

    for delta in deltas:
        label = delta.key.label
        if delta.is_reference:
            reasons.setdefault(label, []).append("лучшая оценка Vina в прогоне")
        elif abs(delta.delta_to_reference) <= 0.5:
            reasons.setdefault(label, []).append(
                "оценка Vina отличается от позы №1 всего на "
                f"{delta.delta_to_reference:+.2f} ккал/моль"
            )

    gained: dict[str, set[str]] = {}
    lost: dict[str, set[str]] = {}
    for change in matrix.changes:
        target = gained if change.change == "gained" else lost
        target.setdefault(change.key.label, set()).add(
            f"{change.residue_label} ({change.interaction_type})"
        )
    for label, items in gained.items():
        reasons.setdefault(label, []).append(
            "появляются контакты: " + ", ".join(sorted(items)[:3])
        )
    for label, items in lost.items():
        reasons.setdefault(label, []).append("исчезают контакты: " + ", ".join(sorted(items)[:3]))

    for finding in findings:
        if finding.pose_label:
            reasons.setdefault(finding.pose_label, []).append(finding.detail)

    options: list[dict] = []
    for label, items in reasons.items():
        pose = by_label.get(label)
        if pose is None:
            continue
        same_pocket = pose.key.pocket_id == reference.pocket_id
        options.append(
            {
                "pose": label,
                "pocket_id": pose.key.pocket_id,
                "pose_id": pose.key.pose_id,
                "is_reference": label == reference_label,
                "same_pocket": same_pocket,
                "contacts": len(pose.contacts),
                "reasons": "; ".join(items[:4]),
                "reason_codes": items[:4],
            }
        )
    # Порядок: сначала позы кармана-победителя и не сама поза №1. Поза из
    # другого кармана стоит в другом участке белка, и её «новые контакты» —
    # это смена участка, а не изменение взаимодействия в том же кармане.
    options.sort(
        key=lambda item: (
            not item["same_pocket"],
            item["is_reference"],
            -len(item["reason_codes"]),
            item["pose_id"],
        )
    )
    return options[:8]


def _conclude(state: WorkflowState, options: list[dict], deltas) -> None:
    if not options:
        return
    chosen = options[0]
    source = "rules"
    rationale = ""
    decision = decide(
        "Какую позу из сравнения стоит проверить подробнее",
        [
            {
                "pose": option["pose"],
                "contacts": option["contacts"],
                "reasons": option["reasons"],
            }
            for option in options
        ],
        criteria=(
            "Интереснее та поза, у которой меняется набор контактов с белком или "
            "обнаружен всплеск в рядах, при близкой оценке Vina. Числа не меняй и не добавляй."
        ),
        context=state.request.text,
    )
    if decision is not None:
        chosen = options[decision.index]
        source = decision.source
        rationale = decision.rationale

    state.add_evidence(
        Evidence(
            kind="hypothesis",
            value=chosen["pose"],
            method="pose_to_inspect",
            source=source,
            context={
                "pocket_id": chosen["pocket_id"],
                "pose_id": chosen["pose_id"],
                "reason_codes": chosen["reason_codes"],
                "rationale": rationale,
                "reference_pose": _reference_label(deltas),
            },
        )
    )


# --------------------------------------------------------------------------
# Точка входа
# --------------------------------------------------------------------------


def compare_poses(state: WorkflowState, csv_path: Path) -> WorkflowState:
    """Сравнение поз одного прогона: графики, таблицы свойств и вывод."""
    with state.agent_step("pose_comparison_agent", "compare_poses"):
        csv_path = Path(csv_path)
        run_dir = csv_path.parent
        rows = load_poses_csv(csv_path)
        if not rows:
            state.add_error(f"Таблица поз пуста: {csv_path}")
            return state

        meta = _load_run_meta(run_dir)
        protein_atoms, receptor_path = _load_receptor(run_dir, meta)
        if not protein_atoms:
            state.add_error(
                "Не найдены координаты белка для разбора контактов. "
                "Нужен receptor.pdb или receptor_pdb_path в run_meta.json."
            )
            return state
        state.observations["receptor_path"] = receptor_path

        analyses = _build_analyses(rows, protein_atoms, run_dir, state)
        if not analyses:
            state.add_error("Ни одного файла позы не удалось прочитать — сравнивать нечего.")
            return state
        if len(analyses) < 2:
            state.limitations.append(
                "В таблице одна поза: сравнение поз требует прогона, где Vina "
                "сохранила несколько режимов. Сравнение не выполнялось."
            )
            state.observations["pose_comparison"] = {"poses": 1, "compared": False}
            return state

        reference = _reference_key(meta, analyses)
        # На графики и в факты идёт читаемая выборка, в таблицы — все позы.
        selected = select_for_comparison(analyses, reference, max_poses=MAX_PLOT_POSES)
        deltas_all = vina_deltas(analyses, reference)
        deltas = vina_deltas(selected, reference)
        matrix = contact_matrix(selected, reference)
        reference_item = next(item for item in analyses if item.key == reference)
        pairs = key_atom_pairs(reference_item, limit=MAX_KEY_PAIRS)
        series = distance_series(analyses, protein_atoms, pairs)
        findings = collect_findings(deltas, series, selected)

        _add_score_evidence(state, deltas, selected)
        _add_rmsd_evidence(state, selected)
        _add_contact_evidence(state, selected, matrix)
        _add_distance_evidence(state, series, selected)
        _add_findings_evidence(state, findings)
        _conclude(state, _candidates(selected, deltas, matrix, findings, reference), deltas)

        if len(selected) < len(analyses):
            state.limitations.append(
                f"На графиках показаны {len(selected)} поз из {len(analyses)}: "
                f"карман позы №1 целиком и лучшая поза каждого из остальных карманов. "
                "Полные таблицы содержат все позы."
            )

        contacts_csv = write_csv(
            run_dir / "pose_contacts.csv",
            pose_metrics.CONTACT_CSV_FIELDS,
            contact_rows(analyses),
        )
        state_fields, state_rows = protein_state_rows(analyses, deltas_all, series)
        protein_csv = write_csv(run_dir / "protein_state.csv", state_fields, state_rows)
        state.artifacts["pose_contacts_csv"] = str(contacts_csv)
        state.artifacts["protein_state_csv"] = str(protein_csv)

        pose_labels = [delta.key.label for delta in deltas]
        state.artifacts["vina_by_pose_png"] = str(
            write_vina_by_pose(deltas, run_dir / "vina_by_pose.png")
        )
        state.artifacts["contact_map_png"] = str(
            write_contact_map(matrix, run_dir / "contact_map.png")
        )
        state.artifacts["key_distances_png"] = str(
            write_distance_series(series, pose_labels, run_dir / "key_distances.png")
        )

        state.observations["pose_comparison"] = {
            "poses": len(analyses),
            "pockets": len({item.key.pocket_id for item in analyses}),
            "compared": True,
            "reference_pose": reference.label,
            "contacts_total": sum(len(item.contacts) for item in analyses),
            "changed_contacts": len(matrix.changes),
            "findings": len(findings),
            "receptor": receptor_path,
        }

        state.limitations.append(
            "Контакты определяются по геометрии тяжёлых атомов и свойствам элементов "
            "без явного расчёта водородов: пороги 3,5 Å для водородной связи, 4,0 Å "
            "для солевого мостика и гидрофобного контакта, 2,8 Å для координации металла."
        )
        state.limitations.append(
            "Сравниваются позы одного прогона при одинаковых настройках скоринга "
            f"({meta.get('scoring', 'vina')}, exhaustiveness={meta.get('exhaustiveness', '?')}); "
            "разница оценок Vina не является разницей измеренной аффинности."
        )
        if any(item.row.rmsd_to_crystal is not None for item in analyses):
            state.limitations.append(
                "RMSD до кристаллической позы показан как контрольная метрика "
                "и в выбор позы не входит."
            )
    return state
