"""Сравнение поз между собой: дельты Vina, карта «поза × остаток», ряды расстояний.

Все числа здесь считаются по координатам и таблице докинга. Языковая модель
в этот модуль не входит: она может только выбрать из готового списка.
"""
from __future__ import annotations

import csv
import statistics
from dataclasses import dataclass, field
from pathlib import Path

from models.docking import PoseRow
from tools.contacts import (
    INTERACTION_TYPES,
    Atom,
    Contact,
    contacts_by_residue,
    count_by_type,
)

# Взаимодействия, по которым выбираются ключевые остатки: от специфичных к общим.
_KEY_INTERACTION_ORDER = (
    "metal",
    "salt_bridge",
    "hydrogen_bond",
    "halogen_bond",
    "pi_stacking",
    "pi_cation",
)

SPIKE_SIGMA = 3.0
MIN_SERIES_LEN = 3


@dataclass(frozen=True)
class PoseKey:
    """Тройка, которой сценарий 2 и сценарий 3 связывают свои таблицы."""

    ligand_id: str
    pocket_id: str
    pose_id: int

    @property
    def label(self) -> str:
        return f"{self.pocket_id}/поза {self.pose_id}"

    @classmethod
    def from_row(cls, row: PoseRow) -> "PoseKey":
        return cls(ligand_id=row.ligand_id, pocket_id=row.pocket_id, pose_id=row.pose_id)


@dataclass
class PoseAnalysis:
    """Одна поза: строка таблицы докинга, её атомы и найденные контакты."""

    key: PoseKey
    row: PoseRow
    ligand_atoms: list[Atom]
    contacts: list[Contact]

    @property
    def vina_score(self) -> float:
        return self.row.vina_score

    def residues(self) -> dict[str, set[str]]:
        return contacts_by_residue(self.contacts)

    def counts(self) -> dict[str, int]:
        return count_by_type(self.contacts)


@dataclass
class VinaDelta:
    key: PoseKey
    vina_score: float
    delta_to_reference: float
    is_reference: bool


@dataclass
class ContactChange:
    """Контакт, появившийся или исчезнувший относительно базовой позы кармана."""

    key: PoseKey
    residue_label: str
    interaction_type: str
    change: str  # "gained" | "lost"
    baseline: str = ""


@dataclass
class ContactMatrix:
    """Карта «поза × остаток»: в какой позе какие взаимодействия есть."""

    poses: list[PoseKey]
    residues: list[str]
    cells: dict[tuple[str, str], set[str]] = field(default_factory=dict)
    reference: PoseKey | None = None
    changes: list[ContactChange] = field(default_factory=list)
    baselines: dict[str, PoseKey] = field(default_factory=dict)

    def types_at(self, pose: PoseKey, residue: str) -> set[str]:
        return self.cells.get((pose.label, residue), set())

    def baseline_for(self, pose: PoseKey) -> PoseKey | None:
        return self.baselines.get(pose.pocket_id)


@dataclass
class DistanceSeries:
    """Одна и та же пара атомов, измеренная во всех позах."""

    ligand_atom: str
    ligand_atom_serial: int
    protein_atom: str
    protein_atom_serial: int
    residue_label: str
    interaction_type: str
    values: dict[str, float | None] = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.residue_label} {self.protein_atom} – лиганд {self.ligand_atom}#{self.ligand_atom_serial}"


@dataclass
class SeriesFinding:
    """Закономерность в ряду: всплеск, тренд или скачок."""

    series: str
    finding: str  # "spike" | "trend_up" | "trend_down" | "jump"
    pose_label: str | None
    value: float
    baseline: float
    deviation: float
    detail: str


# --------------------------------------------------------------------------
# Дельты Vina
# --------------------------------------------------------------------------


def best_pose_key(analyses: list[PoseAnalysis]) -> PoseKey | None:
    """Поза №1 по Vina — минимальная оценка. Используется, если run_meta её не задал."""
    if not analyses:
        return None
    best = min(analyses, key=lambda item: item.vina_score)
    return best.key


def vina_deltas(analyses: list[PoseAnalysis], reference: PoseKey) -> list[VinaDelta]:
    """Оценки всех поз и разница относительно позы отсчёта, ккал/моль."""
    base = next((item.vina_score for item in analyses if item.key == reference), None)
    if base is None:
        base = min(item.vina_score for item in analyses) if analyses else 0.0
    deltas = []
    for item in analyses:
        deltas.append(
            VinaDelta(
                key=item.key,
                vina_score=item.vina_score,
                delta_to_reference=round(item.vina_score - base, 3),
                is_reference=item.key == reference,
            )
        )
    deltas.sort(key=lambda d: (d.key.pocket_id, d.key.pose_id))
    return deltas


# --------------------------------------------------------------------------
# Карта «поза × остаток»
# --------------------------------------------------------------------------


def _residue_sort_key(label: str) -> tuple[str, int]:
    digits = "".join(ch for ch in label if ch.isdigit())
    letters = "".join(ch for ch in label if not ch.isdigit())
    return (letters, int(digits) if digits else 0)


def pocket_baselines(analyses: list[PoseAnalysis]) -> dict[str, PoseKey]:
    """Базовая поза каждого кармана — его первая поза по Vina.

    Контакты сравниваются внутри кармана, а не через весь прогон. Поза из
    другого кармана стоит в другом месте белка, поэтому «исчезли все контакты»
    там означало бы лишь смену участка, а не изменение взаимодействий.
    """
    baselines: dict[str, PoseKey] = {}
    for item in sorted(analyses, key=lambda entry: (entry.key.pocket_id, entry.key.pose_id)):
        baselines.setdefault(item.key.pocket_id, item.key)
    return baselines


def contact_matrix(analyses: list[PoseAnalysis], reference: PoseKey) -> ContactMatrix:
    """Какие взаимодействия появляются и исчезают при переходе между позами.

    Отсчёт для каждой позы — первая поза её собственного кармана; `reference`
    остаётся позой №1 всего прогона и отмечается на карте отдельно.
    """
    cells: dict[tuple[str, str], set[str]] = {}
    residues: set[str] = set()
    for item in analyses:
        for residue, types in item.residues().items():
            residues.add(residue)
            cells[(item.key.label, residue)] = set(types)

    ordered_residues = sorted(residues, key=_residue_sort_key)
    poses = [item.key for item in analyses]

    baselines = pocket_baselines(analyses)
    by_key = {item.key: item for item in analyses}
    changes: list[ContactChange] = []
    for item in analyses:
        baseline_key = baselines.get(item.key.pocket_id)
        if baseline_key is None or baseline_key == item.key:
            continue
        baseline_item = by_key.get(baseline_key)
        if baseline_item is None:
            continue
        before_map = baseline_item.residues()
        current = item.residues()
        for residue in ordered_residues:
            before = before_map.get(residue, set())
            after = current.get(residue, set())
            for interaction in sorted(after - before):
                changes.append(
                    ContactChange(item.key, residue, interaction, "gained", baseline_key.label)
                )
            for interaction in sorted(before - after):
                changes.append(
                    ContactChange(item.key, residue, interaction, "lost", baseline_key.label)
                )

    return ContactMatrix(
        poses=poses,
        residues=ordered_residues,
        cells=cells,
        reference=reference,
        changes=changes,
        baselines=baselines,
    )


def select_for_comparison(
    analyses: list[PoseAnalysis], reference: PoseKey, max_poses: int = 12
) -> list[PoseAnalysis]:
    """Позы, которые идут на графики: карман позы №1 целиком плюс лучшие из прочих.

    Все 27 поз на одной карте нечитаемы, а сравнивать интересно прежде всего
    позы внутри кармана-победителя. Карманы остальных остаются в виде своей
    лучшей позы, чтобы разница между участками была видна.
    """
    if len(analyses) <= max_poses:
        return list(analyses)
    same_pocket = [item for item in analyses if item.key.pocket_id == reference.pocket_id]
    same_pocket.sort(key=lambda item: item.key.pose_id)
    others: dict[str, PoseAnalysis] = {}
    for item in analyses:
        if item.key.pocket_id == reference.pocket_id:
            continue
        current = others.get(item.key.pocket_id)
        if current is None or item.vina_score < current.vina_score:
            others[item.key.pocket_id] = item
    selected = same_pocket[: max_poses - len(others)] + [
        others[pocket] for pocket in sorted(others)
    ]
    selected.sort(key=lambda item: (item.key.pocket_id, item.key.pose_id))
    return selected


# --------------------------------------------------------------------------
# Расстояния до ключевых остатков
# --------------------------------------------------------------------------


def key_atom_pairs(reference: PoseAnalysis, limit: int = 6) -> list[Contact]:
    """Пары атомов из позы отсчёта, которые дальше измеряются во всех позах.

    Берутся самые специфичные взаимодействия и по одному представителю на
    остаток: график расстояний должен быть читаемым, а не полным перечнем.
    """
    order = {name: index for index, name in enumerate(_KEY_INTERACTION_ORDER)}
    candidates = [c for c in reference.contacts if c.interaction_type in order]
    if not candidates:
        candidates = list(reference.contacts)
    candidates.sort(
        key=lambda c: (order.get(c.interaction_type, 99), c.distance_a)
    )
    picked: list[Contact] = []
    seen_residues: set[str] = set()
    for contact in candidates:
        if contact.residue_label in seen_residues:
            continue
        seen_residues.add(contact.residue_label)
        picked.append(contact)
        if len(picked) >= limit:
            break
    return picked


def distance_series(
    analyses: list[PoseAnalysis],
    protein_atoms: list[Atom],
    pairs: list[Contact],
) -> list[DistanceSeries]:
    """Те же пары атомов во всех позах. Отсутствие атома даёт None, а не ноль."""
    protein_by_serial = {atom.serial: atom for atom in protein_atoms}
    series: list[DistanceSeries] = []
    for pair in pairs:
        protein_atom = protein_by_serial.get(pair.protein_atom_serial)
        entry = DistanceSeries(
            ligand_atom=pair.ligand_atom,
            ligand_atom_serial=pair.ligand_atom_serial,
            protein_atom=pair.protein_atom,
            protein_atom_serial=pair.protein_atom_serial,
            residue_label=pair.residue_label,
            interaction_type=pair.interaction_type,
        )
        for item in analyses:
            ligand_atom = next(
                (a for a in item.ligand_atoms if a.serial == pair.ligand_atom_serial),
                None,
            )
            if ligand_atom is None or protein_atom is None:
                entry.values[item.key.label] = None
            else:
                entry.values[item.key.label] = round(ligand_atom.distance_to(protein_atom), 3)
        series.append(entry)
    return series


# --------------------------------------------------------------------------
# Закономерности в рядах
# --------------------------------------------------------------------------


def series_findings(name: str, labels: list[str], values: list[float | None]) -> list[SeriesFinding]:
    """Всплески, монотонные тренды и скачки в одном ряду.

    Всплеск определяется по медиане и медианному абсолютному отклонению:
    среднее и обычная сигма на трёх-девяти точках неустойчивы.
    """
    pairs = [(label, value) for label, value in zip(labels, values) if value is not None]
    if len(pairs) < MIN_SERIES_LEN:
        return []
    numbers = [value for _, value in pairs]
    findings: list[SeriesFinding] = []

    median = statistics.median(numbers)
    deviations = [abs(value - median) for value in numbers]
    mad = statistics.median(deviations)
    scale = mad * 1.4826 if mad > 0 else 0.0
    if scale > 0:
        for label, value in pairs:
            score = abs(value - median) / scale
            if score >= SPIKE_SIGMA:
                findings.append(
                    SeriesFinding(
                        series=name,
                        finding="spike",
                        pose_label=label,
                        value=round(value, 3),
                        baseline=round(median, 3),
                        deviation=round(value - median, 3),
                        detail=(
                            f"{label}: {value:.2f} против медианы {median:.2f} "
                            f"(отклонение {score:.1f} робастных сигм)"
                        ),
                    )
                )

    increasing = all(b > a for a, b in zip(numbers, numbers[1:]))
    decreasing = all(b < a for a, b in zip(numbers, numbers[1:]))
    if increasing or decreasing:
        findings.append(
            SeriesFinding(
                series=name,
                finding="trend_up" if increasing else "trend_down",
                pose_label=None,
                value=round(numbers[-1], 3),
                baseline=round(numbers[0], 3),
                deviation=round(numbers[-1] - numbers[0], 3),
                detail=(
                    f"ряд монотонно {'растёт' if increasing else 'падает'} "
                    f"от {numbers[0]:.2f} до {numbers[-1]:.2f} по {len(numbers)} позам"
                ),
            )
        )

    steps = [abs(b - a) for a, b in zip(numbers, numbers[1:])]
    if len(steps) >= 2:
        median_step = statistics.median(steps)
        if median_step > 0:
            largest = max(range(len(steps)), key=lambda index: steps[index])
            if steps[largest] >= 3.0 * median_step:
                findings.append(
                    SeriesFinding(
                        series=name,
                        finding="jump",
                        pose_label=pairs[largest + 1][0],
                        value=round(numbers[largest + 1], 3),
                        baseline=round(numbers[largest], 3),
                        deviation=round(numbers[largest + 1] - numbers[largest], 3),
                        detail=(
                            f"скачок между {pairs[largest][0]} и {pairs[largest + 1][0]}: "
                            f"{numbers[largest]:.2f} → {numbers[largest + 1]:.2f}"
                        ),
                    )
                )
    return findings


def collect_findings(
    deltas: list[VinaDelta],
    series: list[DistanceSeries],
    analyses: list[PoseAnalysis],
) -> list[SeriesFinding]:
    """Все закономерности по оценкам Vina, расстояниям и числу контактов."""
    labels = [delta.key.label for delta in deltas]
    findings = series_findings(
        "Оценка Vina, ккал/моль", labels, [delta.vina_score for delta in deltas]
    )
    counts_by_label = {item.key.label: len(item.contacts) for item in analyses}
    findings.extend(
        series_findings(
            "Число контактов", labels, [counts_by_label.get(label) for label in labels]
        )
    )
    for entry in series:
        findings.extend(
            series_findings(
                f"Расстояние {entry.label}, Å",
                labels,
                [entry.values.get(label) for label in labels],
            )
        )
    return findings


# --------------------------------------------------------------------------
# Таблицы для записи на диск
# --------------------------------------------------------------------------

CONTACT_CSV_FIELDS = (
    "ligand_id",
    "pocket_id",
    "pose_id",
    "residue_name",
    "residue_seq",
    "chain",
    "interaction_type",
    "ligand_atom",
    "ligand_atom_serial",
    "protein_atom",
    "protein_atom_serial",
    "distance_a",
    "pose_file",
)


def contact_rows(analyses: list[PoseAnalysis]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in analyses:
        for contact in item.contacts:
            rows.append(
                {
                    "ligand_id": item.key.ligand_id,
                    "pocket_id": item.key.pocket_id,
                    "pose_id": item.key.pose_id,
                    "residue_name": contact.residue_name,
                    "residue_seq": contact.residue_seq,
                    "chain": contact.chain,
                    "interaction_type": contact.interaction_type,
                    "ligand_atom": contact.ligand_atom,
                    "ligand_atom_serial": contact.ligand_atom_serial,
                    "protein_atom": contact.protein_atom,
                    "protein_atom_serial": contact.protein_atom_serial,
                    "distance_a": contact.distance_a,
                    "pose_file": item.row.pose_file,
                }
            )
    return rows


def protein_state_rows(
    analyses: list[PoseAnalysis],
    deltas: list[VinaDelta],
    series: list[DistanceSeries],
) -> tuple[list[str], list[dict[str, object]]]:
    """Широкая таблица «поза × свойство» — тот самый набор свойств белка по позе."""
    delta_by_label = {delta.key.label: delta for delta in deltas}
    fields = [
        "ligand_id",
        "pocket_id",
        "pose_id",
        "vina_score",
        "vina_delta_to_reference",
        "rmsd_lb",
        "rmsd_ub",
        "rmsd_to_crystal",
        "contact_count",
        "residue_count",
    ]
    fields.extend(f"n_{name}" for name in INTERACTION_TYPES)
    distance_fields = [f"d_{entry.residue_label}_{entry.protein_atom}" for entry in series]
    fields.extend(distance_fields)

    rows: list[dict[str, object]] = []
    for item in analyses:
        label = item.key.label
        delta = delta_by_label.get(label)
        counts = item.counts()
        row: dict[str, object] = {
            "ligand_id": item.key.ligand_id,
            "pocket_id": item.key.pocket_id,
            "pose_id": item.key.pose_id,
            "vina_score": item.vina_score,
            "vina_delta_to_reference": delta.delta_to_reference if delta else "",
            "rmsd_lb": _blank(getattr(item.row, "rmsd_lb", None)),
            "rmsd_ub": _blank(getattr(item.row, "rmsd_ub", None)),
            "rmsd_to_crystal": _blank(item.row.rmsd_to_crystal),
            "contact_count": len(item.contacts),
            "residue_count": len(item.residues()),
        }
        for name in INTERACTION_TYPES:
            row[f"n_{name}"] = counts.get(name, 0)
        for entry, column in zip(series, distance_fields):
            row[column] = _blank(entry.values.get(label))
        rows.append(row)
    return fields, rows


def _blank(value: float | None) -> object:
    """Нет данных — пустая ячейка. Ноль означал бы измеренный ноль."""
    return "" if value is None else value


def write_csv(path: Path, fields: list[str] | tuple[str, ...], rows: list[dict[str, object]]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
