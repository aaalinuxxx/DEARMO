"""Геометрия позы: проверки в духе PoseBusters, посчитанные здесь.

Пакет posebusters в окружении не установлен, а сценарию 4 его проверки нужны
целиком. Главное требование кейса — вердикт обязан называть пару атомов,
измеренное расстояние и порог, а не только статус «прошла / не прошла».
Инструмент, возвращающий один флаг, этого не даёт, поэтому каждая проверка
считается здесь и складывает то, что реально померила, в `Measurement`.

Откуда берутся эталонные числа:

* радиусы ван-дер-Ваальса — из таблицы элементов RDKit;
* допустимые длины связей и валентные углы — из матрицы границ дистанционной
  геометрии, построенной RDKit для той же самой молекулы;
* эталон координации металла — измеряется на кристаллической структуре, если
  она передана.

Ни одно из этих значений не зашито в код списком «правильных» ответов: иначе
проверка подтверждала бы сама себя. Единственные константы здесь — пороги,
то есть правила приёмки, и они собраны вверху файла под своими именами.

Парсер координат не дублируется: `Atom`, `read_structure` и `heavy_atoms`
берутся из `tools.contacts`. В PDBQT колонки 77-78 заняты типом AutoDock,
а не элементом, и файл позы от Vina бывает многомодельным — и то и другое
там уже разобрано.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from tools.contacts import Atom, heavy_atoms, read_structure

try:  # pragma: no cover - окружение без C++ RDKit
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom, rdMolTransforms

    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    rdDistGeom = None  # type: ignore[assignment]
    rdMolTransforms = None  # type: ignore[assignment]
    _RDKIT_AVAILABLE = False


# --------------------------------------------------------------------------
# Пороги приёмки. Это правила, а не измерения: измерения считаются по данным.
# --------------------------------------------------------------------------

#: Порог PoseBusters для пары белок-лиганд: доля от суммы радиусов
#: ван-дер-Ваальса двух тяжёлых атомов. Ближе — столкновение.
PROTEIN_LIGAND_CLASH_FRACTION = 0.75

#: Тот же принцип для несвязанных атомов внутри лиганда. Порог мягче:
#: внутри молекулы сближение ограничено ещё и её собственной связностью.
INTERNAL_CLASH_FRACTION = 0.8

#: Сколько связей должно разделять атомы лиганда, чтобы считать пару
#: несвязанной. Соседи через одну и две связи сближены самой геометрией
#: валентных углов, и сравнивать их с суммой радиусов бессмысленно.
INTERNAL_CLASH_MIN_BOND_SEPARATION = 4

#: Допуск к границам дистанционной геометрии. Границы задают диапазон,
#: достижимый для этой молекулы; докинг слегка выходит за него без потери
#: химического смысла, поэтому диапазон расширяется на эту долю.
DISTANCE_GEOMETRY_TOLERANCE = 0.25

#: Максимальное отклонение атома ароматического кольца от плоскости кольца.
RING_PLANARITY_TOLERANCE_A = 0.25

#: Расстояние, внутри которого атом лиганда считается координирующим металл.
METAL_COORDINATION_MAX_A = 2.8

#: Допустимое отклонение расстояния до металла от кристаллического значения.
METAL_REFERENCE_TOLERANCE_A = 0.5

#: Запас к габаритному ящику лиганда при отборе атомов белка. Перебор пар
#: по всему белку (тысячи атомов) не нужен: столкновение всегда рядом.
_NEIGHBOR_BOX_MARGIN_A = 6.0

#: Элементы лиганда, способные координировать металл.
_METAL_DONOR_ELEMENTS = frozenset({"N", "O", "S", "SE", "P"})

#: Запасные радиусы ван-дер-Ваальса (Å) на случай, если RDKit недоступен.
#: Значения Бонди; используются только чтобы проверка столкновений не
#: отключалась вместе с RDKit, и факт замены попадает в примечания отчёта.
_FALLBACK_VDW_A = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90, "B": 1.92,
    "SI": 2.10, "ZN": 1.39, "MG": 1.73, "MN": 2.05, "FE": 2.05, "CU": 1.40,
    "CO": 2.00, "NI": 1.63, "CA": 2.31, "NA": 2.27, "K": 2.75,
}

#: Статусы проверки. «Не применимо» — полноценный исход, а не неудача:
#: у молекулы без заданной стереохимии нечего проверять, и считать это
#: провалом значило бы занижать процент пройденности.
PASSED = "passed"
FAILED = "failed"
NOT_APPLICABLE = "not_applicable"

#: Вердикт по позе целиком.
VERDICT_PASSED = "poza_prinyata"
VERDICT_WARNINGS = "poza_prinyata_s_zamechaniyami"
VERDICT_REJECTED = "poza_otklonena"

_VERDICT_TEXT = {
    VERDICT_PASSED: "поза принята",
    VERDICT_WARNINGS: "поза принята с замечаниями",
    VERDICT_REJECTED: "поза отклонена",
}


# --------------------------------------------------------------------------
# Результаты
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Measurement:
    """Одно измерение, сохранённое на этапе расчёта.

    Кейс требует объяснимого вердикта: мало сказать «не прошла», нужно
    показать, что именно померили, какое значение получили и какой порог
    нарушен. Поэтому измерение хранится даже тогда, когда проверка пройдена.
    """

    label: str
    value: float | str | None
    unit: str = "Å"
    low: float | None = None
    high: float | None = None
    passed: bool = True
    detail: str = ""

    @property
    def range_text(self) -> str:
        """Допустимый диапазон словами. Односторонний порог тоже читается."""
        if self.low is not None and self.high is not None:
            return f"{self.low:.3g}–{self.high:.3g} {self.unit}"
        if self.low is not None:
            return f"не менее {self.low:.3g} {self.unit}"
        if self.high is not None:
            return f"не более {self.high:.3g} {self.unit}"
        return "диапазон не задан"

    @property
    def value_text(self) -> str:
        if self.value is None:
            return "не измерено"
        if isinstance(self.value, str):
            return self.value
        return f"{self.value:.3g} {self.unit}".strip()

    def as_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "value": self.value,
            "unit": self.unit,
            "low": self.low,
            "high": self.high,
            "passed": self.passed,
            "detail": self.detail,
            "range": self.range_text,
        }

    def __str__(self) -> str:
        head = f"{self.label}: {self.value_text}"
        if self.low is not None or self.high is not None:
            head += f" (допустимо {self.range_text})"
        return f"{head}. {self.detail}".strip().rstrip(".") + "."


@dataclass
class GeometryCheck:
    """Одна проверка из таблицы сценария 4."""

    name: str
    title: str
    status: str
    summary: str
    measurements: list[Measurement] = field(default_factory=list)
    #: Критическая проверка отменяет общий вердикт независимо от того,
    #: сколько других проверок пройдено. Так столкновение нельзя «перевесить»
    #: четырьмя удачными проверками из пяти.
    critical: bool = False

    @property
    def applicable(self) -> bool:
        return self.status != NOT_APPLICABLE

    @property
    def passed(self) -> bool:
        return self.status == PASSED

    @property
    def failures(self) -> list[Measurement]:
        return [item for item in self.measurements if not item.passed]

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "title": self.title,
            "status": self.status,
            "summary": self.summary,
            "critical": self.critical,
            "measurements": [item.as_dict() for item in self.measurements],
        }


@dataclass
class GeometryReport:
    """Итог сценария 4: вердикт, счёт и список критических нарушений."""

    checks: list[GeometryCheck] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def applicable(self) -> list[GeometryCheck]:
        return [check for check in self.checks if check.applicable]

    @property
    def passed_checks(self) -> list[GeometryCheck]:
        return [check for check in self.checks if check.passed]

    @property
    def score_text(self) -> str:
        """«4 из 5» — пройдено из применимых, а не из всех написанных."""
        return f"{len(self.passed_checks)} из {len(self.applicable)}"

    @property
    def percent_passed(self) -> float:
        """Процент пройденности считается от применимых проверок.

        Неприменимые в знаменатель не входят: молекула без хиральных центров
        не должна получать штраф за то, что их не у чего проверять.
        """
        applicable = self.applicable
        if not applicable:
            return 0.0
        return round(100.0 * len(self.passed_checks) / len(applicable), 1)

    @property
    def critical_violations(self) -> list[str]:
        violations: list[str] = []
        for check in self.checks:
            if not (check.critical and check.status == FAILED):
                continue
            if check.failures:
                violations.extend(f"{check.title}: {item}" for item in check.failures)
            else:
                # Критическая проверка может провалиться, не дойдя до измерений:
                # например, файла позы нет. Без этой ветки такой провал не попадал
                # бы в список нарушений и вердикт смягчался бы до «с замечаниями».
                violations.append(f"{check.title}: {check.summary}")
        return violations

    @property
    def verdict(self) -> str:
        if self.critical_violations:
            return VERDICT_REJECTED
        if any(check.status == FAILED for check in self.checks):
            return VERDICT_WARNINGS
        return VERDICT_PASSED

    @property
    def verdict_text(self) -> str:
        return _VERDICT_TEXT[self.verdict]

    def check(self, name: str) -> GeometryCheck | None:
        for item in self.checks:
            if item.name == name:
                return item
        return None

    def as_dict(self) -> dict[str, Any]:
        return {
            "verdict": self.verdict,
            "verdict_text": self.verdict_text,
            "score": self.score_text,
            "percent_passed": self.percent_passed,
            "critical_violations": self.critical_violations,
            "checks": [check.as_dict() for check in self.checks],
            "notes": list(self.notes),
        }


# --------------------------------------------------------------------------
# Радиусы и вспомогательная геометрия
# --------------------------------------------------------------------------


def vdw_radius(element: str) -> float:
    """Радиус ван-дер-Ваальса элемента из таблицы RDKit, Å."""
    symbol = (element or "C").strip()
    if _RDKIT_AVAILABLE:
        try:
            table = Chem.GetPeriodicTable()
            radius = table.GetRvdw(symbol.capitalize())
            if radius and radius > 0:
                return float(radius)
        except Exception:  # noqa: BLE001 — неизвестный символ, идём в запасную таблицу
            pass
    return _FALLBACK_VDW_A.get(symbol.upper(), 1.70)


def clash_threshold(first: str, second: str, fraction: float) -> float:
    """Порог сближения пары элементов: доля от суммы их радиусов, Å."""
    return fraction * (vdw_radius(first) + vdw_radius(second))


def _bounding_box(atoms: Iterable[Atom], margin: float) -> tuple[float, float, float, float, float, float]:
    xs = [atom.x for atom in atoms]
    ys = [atom.y for atom in atoms]
    zs = [atom.z for atom in atoms]
    return (
        min(xs) - margin, max(xs) + margin,
        min(ys) - margin, max(ys) + margin,
        min(zs) - margin, max(zs) + margin,
    )


def nearby_protein_atoms(ligand: list[Atom], protein: list[Atom], margin: float = _NEIGHBOR_BOX_MARGIN_A) -> list[Atom]:
    """Атомы белка в габаритном ящике лиганда.

    Перебор пар по всему белку (у карбоангидразы это тысячи атомов) считает
    то же самое, но заметно дольше: столкновение по определению находится
    рядом с лигандом.
    """
    if not ligand or not protein:
        return []
    x0, x1, y0, y1, z0, z1 = _bounding_box(ligand, margin)
    return [
        atom
        for atom in protein
        if x0 <= atom.x <= x1 and y0 <= atom.y <= y1 and z0 <= atom.z <= z1
    ]


def _angle_deg(a: float, b: float, c: float) -> float | None:
    """Угол между сторонами a и b треугольника с третьей стороной c."""
    if a <= 0 or b <= 0:
        return None
    cosine = (a * a + b * b - c * c) / (2 * a * b)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _plane_deviation(points: list[tuple[float, float, float]]) -> float:
    """Максимальное отклонение точек от их плоскости наименьших квадратов, Å."""
    count = len(points)
    cx = sum(p[0] for p in points) / count
    cy = sum(p[1] for p in points) / count
    cz = sum(p[2] for p in points) / count
    # Нормаль методом Ньюэлла: устойчива и не требует numpy.
    nx = ny = nz = 0.0
    for index in range(count):
        ax, ay, az = points[index]
        bx, by, bz = points[(index + 1) % count]
        nx += (ay - by) * (az + bz)
        ny += (az - bz) * (ax + bx)
        nz += (ax - bx) * (ay + by)
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-9:
        return float("inf")
    nx, ny, nz = nx / norm, ny / norm, nz / norm
    return max(
        abs((p[0] - cx) * nx + (p[1] - cy) * ny + (p[2] - cz) * nz)
        for p in points
    )


# --------------------------------------------------------------------------
# Молекула позы
# --------------------------------------------------------------------------


def _pdb_block(atoms: list[Atom], resname: str = "LIG") -> str:
    """PDB-блок с координатами позы.

    Имена атомов перенумеровываются: в PDBQT они повторяются, а RDKit по
    повторяющимся именам в одном остатке строит неверную связность.
    """
    lines = []
    for index, atom in enumerate(atoms, start=1):
        element = atom.element.upper()[:2]
        name = f"{element}{index}"[:4]
        # Колонки PDB: имя атома 13-16, altLoc 17, остаток 18-20, цепь 22,
        # номер остатка 23-26, координаты с 31, элемент 77-78. Для элемента
        # из одной буквы имя начинается с 14-й колонки.
        name_field = f" {name:<3}" if len(element) == 1 else f"{name:<4}"
        lines.append(
            f"HETATM{index:>5} {name_field} {resname:>3} A{1:>4}    "
            f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {element:>2}"
        )
    lines.append("END")
    return "\n".join(lines)


def build_pose_molecule(atoms: list[Atom], smiles: str | None) -> tuple[Any | None, str | None]:
    """Молекула RDKit с координатами позы и порядками связей из SMILES прогона.

    Из PDBQT порядки связей не восстанавливаются: там нет ни кратностей, ни
    ароматичности. Поэтому связность берётся по расстояниям, а химия —
    из SMILES того же прогона (run_meta.json). Без этого шаблона проверки
    стереохимии и планарности колец проверяли бы не ту молекулу.

    Возвращает (молекула, примечание). Примечание не пустое, если молекулу
    удалось собрать лишь частично — оно идёт в отчёт, а не теряется.
    """
    if not _RDKIT_AVAILABLE:
        return None, "RDKit недоступен: проверки по молекуле не выполнялись."
    heavy = heavy_atoms(atoms)
    if len(heavy) < 2:
        return None, f"В позе {len(heavy)} тяжёлых атомов: молекулу собрать нельзя."

    raw = Chem.MolFromPDBBlock(_pdb_block(heavy), sanitize=False, removeHs=False, proximityBonding=True)
    if raw is None:
        return None, "RDKit не прочитал координаты позы."
    try:
        Chem.SanitizeMol(raw)
    except Exception as exc:  # noqa: BLE001 — связность по расстояниям бывает невалентной
        return None, f"Связность позы не прошла проверку валентностей: {exc}"

    if not smiles:
        return raw, "SMILES прогона не передан: порядки связей и ароматичность не восстановлены."

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return raw, f"SMILES прогона не разобран ({smiles[:40]}): порядки связей не восстановлены."
    template = Chem.RemoveHs(template)
    if template.GetNumAtoms() != raw.GetNumAtoms():
        return raw, (
            f"В позе {raw.GetNumAtoms()} тяжёлых атомов, в SMILES прогона "
            f"{template.GetNumAtoms()}: порядки связей не восстановлены."
        )
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
    except Exception as exc:  # noqa: BLE001 — поза не совпала с шаблоном
        return raw, f"Поза не легла на SMILES прогона ({exc}): порядки связей не восстановлены."
    return mol, None


# --------------------------------------------------------------------------
# Проверки
# --------------------------------------------------------------------------


def check_protein_ligand_clashes(
    ligand: list[Atom],
    protein: list[Atom],
    *,
    fraction: float = PROTEIN_LIGAND_CLASH_FRACTION,
) -> GeometryCheck:
    """Слишком близкие пары тяжёлых атомов белка и лиганда.

    Атомы металлов из перебора исключены: координационная связь короче любого
    порога по ван-дер-Ваальсу, и иначе правильно скоординированный цинк
    всегда объявлялся бы столкновением. Геометрию металла проверяет
    `check_metal_coordination`.
    """
    title = "Столкновения лиганда с белком"
    ligand_heavy = heavy_atoms(ligand)
    protein_heavy = [atom for atom in heavy_atoms(protein) if not atom.is_metal]
    if not ligand_heavy or not protein_heavy:
        return GeometryCheck(
            name="protein_ligand_clash",
            title=title,
            status=NOT_APPLICABLE,
            summary="Нет координат белка или лиганда — сравнивать нечего.",
            critical=True,
        )

    candidates = nearby_protein_atoms(ligand_heavy, protein_heavy)
    violations: list[Measurement] = []
    closest: tuple[float, Atom, Atom, float] | None = None
    for lig in ligand_heavy:
        for prot in candidates:
            distance = lig.distance_to(prot)
            threshold = clash_threshold(lig.element, prot.element, fraction)
            margin = distance - threshold
            if closest is None or margin < closest[0]:
                closest = (margin, lig, prot, distance)
            if distance < threshold:
                violations.append(
                    Measurement(
                        label=f"{lig.element}{lig.serial} — {prot.residue_label} {prot.name}",
                        value=round(distance, 3),
                        low=round(threshold, 3),
                        passed=False,
                        detail=(
                            f"Расстояние меньше {fraction:g} суммы радиусов "
                            f"ван-дер-Ваальса ({vdw_radius(lig.element):.2f} + "
                            f"{vdw_radius(prot.element):.2f} Å)"
                        ),
                    )
                )

    measurements = sorted(violations, key=lambda item: item.value or 0.0)
    if not violations and closest is not None:
        # Пройденная проверка без единого числа — это голое «ок», которое не с
        # чем сопоставить, поэтому самый тесный контакт сохраняется и тогда.
        # Если нарушения есть, самое тесное из них и есть этот контакт, и
        # повторять его отдельной строкой значило бы задвоить вердикт.
        margin, lig, prot, distance = closest
        measurements.append(
            Measurement(
                label=f"Самый тесный контакт: {lig.element}{lig.serial} — {prot.residue_label} {prot.name}",
                value=round(distance, 3),
                low=round(clash_threshold(lig.element, prot.element, fraction), 3),
                passed=True,
                detail=f"Запас до порога {margin:+.3f} Å",
            )
        )

    if violations:
        worst = measurements[0]
        return GeometryCheck(
            name="protein_ligand_clash",
            title=title,
            status=FAILED,
            summary=f"Столкновений: {len(violations)}. Самое тесное — {worst}",
            measurements=measurements,
            critical=True,
        )
    return GeometryCheck(
        name="protein_ligand_clash",
        title=title,
        status=PASSED,
        summary=f"Столкновений нет, проверено пар атомов: {len(ligand_heavy) * len(candidates)}.",
        measurements=measurements,
        critical=True,
    )


def check_internal_clashes(mol: Any, *, fraction: float = INTERNAL_CLASH_FRACTION) -> GeometryCheck:
    """Чрезмерное сближение несвязанных атомов внутри лиганда."""
    title = "Столкновения внутри лиганда"
    if mol is None or not _RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
        return GeometryCheck(
            name="internal_clash",
            title=title,
            status=NOT_APPLICABLE,
            summary="Молекула позы не собрана — внутренние расстояния не считались.",
            critical=True,
        )

    conf = mol.GetConformer()
    topological = Chem.GetDistanceMatrix(mol)
    violations: list[Measurement] = []
    closest: tuple[float, str, float, float] | None = None
    count = mol.GetNumAtoms()
    for i in range(count):
        for j in range(i + 1, count):
            if topological[i][j] < INTERNAL_CLASH_MIN_BOND_SEPARATION:
                continue
            first = mol.GetAtomWithIdx(i)
            second = mol.GetAtomWithIdx(j)
            distance = rdMolTransforms.GetBondLength(conf, i, j)
            threshold = clash_threshold(first.GetSymbol(), second.GetSymbol(), fraction)
            label = f"{first.GetSymbol()}{i + 1} — {second.GetSymbol()}{j + 1}"
            margin = distance - threshold
            if closest is None or margin < closest[0]:
                closest = (margin, label, distance, threshold)
            if distance < threshold:
                violations.append(
                    Measurement(
                        label=label,
                        value=round(distance, 3),
                        low=round(threshold, 3),
                        passed=False,
                        detail=(
                            f"Атомы разделены {int(topological[i][j])} связями, "
                            f"но сближены сильнее {fraction:g} суммы радиусов"
                        ),
                    )
                )

    if closest is None:
        return GeometryCheck(
            name="internal_clash",
            title=title,
            status=NOT_APPLICABLE,
            summary=(
                "В молекуле нет пар атомов, разделённых "
                f"{INTERNAL_CLASH_MIN_BOND_SEPARATION} и более связями."
            ),
            critical=True,
        )

    margin, label, distance, threshold = closest
    measurements = sorted(violations, key=lambda item: item.value or 0.0)
    if not violations:
        # Как и для контактов с белком: число показывается и при успехе,
        # но при нарушении оно уже стоит первой строкой и не дублируется.
        measurements.append(
            Measurement(
                label=f"Самая тесная несвязанная пара: {label}",
                value=round(distance, 3),
                low=round(threshold, 3),
                passed=True,
                detail=f"Запас до порога {margin:+.3f} Å",
            )
        )
    if violations:
        return GeometryCheck(
            name="internal_clash",
            title=title,
            status=FAILED,
            summary=f"Несвязанных пар ближе порога: {len(violations)}.",
            measurements=measurements,
            critical=True,
        )
    return GeometryCheck(
        name="internal_clash",
        title=title,
        status=PASSED,
        summary="Несвязанные атомы лиганда не сближены сверх нормы.",
        measurements=measurements,
        critical=True,
    )


def _bounds_matrix(mol: Any) -> Any | None:
    try:
        return rdDistGeom.GetMoleculeBoundsMatrix(mol)
    except Exception:  # noqa: BLE001 — границы не строятся для нестандартной валентности
        return None


def _bounds_for(bounds: Any, i: int, j: int) -> tuple[float, float]:
    """Нижняя и верхняя границы расстояния между атомами i и j.

    В матрице RDKit верхний треугольник хранит верхние границы, нижний —
    нижние, поэтому индексы нельзя брать в произвольном порядке.
    """
    low, high = (i, j) if i < j else (j, i)
    return float(bounds[high][low]), float(bounds[low][high])


def check_bond_lengths(mol: Any, *, tolerance: float = DISTANCE_GEOMETRY_TOLERANCE) -> GeometryCheck:
    """Длины связей лиганда против границ дистанционной геометрии."""
    title = "Длины связей"
    if mol is None or not _RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
        return GeometryCheck(
            name="bond_lengths",
            title=title,
            status=NOT_APPLICABLE,
            summary="Молекула позы не собрана — длины связей не считались.",
        )
    bounds = _bounds_matrix(mol)
    if bounds is None:
        return GeometryCheck(
            name="bond_lengths",
            title=title,
            status=NOT_APPLICABLE,
            summary="RDKit не построил границы дистанционной геометрии для этой молекулы.",
        )

    conf = mol.GetConformer()
    measurements: list[Measurement] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        low, high = _bounds_for(bounds, i, j)
        low *= 1 - tolerance
        high *= 1 + tolerance
        length = rdMolTransforms.GetBondLength(conf, i, j)
        ok = low <= length <= high
        label = (
            f"{mol.GetAtomWithIdx(i).GetSymbol()}{i + 1}–"
            f"{mol.GetAtomWithIdx(j).GetSymbol()}{j + 1}"
        )
        measurements.append(
            Measurement(
                label=label,
                value=round(length, 3),
                low=round(low, 3),
                high=round(high, 3),
                passed=ok,
                detail="" if ok else "Длина связи вне допустимого диапазона",
            )
        )

    bad = [item for item in measurements if not item.passed]
    if not measurements:
        return GeometryCheck(
            name="bond_lengths",
            title=title,
            status=NOT_APPLICABLE,
            summary="У молекулы нет связей для проверки.",
        )
    if bad:
        return GeometryCheck(
            name="bond_lengths",
            title=title,
            status=FAILED,
            summary=f"Вне диапазона {len(bad)} связей из {len(measurements)}.",
            measurements=bad + [item for item in measurements if item.passed],
        )
    return GeometryCheck(
        name="bond_lengths",
        title=title,
        status=PASSED,
        summary=f"Все {len(measurements)} связей в допустимых диапазонах.",
        measurements=measurements,
    )


def check_bond_angles(mol: Any, *, tolerance: float = DISTANCE_GEOMETRY_TOLERANCE) -> GeometryCheck:
    """Валентные углы лиганда против границ дистанционной геометрии.

    Границы RDKit заданы на расстояниях, а не на углах, поэтому проверка
    идёт по расстоянию между крайними атомами, а в отчёт пишется угол:
    по нему нарушение читается человеком, а не только программой.
    """
    title = "Валентные углы"
    if mol is None or not _RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
        return GeometryCheck(
            name="bond_angles",
            title=title,
            status=NOT_APPLICABLE,
            summary="Молекула позы не собрана — углы не считались.",
        )
    bounds = _bounds_matrix(mol)
    if bounds is None:
        return GeometryCheck(
            name="bond_angles",
            title=title,
            status=NOT_APPLICABLE,
            summary="RDKit не построил границы дистанционной геометрии для этой молекулы.",
        )

    conf = mol.GetConformer()
    measurements: list[Measurement] = []
    for atom in mol.GetAtoms():
        neighbors = [nb.GetIdx() for nb in atom.GetNeighbors()]
        if len(neighbors) < 2:
            continue
        j = atom.GetIdx()
        for first in range(len(neighbors)):
            for second in range(first + 1, len(neighbors)):
                i, k = neighbors[first], neighbors[second]
                side_a = rdMolTransforms.GetBondLength(conf, i, j)
                side_b = rdMolTransforms.GetBondLength(conf, j, k)
                span = rdMolTransforms.GetBondLength(conf, i, k)
                low, high = _bounds_for(bounds, i, k)
                low *= 1 - tolerance
                high *= 1 + tolerance
                angle = _angle_deg(side_a, side_b, span)
                # Большему расстоянию между крайними атомами отвечает больший
                # угол, поэтому верхняя граница расстояния даёт верхнюю границу
                # угла. Расстояние вне треугольника невозможно, отсюда зажимы.
                angle_max = _angle_deg(side_a, side_b, min(high, side_a + side_b))
                angle_min = _angle_deg(side_a, side_b, max(low, abs(side_a - side_b)))
                ok = low <= span <= high
                label = (
                    f"{mol.GetAtomWithIdx(i).GetSymbol()}{i + 1}–"
                    f"{atom.GetSymbol()}{j + 1}–"
                    f"{mol.GetAtomWithIdx(k).GetSymbol()}{k + 1}"
                )
                measurements.append(
                    Measurement(
                        label=label,
                        value=round(angle, 1) if angle is not None else None,
                        unit="°",
                        low=round(angle_min, 1) if angle_min is not None else None,
                        high=round(angle_max, 1) if angle_max is not None else None,
                        passed=ok,
                        detail=(
                            ""
                            if ok
                            else f"Расстояние между крайними атомами {span:.3f} Å вне "
                            f"диапазона {low:.3f}–{high:.3f} Å"
                        ),
                    )
                )

    if not measurements:
        return GeometryCheck(
            name="bond_angles",
            title=title,
            status=NOT_APPLICABLE,
            summary="У молекулы нет валентных углов для проверки.",
        )
    bad = [item for item in measurements if not item.passed]
    if bad:
        return GeometryCheck(
            name="bond_angles",
            title=title,
            status=FAILED,
            summary=f"Вне диапазона {len(bad)} углов из {len(measurements)}.",
            measurements=bad + [item for item in measurements if item.passed],
        )
    return GeometryCheck(
        name="bond_angles",
        title=title,
        status=PASSED,
        summary=f"Все {len(measurements)} валентных углов в допустимых диапазонах.",
        measurements=measurements,
    )


def check_stereochemistry(mol: Any, smiles: str | None) -> GeometryCheck:
    """Сохранилась ли заданная конфигурация хиральных центров и двойных связей.

    Проверяется только то, что было задано: если в SMILES прогона стереохимия
    не указана, проверять нечего и результат — «не применимо», а не провал.
    """
    title = "Стереохимия"
    if mol is None or not _RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=NOT_APPLICABLE,
            summary="Молекула позы не собрана — стереохимию не проверяли.",
        )
    if not smiles:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=NOT_APPLICABLE,
            summary="SMILES прогона не передан: заданной стереохимии нет — не применимо.",
        )
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=NOT_APPLICABLE,
            summary=f"SMILES прогона не разобран ({smiles[:40]}): стереохимию не проверяли.",
        )
    template = Chem.RemoveHs(template)
    Chem.AssignStereochemistry(template, cleanIt=True, force=True)

    expected_centers = {
        index: code
        for index, code in Chem.FindMolChiralCenters(
            template, includeUnassigned=False, useLegacyImplementation=False
        )
    }
    expected_bonds = {
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()): str(bond.GetStereo())
        for bond in template.GetBonds()
        if bond.GetStereo() != Chem.BondStereo.STEREONONE
    }
    if not expected_centers and not expected_bonds:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=NOT_APPLICABLE,
            summary="В молекуле не задана стереохимия — не применимо.",
        )

    match = mol.GetSubstructMatch(template)
    if not match:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=NOT_APPLICABLE,
            summary="Поза не сопоставилась с шаблоном SMILES: стереохимию не проверяли.",
        )

    pose = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(pose)
    actual_centers = {
        index: code
        for index, code in Chem.FindMolChiralCenters(
            pose, includeUnassigned=True, useLegacyImplementation=False
        )
    }

    measurements: list[Measurement] = []
    for template_index, expected in expected_centers.items():
        pose_index = match[template_index]
        actual = actual_centers.get(pose_index, "не определена")
        ok = actual == expected
        measurements.append(
            Measurement(
                label=f"Хиральный центр {pose.GetAtomWithIdx(pose_index).GetSymbol()}{pose_index + 1}",
                value=actual,
                unit="",
                passed=ok,
                detail=f"Задано {expected}" + ("" if ok else f", в позе {actual}"),
            )
        )

    for (begin, end), expected in expected_bonds.items():
        pose_bond = pose.GetBondBetweenAtoms(match[begin], match[end])
        actual = str(pose_bond.GetStereo()) if pose_bond is not None else "нет связи"
        # STEREOCIS/STEREOZ и STEREOTRANS/STEREOE обозначают одно и то же.
        ok = actual == expected or {actual[-1], expected[-1]} in ({"Z", "S"}, {"E", "S"})
        measurements.append(
            Measurement(
                label=f"Двойная связь {match[begin] + 1}={match[end] + 1}",
                value=actual,
                unit="",
                passed=ok,
                detail=f"Задано {expected}",
            )
        )

    bad = [item for item in measurements if not item.passed]
    if bad:
        return GeometryCheck(
            name="stereochemistry",
            title=title,
            status=FAILED,
            summary=f"Изменилась конфигурация: {len(bad)} из {len(measurements)}.",
            measurements=bad + [item for item in measurements if item.passed],
            critical=True,
        )
    return GeometryCheck(
        name="stereochemistry",
        title=title,
        status=PASSED,
        summary=f"Заданная стереохимия сохранена ({len(measurements)} элементов).",
        measurements=measurements,
        critical=True,
    )


def check_ring_planarity(mol: Any, *, tolerance: float = RING_PLANARITY_TOLERANCE_A) -> GeometryCheck:
    """Сохранилась ли планарность ароматических колец.

    Кольца берутся из молекулы, а не из графа расстояний: ароматичность — это
    свойство химии, а не геометрии, и определять её по тому же расстоянию,
    которое потом проверяется, значило бы проверять круговую величину.
    """
    title = "Планарность колец"
    if mol is None or not _RDKIT_AVAILABLE or mol.GetNumConformers() == 0:
        return GeometryCheck(
            name="ring_planarity",
            title=title,
            status=NOT_APPLICABLE,
            summary="Молекула позы не собрана — планарность не считалась.",
        )
    conf = mol.GetConformer()
    measurements: list[Measurement] = []
    for ring in mol.GetRingInfo().AtomRings():
        if not all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            continue
        points = [
            (conf.GetAtomPosition(index).x, conf.GetAtomPosition(index).y, conf.GetAtomPosition(index).z)
            for index in ring
        ]
        deviation = _plane_deviation(points)
        label = "Кольцо " + "–".join(
            f"{mol.GetAtomWithIdx(index).GetSymbol()}{index + 1}" for index in ring
        )
        measurements.append(
            Measurement(
                label=label,
                value=round(deviation, 3),
                high=tolerance,
                passed=deviation <= tolerance,
                detail="" if deviation <= tolerance else "Кольцо деформировано",
            )
        )

    if not measurements:
        return GeometryCheck(
            name="ring_planarity",
            title=title,
            status=NOT_APPLICABLE,
            summary="В молекуле нет ароматических колец — не применимо.",
        )
    bad = [item for item in measurements if not item.passed]
    if bad:
        return GeometryCheck(
            name="ring_planarity",
            title=title,
            status=FAILED,
            summary=f"Деформировано колец: {len(bad)} из {len(measurements)}.",
            measurements=bad + [item for item in measurements if item.passed],
        )
    return GeometryCheck(
        name="ring_planarity",
        title=title,
        status=PASSED,
        summary=f"Все {len(measurements)} ароматических колец плоские.",
        measurements=measurements,
    )


def check_metal_coordination(
    ligand: list[Atom],
    protein: list[Atom],
    *,
    reference_ligand: list[Atom] | None = None,
    max_distance: float = METAL_COORDINATION_MAX_A,
    tolerance: float = METAL_REFERENCE_TOLERANCE_A,
) -> GeometryCheck:
    """Геометрия контакта лиганда с ионом металла, если металл в структуре есть.

    Эталон не задаётся числом в коде: если передан кристаллический лиганд,
    то же расстояние измеряется на нём, и в отчёт идёт отклонение.
    """
    title = "Координация металла"
    metals = [atom for atom in protein if atom.is_metal]
    if not metals:
        return GeometryCheck(
            name="metal_coordination",
            title=title,
            status=NOT_APPLICABLE,
            summary="В структуре нет ионов металла — не применимо.",
        )
    donors = [
        atom
        for atom in heavy_atoms(ligand)
        if atom.element.upper() in _METAL_DONOR_ELEMENTS
    ]
    if not donors:
        return GeometryCheck(
            name="metal_coordination",
            title=title,
            status=NOT_APPLICABLE,
            summary="В лиганде нет атомов, способных координировать металл — не применимо.",
        )

    def _closest(candidates: list[Atom]) -> tuple[Atom, Atom, float] | None:
        best: tuple[Atom, Atom, float] | None = None
        for donor in candidates:
            for metal in metals:
                distance = donor.distance_to(metal)
                if best is None or distance < best[2]:
                    best = (donor, metal, distance)
        return best

    closest = _closest(donors)
    if closest is None:  # недостижимо: доноры и металлы выше уже непусты
        return GeometryCheck(
            name="metal_coordination",
            title=title,
            status=NOT_APPLICABLE,
            summary="Пара «донор — металл» не найдена.",
        )
    donor, metal, distance = closest
    metal_label = metal.element.upper()

    if distance > max_distance:
        # Лиганд до металла не достаёт. Это не провал геометрии, а отсутствие
        # координации: число сохраняется, чтобы человек мог судить сам.
        return GeometryCheck(
            name="metal_coordination",
            title=title,
            status=NOT_APPLICABLE,
            summary=(
                f"Лиганд не координирует {metal_label}: ближайший "
                f"{donor.element}{donor.serial} в {distance:.2f} Å, "
                f"порог координации {max_distance:g} Å."
            ),
            measurements=[
                Measurement(
                    label=f"{donor.element}{donor.serial}–{metal_label}",
                    value=round(distance, 3),
                    high=max_distance,
                    passed=True,
                    detail="Координация не обнаружена",
                )
            ],
        )

    measurements = [
        Measurement(
            label=f"{donor.element}{donor.serial}–{metal_label} ({metal.residue_label})",
            value=round(distance, 3),
            high=max_distance,
            passed=True,
            detail="Расстояние координации",
        )
    ]

    status = PASSED
    summary = f"{donor.element}–{metal_label} = {distance:.2f} Å."
    if reference_ligand:
        reference_donors = [
            atom
            for atom in heavy_atoms(reference_ligand)
            if atom.element.upper() in _METAL_DONOR_ELEMENTS
        ]
        reference = _closest(reference_donors) if reference_donors else None
        if reference is not None:
            deviation = distance - reference[2]
            ok = abs(deviation) <= tolerance
            measurements.append(
                Measurement(
                    label=f"Отклонение от кристаллического {metal_label}-контакта",
                    value=round(deviation, 3),
                    low=-tolerance,
                    high=tolerance,
                    passed=ok,
                    detail=f"В кристалле {reference[2]:.2f} Å",
                )
            )
            status = PASSED if ok else FAILED
            summary += f" В кристалле {reference[2]:.2f} Å, отклонение {deviation:+.2f} Å."
        else:
            summary += " Эталонный контакт в кристалле не найден."

    return GeometryCheck(
        name="metal_coordination",
        title=title,
        status=status,
        summary=summary,
        measurements=measurements,
    )
