"""RMSD между позами лиганда с учётом симметрии молекулы.

До этого модуля RMSD в проекте не вычислялся нигде: числа в poses.csv только
читались. Значит, значение 0,356 Å для 3PTB ничем не производилось и проверить
его было нельзя.

Два правила, от которых зависит смысл результата:
1. Координаты не выравниваются. Поза и кристалл уже стоят в одной системе
   координат белка, и наложение структур превратило бы «попали ли мы в карман»
   в «похожа ли форма молекулы» — а это совсем другой вопрос.
2. Симметрия учитывается. У бензамидина два азота амидиновой группы
   неразличимы, и нумерация атомов в докированной позе не обязана совпадать с
   кристаллической. Без перебора автоморфизмов правильная поза получила бы
   завышенный RMSD просто из-за порядка атомов.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

# Элементы, которые не участвуют в RMSD тяжёлых атомов.
_IGNORED_ELEMENTS = {"H", "D"}
# Псевдоатомы AutoDock, не являющиеся настоящими атомами лиганда.
_IGNORED_AD_TYPES = {"G0", "G1", "G2", "G3"}


@dataclass(frozen=True)
class Atom:
    element: str
    x: float
    y: float
    z: float


class RmsdError(RuntimeError):
    """RMSD посчитать нельзя. Поднимается вместо возврата правдоподобного числа."""


# Типы атомов AutoDock → химический элемент. Тип — не элемент: «A» это
# ароматический углерод, «NA» — акцепторный азот, «OA» — акцепторный кислород.
_AD_TO_ELEMENT = {
    "A": "C", "C": "C", "CG0": "C", "CG1": "C", "CG2": "C", "CG3": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P", "F": "F", "I": "I",
    "CL": "CL", "BR": "BR",
    "H": "H", "HD": "H", "HS": "H",
    "ZN": "ZN", "MG": "MG", "MN": "MN", "FE": "FE", "CU": "CU", "CA": "CA",
    "CO": "CO", "NI": "NI", "K": "K",
}


def _element_from_line(line: str) -> str:
    """Элемент атома из строки PDB или PDBQT.

    В PDBQT колонки 77–78 заняты типом AutoDock, а не элементом. Неизвестный
    тип не угадывается: пустая строка отбрасывает атом, и несовпадение состава
    потом поднимет явную ошибку вместо неверного числа.
    """
    # Тип AutoDock в PDBQT стоит в колонках 78–79, а элемент в PDB — в 77–78.
    # Срез [76:78] отрезал вторую букву типа: «Cl» читался как углерод «C», а
    # «Br» как бор, и RMSD любого хлорсодержащего лиганда (гефитиниб) не считался.
    # Первый токен хвоста строки одинаково верен для обоих форматов.
    tokens = line[76:].split() if len(line) > 76 else []
    tail = tokens[0].upper() if tokens else ""
    if tail in _AD_TO_ELEMENT:
        return _AD_TO_ELEMENT[tail]
    if tail:
        return tail
    name = line[12:16].strip().upper()
    if name[:2] in ("CL", "BR"):
        return name[:2]
    return name[:1] if name else ""


def parse_ligand_atoms(text: str, resname: str | None = None) -> list[Atom]:
    """Тяжёлые атомы лиганда из PDB/PDBQT. Только первая модель."""
    atoms: list[Atom] = []
    for raw in text.splitlines():
        if raw.startswith("ENDMDL"):
            break
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        if resname is not None and raw[17:20].strip().upper() != resname.upper():
            continue
        element = _element_from_line(raw)
        if not element or element in _IGNORED_ELEMENTS:
            continue
        try:
            atoms.append(Atom(element, float(raw[30:38]), float(raw[38:46]), float(raw[46:54])))
        except ValueError:
            continue
    return atoms


def read_pose_atoms(path: Path | str) -> list[Atom]:
    return parse_ligand_atoms(Path(path).read_text(encoding="utf-8", errors="replace"))


_METAL_CODES = {"ZN", "MG", "MN", "FE", "CU", "CO", "NI", "CA"}


def _copies(pdb_text: str, comp_id: str, chain: str | None) -> dict[tuple[str, str], list[Atom]]:
    """Кристаллические копии лиганда, сгруппированные по (цепь, номер остатка).

    В одной структуре лиганд часто лежит в нескольких местах: в 3HS4 три копии
    ацетазоламида, из них только одна в активном центре. Склеивать их в один
    эталон нельзя — получится «молекула» из 39 атомов вместо 13.
    """
    groups: dict[tuple[str, str], list[Atom]] = {}
    wanted = comp_id.upper()
    for raw in pdb_text.splitlines():
        if raw.startswith("ENDMDL"):
            break
        if not raw.startswith(("ATOM", "HETATM")) or raw[17:20].strip().upper() != wanted:
            continue
        if chain and raw[21:22].strip().upper() != chain.strip().upper():
            continue
        element = _element_from_line(raw)
        if not element or element in _IGNORED_ELEMENTS:
            continue
        try:
            atom = Atom(element, float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))
        except ValueError:
            continue
        key = (raw[21:22].strip(), raw[22:27].strip())
        groups.setdefault(key, []).append(atom)
    return groups


def _metal_positions(pdb_text: str, chain: str | None) -> list[Atom]:
    metals: list[Atom] = []
    for raw in pdb_text.splitlines():
        if not raw.startswith("HETATM") or raw[17:20].strip().upper() not in _METAL_CODES:
            continue
        if chain and raw[21:22].strip().upper() != chain.strip().upper():
            continue
        try:
            metals.append(
                Atom(raw[17:20].strip().upper(), float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))
            )
        except ValueError:
            continue
    return metals


def select_crystal_reference(
    pdb_text: str, comp_id: str, chain: str | None = None
) -> tuple[list[Atom], str]:
    """Эталонная копия лиганда и пояснение, почему выбрана именно она.

    Правило не смотрит на докированную позу: иначе эталоном становилась бы
    копия, на которую поза похожа сильнее всего, и метрика подгонялась бы под
    результат. Выбор: копия ближе всего к иону металла (активный центр
    металлофермента); без металла — первая копия цепи по номеру остатка.
    """
    copies = _copies(pdb_text, comp_id, chain)
    if not copies:
        raise RmsdError(f"В структуре нет лиганда {comp_id}" + (f" в цепи {chain}" if chain else ""))
    if len(copies) == 1:
        (key, atoms), = copies.items()
        return atoms, f"единственная копия {comp_id} {key[0]}{key[1]}"

    metals = _metal_positions(pdb_text, chain)
    if metals:
        def distance(atoms: list[Atom]) -> float:
            return min(
                math.sqrt((m.x - a.x) ** 2 + (m.y - a.y) ** 2 + (m.z - a.z) ** 2)
                for m in metals
                for a in atoms
            )

        key = min(copies, key=lambda k: distance(copies[k]))
        return copies[key], (
            f"копия {comp_id} {key[0]}{key[1]}, ближайшая к иону металла "
            f"({distance(copies[key]):.2f} Å); всего копий {len(copies)}"
        )

    def order(k: tuple[str, str]) -> tuple[str, int]:
        digits = "".join(ch for ch in k[1] if ch.isdigit() or ch == "-")
        return (k[0], int(digits) if digits not in ("", "-") else 0)

    key = min(copies, key=order)
    return copies[key], f"первая копия {comp_id} {key[0]}{key[1]} из {len(copies)} (металла в структуре нет)"


def crystal_ligand_atoms(pdb_text: str, comp_id: str, chain: str | None = None) -> list[Atom]:
    """Атомы эталонной копии кристаллического лиганда."""
    return select_crystal_reference(pdb_text, comp_id, chain)[0]


def _plain_rmsd(a: list[Atom], b: list[Atom]) -> float:
    total = 0.0
    for p, q in zip(a, b):
        total += (p.x - q.x) ** 2 + (p.y - q.y) ** 2 + (p.z - q.z) ** 2
    return math.sqrt(total / len(a))


def _composition(atoms: list[Atom]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for atom in atoms:
        counts[atom.element] = counts.get(atom.element, 0) + 1
    return counts


def _rdkit_symmetry_rmsd(probe_block: str, ref_block: str) -> float | None:
    """RMSD по лучшему из автоморфизмов молекулы. None, если RDKit не справился.

    Связи достраиваются по расстояниям, а их порядок приводится к одинарному:
    сравниваются геометрия и топология, а не способ записи ароматичности,
    который у докированной позы и у кристалла различается.
    """
    try:
        from rdkit import Chem, RDLogger
    except ImportError:
        return None

    RDLogger.DisableLog("rdApp.*")

    def _mol(block: str):
        mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True, proximityBonding=True)
        if mol is None:
            return None
        editable = Chem.RWMol(mol)
        for bond in editable.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
        for atom in editable.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
        result = editable.GetMol()
        try:
            Chem.SanitizeMol(result, Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        except Exception:  # noqa: BLE001 — санитизация здесь необязательна
            pass
        return result

    probe, ref = _mol(probe_block), _mol(ref_block)
    if probe is None or ref is None:
        return None
    if probe.GetNumAtoms() != ref.GetNumAtoms() or probe.GetNumAtoms() == 0:
        return None

    matches = probe.GetSubstructMatches(ref, uniquify=False, useChirality=False, maxMatches=5000)
    if not matches:
        return None

    probe_conf, ref_conf = probe.GetConformer(), ref.GetConformer()
    best: float | None = None
    for match in matches:
        total = 0.0
        for ref_idx, probe_idx in enumerate(match):
            p = probe_conf.GetAtomPosition(probe_idx)
            q = ref_conf.GetAtomPosition(ref_idx)
            total += (p.x - q.x) ** 2 + (p.y - q.y) ** 2 + (p.z - q.z) ** 2
        value = math.sqrt(total / len(match))
        if best is None or value < best:
            best = value
    return best


def _atoms_to_pdb_block(atoms: list[Atom], resname: str = "LIG") -> str:
    lines = []
    for index, atom in enumerate(atoms, start=1):
        name = f"{atom.element}{index}"[:4]
        lines.append(
            f"HETATM{index:5d} {name:<4} {resname:>3} A   1    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {atom.element:>2}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def symmetry_rmsd(probe: list[Atom], reference: list[Atom]) -> float:
    """RMSD без выравнивания, с перебором симметрично эквивалентных нумераций."""
    if not probe or not reference:
        raise RmsdError("Нет атомов для расчёта RMSD")
    if len(probe) != len(reference):
        raise RmsdError(
            f"Разное число тяжёлых атомов: поза {len(probe)}, эталон {len(reference)}. "
            "RMSD между разными молекулами не определён."
        )
    if _composition(probe) != _composition(reference):
        raise RmsdError(
            f"Разный состав: поза {_composition(probe)}, эталон {_composition(reference)}."
        )

    value = _rdkit_symmetry_rmsd(_atoms_to_pdb_block(probe), _atoms_to_pdb_block(reference))
    if value is not None:
        return value
    # Без RDKit остаётся сравнение в исходном порядке атомов. Оно корректно
    # только если нумерация совпадает, поэтому результат может быть завышен —
    # но не занижен, то есть порог «≤ 2 Å» так не обойти.
    return _plain_rmsd(probe, reference)


def pose_rmsd(probe_path: Path | str, reference_path: Path | str) -> float:
    """RMSD между двумя файлами поз одного лиганда."""
    return symmetry_rmsd(read_pose_atoms(probe_path), read_pose_atoms(reference_path))


def rmsd_to_crystal(
    pose_path: Path | str, pdb_text: str, comp_id: str, chain: str | None = None
) -> float:
    """RMSD позы до эталонной копии кристаллического лиганда в системе координат белка."""
    return symmetry_rmsd(read_pose_atoms(pose_path), crystal_ligand_atoms(pdb_text, comp_id, chain))


def centroid(atoms: list[Atom]) -> tuple[float, float, float]:
    n = len(atoms)
    if not n:
        raise RmsdError("Пустой список атомов")
    return (
        sum(a.x for a in atoms) / n,
        sum(a.y for a in atoms) / n,
        sum(a.z for a in atoms) / n,
    )


def min_distance_to_atoms(point: tuple[float, float, float], atoms: list[Atom]) -> float:
    """Расстояние от точки до ближайшего тяжёлого атома.

    По этой величине проверяется контрольный критерий защиты: карман считается
    верным, если его центр лежит в пределах 4 Å от кристаллического лиганда.
    """
    if not atoms:
        raise RmsdError("Пустой список атомов")
    px, py, pz = point
    return min(
        math.sqrt((px - a.x) ** 2 + (py - a.y) ** 2 + (pz - a.z) ** 2) for a in atoms
    )
