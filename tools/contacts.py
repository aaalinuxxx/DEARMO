"""Контакты белок–лиганд по геометрии и свойствам атомов.

Подход PLIP, но без внешней зависимости: типы взаимодействий определяются
по элементам, принадлежности к заряженной группе и к ароматическому кольцу,
а решение принимается по расстоянию между тяжёлыми атомами.

Каждый контакт несёт имена и номера обоих атомов и измеренное расстояние,
чтобы любое число с графика сценария 3 можно было сверить по координатам.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

# --- Параметры критериев. Собраны в одном месте: их спрашивают на защите. ---

HYDROPHOBIC_MAX_A = 4.0
HBOND_MAX_A = 3.5
SALT_BRIDGE_MAX_A = 4.0
METAL_MAX_A = 2.8
HALOGEN_MAX_A = 3.5
PI_STACK_MAX_A = 5.5
PI_STACK_PARALLEL_DEG = 30.0
PI_STACK_TSHAPE_MIN_DEG = 50.0
PI_CATION_MAX_A = 6.0
RING_PLANARITY_MAX_A = 0.25

INTERACTION_TYPES = (
    "metal",
    "salt_bridge",
    "hydrogen_bond",
    "halogen_bond",
    "pi_stacking",
    "pi_cation",
    "hydrophobic",
)

# Приоритет при совпадении одной пары атомов по нескольким критериям:
# специфичное взаимодействие вытесняет неспецифичное.
_PRIORITY = {name: index for index, name in enumerate(INTERACTION_TYPES)}

_METALS = {"ZN", "MG", "MN", "FE", "CU", "CO", "NI", "CA", "NA", "K", "CD", "HG"}
_HALOGENS = {"F", "CL", "BR", "I"}
_POLAR = {"N", "O", "S"}

_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SE": 1.20,
    "ZN": 1.22, "MG": 1.41, "MN": 1.39, "FE": 1.32, "CU": 1.32,
    "CO": 1.26, "NI": 1.24, "CA": 1.76, "NA": 1.66, "K": 2.03,
}
_BOND_TOLERANCE_A = 0.45

# Типы AutoDock → элемент. В PDBQT колонка элемента занята типом AD,
# поэтому без этой карты металл читается как углерод.
_AD_TO_ELEMENT = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "O": "O", "OA": "O",
    "OS": "O", "S": "S", "SA": "S", "H": "H", "HD": "H", "HS": "H",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I", "P": "P", "SE": "SE",
    "ZN": "ZN", "MG": "MG", "MN": "MN", "FE": "FE", "CU": "CU",
    "CO": "CO", "NI": "NI", "CA": "CA",
}

# Заряженные группы белка: остаток → (атомы, знак формального заряда).
_PROTEIN_CHARGED = {
    "ASP": (("OD1", "OD2"), -1),
    "GLU": (("OE1", "OE2"), -1),
    "ARG": (("NE", "NH1", "NH2"), +1),
    "LYS": (("NZ",), +1),
    "HIS": (("ND1", "NE2"), +1),
}

# Ароматические кольца белка по именам атомов остатка.
_PROTEIN_RINGS: dict[str, tuple[tuple[str, ...], ...]] = {
    "PHE": (("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),),
    "TYR": (("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),),
    "HIS": (("CG", "ND1", "CD2", "CE1", "NE2"),),
    "TRP": (
        ("CG", "CD1", "NE1", "CE2", "CD2"),
        ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    ),
}

# Доноры водородной связи в белке по имени атома (остальные N/O считаем акцепторами).
_BACKBONE_DONOR = {"N"}
_BACKBONE_ACCEPTOR = {"O", "OXT"}


@dataclass(frozen=True)
class Atom:
    """Один атом из PDB или PDBQT. Координаты в Å."""

    serial: int
    name: str
    element: str
    res_name: str
    res_seq: int
    chain: str
    x: float
    y: float
    z: float
    is_hetatm: bool = False

    @property
    def residue_label(self) -> str:
        """Подпись вида Asp189 — так остаток называется в отчёте и на графике."""
        name = self.res_name.capitalize() if self.res_name.isalpha() else self.res_name
        return f"{name}{self.res_seq}"

    @property
    def is_metal(self) -> bool:
        return self.element.upper() in _METALS

    @property
    def is_heavy(self) -> bool:
        return self.element.upper() != "H"

    def distance_to(self, other: "Atom") -> float:
        return math.dist((self.x, self.y, self.z), (other.x, other.y, other.z))


@dataclass(frozen=True)
class Contact:
    """Одно взаимодействие лиганда с белком."""

    interaction_type: str
    distance_a: float
    ligand_atom: str
    ligand_atom_serial: int
    ligand_element: str
    protein_atom: str
    protein_atom_serial: int
    residue_name: str
    residue_seq: int
    chain: str

    @property
    def residue_label(self) -> str:
        name = self.residue_name.capitalize() if self.residue_name.isalpha() else self.residue_name
        return f"{name}{self.residue_seq}"

    @property
    def pair_key(self) -> tuple[int, int]:
        """Пара атомов — то, что измеряется одинаково во всех позах."""
        return (self.ligand_atom_serial, self.protein_atom_serial)


@dataclass
class Ring:
    atoms: list[Atom]
    centroid: tuple[float, float, float]
    normal: tuple[float, float, float]
    label: str = ""


# --------------------------------------------------------------------------
# Чтение координат
# --------------------------------------------------------------------------


def _element_from_name(raw_name: str) -> str:
    name = raw_name.strip()
    if not name:
        return "C"
    head = name.lstrip("0123456789")
    if len(head) >= 2 and head[:2].upper() in _COVALENT_RADII:
        return head[:2].upper()
    return head[:1].upper() if head else "C"


def _parse_atom_line(line: str, *, pdbqt: bool) -> Atom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except (ValueError, IndexError):
        return None
    try:
        serial = int(line[6:11])
    except ValueError:
        serial = 0
    name = line[12:16].strip()
    res_name = line[17:20].strip()
    chain = line[21:22].strip() or "A"
    raw_seq = line[22:26].strip()
    try:
        res_seq = int(raw_seq)
    except ValueError:
        res_seq = 0

    if pdbqt:
        ad_type = line[77:79].strip().upper() if len(line) >= 78 else ""
        element = _AD_TO_ELEMENT.get(ad_type, "")
        if not element:
            element = _element_from_name(name)
    else:
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = _element_from_name(name)
    return Atom(
        serial=serial,
        name=name,
        element=element.upper(),
        res_name=res_name,
        res_seq=res_seq,
        chain=chain,
        x=x,
        y=y,
        z=z,
        is_hetatm=line.startswith("HETATM"),
    )


def parse_atoms(text: str, *, pdbqt: bool, model: int | None = None) -> list[Atom]:
    """Атомы из текста PDB/PDBQT. `model` выбирает MODEL N в многомодельном файле."""
    atoms: list[Atom] = []
    current_model = 0
    inside_wanted = model is None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current_model += 1
            parts = line.split()
            if len(parts) > 1:
                try:
                    current_model = int(parts[1])
                except ValueError:
                    pass
            inside_wanted = model is None or current_model == model
            continue
        if line.startswith("ENDMDL"):
            inside_wanted = model is None
            continue
        if not inside_wanted:
            continue
        atom = _parse_atom_line(line, pdbqt=pdbqt)
        if atom is not None:
            atoms.append(atom)
    return atoms


def count_models(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.startswith("MODEL"))


def read_structure(path: Path | str, *, model: int | None = None) -> list[Atom]:
    """Атомы из файла. Формат определяется по расширению, модель — по номеру."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    is_pdbqt = path.suffix.lower() == ".pdbqt"
    if model is None and is_pdbqt and count_models(text) > 1:
        model = 1
    return parse_atoms(text, pdbqt=is_pdbqt, model=model)


def heavy_atoms(atoms: list[Atom]) -> list[Atom]:
    return [atom for atom in atoms if atom.is_heavy]


# --------------------------------------------------------------------------
# Свойства атомов: связи, заряженные группы, кольца
# --------------------------------------------------------------------------


def bond_graph(atoms: list[Atom]) -> dict[int, set[int]]:
    """Связи по сумме ковалентных радиусов. Индексы — позиции в списке."""
    graph: dict[int, set[int]] = {index: set() for index in range(len(atoms))}
    for i in range(len(atoms)):
        ri = _COVALENT_RADII.get(atoms[i].element, 0.77)
        for j in range(i + 1, len(atoms)):
            rj = _COVALENT_RADII.get(atoms[j].element, 0.77)
            dist = atoms[i].distance_to(atoms[j])
            if 0.4 < dist <= ri + rj + _BOND_TOLERANCE_A:
                graph[i].add(j)
                graph[j].add(i)
    return graph


def ligand_charged_atoms(atoms: list[Atom]) -> dict[int, int]:
    """Формальный заряд групп лиганда, распознанных по графу связей.

    Распознаются две группы, которых хватает для демонстрационных пар:
    амидин/гуанидин (углерод с двумя и более азотами) — положительная,
    карбоксилат (углерод ровно с двумя кислородами) — отрицательная.
    Остальные атомы считаются незаряженными: лучше пропустить солевой мостик,
    чем объявить несуществующий.
    """
    heavy = [atom for atom in atoms if atom.is_heavy]
    graph = bond_graph(heavy)
    charges: dict[int, int] = {}
    for index, atom in enumerate(heavy):
        if atom.element != "C":
            continue
        nitrogens = [n for n in graph[index] if heavy[n].element == "N"]
        oxygens = [n for n in graph[index] if heavy[n].element == "O"]
        if len(nitrogens) >= 2:
            for n in nitrogens:
                if len(graph[n]) <= 2:
                    charges[heavy[n].serial] = +1
        if len(oxygens) == 2 and all(len(graph[o]) == 1 for o in oxygens):
            for o in oxygens:
                charges[heavy[o].serial] = -1
    return charges


def protein_charge(atom: Atom) -> int:
    entry = _PROTEIN_CHARGED.get(atom.res_name.upper())
    if not entry:
        return 0
    names, sign = entry
    return sign if atom.name.upper() in names else 0


def _ring_plane(atoms: list[Atom]) -> tuple[tuple[float, float, float], tuple[float, float, float], float]:
    """Центроид, нормаль (метод Ньюэлла) и максимальное отклонение от плоскости."""
    count = len(atoms)
    cx = sum(a.x for a in atoms) / count
    cy = sum(a.y for a in atoms) / count
    cz = sum(a.z for a in atoms) / count
    nx = ny = nz = 0.0
    for index in range(count):
        a = atoms[index]
        b = atoms[(index + 1) % count]
        nx += (a.y - b.y) * (a.z + b.z)
        ny += (a.z - b.z) * (a.x + b.x)
        nz += (a.x - b.x) * (a.y + b.y)
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-9:
        return (cx, cy, cz), (0.0, 0.0, 1.0), float("inf")
    normal = (nx / norm, ny / norm, nz / norm)
    deviation = max(
        abs((a.x - cx) * normal[0] + (a.y - cy) * normal[1] + (a.z - cz) * normal[2])
        for a in atoms
    )
    return (cx, cy, cz), normal, deviation


def _find_cycles(graph: dict[int, set[int]], sizes: tuple[int, ...] = (5, 6)) -> list[list[int]]:
    """Простые циклы длиной 5-6. Лиганды малы, полный обход допустим."""
    found: set[tuple[int, ...]] = set()
    cycles: list[list[int]] = []
    max_size = max(sizes)

    def walk(start: int, current: int, path: list[int], visited: set[int]) -> None:
        for neighbour in graph[current]:
            if neighbour == start and len(path) in sizes:
                key = tuple(sorted(path))
                if key not in found:
                    found.add(key)
                    cycles.append(list(path))
                continue
            if neighbour in visited or neighbour < start or len(path) >= max_size:
                continue
            visited.add(neighbour)
            path.append(neighbour)
            walk(start, neighbour, path, visited)
            path.pop()
            visited.discard(neighbour)

    for start in sorted(graph):
        walk(start, start, [start], {start})
    return cycles


def ligand_rings(atoms: list[Atom]) -> list[Ring]:
    """Плоские кольца лиганда из графа связей — кандидаты в ароматические."""
    heavy = [atom for atom in atoms if atom.is_heavy]
    if len(heavy) < 5:
        return []
    graph = bond_graph(heavy)
    rings: list[Ring] = []
    for cycle in _find_cycles(graph):
        members = [heavy[index] for index in cycle]
        if any(atom.element not in ("C", "N", "O", "S") for atom in members):
            continue
        centroid, normal, deviation = _ring_plane(members)
        if deviation > RING_PLANARITY_MAX_A:
            continue
        rings.append(Ring(atoms=members, centroid=centroid, normal=normal, label="ligand_ring"))
    return rings


def protein_rings(atoms: list[Atom]) -> list[Ring]:
    """Ароматические кольца белка по именам атомов остатков."""
    by_residue: dict[tuple[str, int, str], dict[str, Atom]] = {}
    for atom in atoms:
        key = (atom.chain, atom.res_seq, atom.res_name.upper())
        by_residue.setdefault(key, {})[atom.name.upper()] = atom
    rings: list[Ring] = []
    for (chain, res_seq, res_name), named in by_residue.items():
        for ring_names in _PROTEIN_RINGS.get(res_name, ()):
            members = [named[name] for name in ring_names if name in named]
            if len(members) != len(ring_names):
                continue
            centroid, normal, _ = _ring_plane(members)
            label = f"{res_name.capitalize()}{res_seq}"
            rings.append(Ring(atoms=members, centroid=centroid, normal=normal, label=label))
    return rings


def _ring_angle_deg(first: Ring, second: Ring) -> float:
    dot = sum(a * b for a, b in zip(first.normal, second.normal))
    dot = max(-1.0, min(1.0, abs(dot)))
    return math.degrees(math.acos(dot))


# --------------------------------------------------------------------------
# Определение контактов
# --------------------------------------------------------------------------


def _is_hbond_capable(atom: Atom) -> bool:
    return atom.element in _POLAR


def _classify_pair(
    ligand_atom: Atom,
    protein_atom: Atom,
    distance: float,
    ligand_charges: dict[int, int],
) -> str | None:
    lig_element = ligand_atom.element
    prot_element = protein_atom.element

    if protein_atom.is_metal and lig_element in _POLAR and distance <= METAL_MAX_A:
        return "metal"
    if ligand_atom.is_metal and prot_element in _POLAR and distance <= METAL_MAX_A:
        return "metal"

    lig_charge = ligand_charges.get(ligand_atom.serial, 0)
    prot_charge = protein_charge(protein_atom)
    if lig_charge and prot_charge and lig_charge * prot_charge < 0 and distance <= SALT_BRIDGE_MAX_A:
        return "salt_bridge"

    if (
        _is_hbond_capable(ligand_atom)
        and _is_hbond_capable(protein_atom)
        and distance <= HBOND_MAX_A
        and not (lig_element == "S" and prot_element == "S")
    ):
        return "hydrogen_bond"

    if lig_element in _HALOGENS and prot_element in _POLAR and distance <= HALOGEN_MAX_A:
        return "halogen_bond"

    if lig_element == "C" and prot_element == "C" and distance <= HYDROPHOBIC_MAX_A:
        return "hydrophobic"
    return None


def _keep_closest(contacts: list[Contact]) -> list[Contact]:
    """Один контакт на пару «кольцо — остаток».

    Без этого один остаток His даёт два π-катиона через ND1 и NE2, а остаток
    Arg — три через NE, NH1 и NH2, и счётчик контактов раздувается втрое
    там, где взаимодействие ровно одно.
    """
    best: dict[tuple[int, str, str], Contact] = {}
    for contact in contacts:
        key = (contact.ligand_atom_serial, contact.residue_label, contact.interaction_type)
        previous = best.get(key)
        if previous is None or contact.distance_a < previous.distance_a:
            best[key] = contact
    return list(best.values())


def _ring_contacts(
    ligand_atoms: list[Atom],
    protein_atoms: list[Atom],
    ligand_charges: dict[int, int],
) -> list[Contact]:
    contacts: list[Contact] = []
    lig_rings = ligand_rings(ligand_atoms)
    prot_rings = protein_rings(protein_atoms)

    for lig_ring in lig_rings:
        lig_anchor = lig_ring.atoms[0]
        for prot_ring in prot_rings:
            distance = math.dist(lig_ring.centroid, prot_ring.centroid)
            if distance > PI_STACK_MAX_A:
                continue
            angle = _ring_angle_deg(lig_ring, prot_ring)
            if angle > PI_STACK_PARALLEL_DEG and angle < PI_STACK_TSHAPE_MIN_DEG:
                continue
            anchor = prot_ring.atoms[0]
            contacts.append(
                Contact(
                    interaction_type="pi_stacking",
                    distance_a=round(distance, 3),
                    ligand_atom=f"ring({lig_anchor.name})",
                    ligand_atom_serial=lig_anchor.serial,
                    ligand_element="ring",
                    protein_atom=f"ring({anchor.name})",
                    protein_atom_serial=anchor.serial,
                    residue_name=anchor.res_name,
                    residue_seq=anchor.res_seq,
                    chain=anchor.chain,
                )
            )
        for protein_atom in protein_atoms:
            if protein_charge(protein_atom) <= 0:
                continue
            distance = math.dist(
                lig_ring.centroid, (protein_atom.x, protein_atom.y, protein_atom.z)
            )
            if distance <= PI_CATION_MAX_A:
                contacts.append(
                    Contact(
                        interaction_type="pi_cation",
                        distance_a=round(distance, 3),
                        ligand_atom=f"ring({lig_anchor.name})",
                        ligand_atom_serial=lig_anchor.serial,
                        ligand_element="ring",
                        protein_atom=protein_atom.name,
                        protein_atom_serial=protein_atom.serial,
                        residue_name=protein_atom.res_name,
                        residue_seq=protein_atom.res_seq,
                        chain=protein_atom.chain,
                    )
                )

    for prot_ring in prot_rings:
        anchor = prot_ring.atoms[0]
        for ligand_atom in ligand_atoms:
            if ligand_charges.get(ligand_atom.serial, 0) <= 0:
                continue
            distance = math.dist(
                prot_ring.centroid, (ligand_atom.x, ligand_atom.y, ligand_atom.z)
            )
            if distance <= PI_CATION_MAX_A:
                contacts.append(
                    Contact(
                        interaction_type="pi_cation",
                        distance_a=round(distance, 3),
                        ligand_atom=ligand_atom.name,
                        ligand_atom_serial=ligand_atom.serial,
                        ligand_element=ligand_atom.element,
                        protein_atom=f"ring({anchor.name})",
                        protein_atom_serial=anchor.serial,
                        residue_name=anchor.res_name,
                        residue_seq=anchor.res_seq,
                        chain=anchor.chain,
                    )
                )
    return _keep_closest(contacts)


def detect_contacts(
    ligand_atoms: list[Atom],
    protein_atoms: list[Atom],
    *,
    include_rings: bool = True,
) -> list[Contact]:
    """Контакты одной позы. На одну пару атомов остаётся один, самый специфичный тип."""
    ligand_heavy = heavy_atoms(ligand_atoms)
    protein_heavy = heavy_atoms(protein_atoms)
    if not ligand_heavy or not protein_heavy:
        return []
    charges = ligand_charged_atoms(ligand_atoms)

    # Отсечение по кубу вокруг лиганда: белок большой, перебор пар без него дорог.
    margin = max(PI_CATION_MAX_A, HYDROPHOBIC_MAX_A) + 1.0
    min_x = min(a.x for a in ligand_heavy) - margin
    max_x = max(a.x for a in ligand_heavy) + margin
    min_y = min(a.y for a in ligand_heavy) - margin
    max_y = max(a.y for a in ligand_heavy) + margin
    min_z = min(a.z for a in ligand_heavy) - margin
    max_z = max(a.z for a in ligand_heavy) + margin
    nearby = [
        atom
        for atom in protein_heavy
        if min_x <= atom.x <= max_x and min_y <= atom.y <= max_y and min_z <= atom.z <= max_z
    ]

    best: dict[tuple[int, int], Contact] = {}
    for ligand_atom in ligand_heavy:
        for protein_atom in nearby:
            distance = ligand_atom.distance_to(protein_atom)
            if distance > HYDROPHOBIC_MAX_A and distance > HBOND_MAX_A:
                continue
            kind = _classify_pair(ligand_atom, protein_atom, distance, charges)
            if kind is None:
                continue
            contact = Contact(
                interaction_type=kind,
                distance_a=round(distance, 3),
                ligand_atom=ligand_atom.name,
                ligand_atom_serial=ligand_atom.serial,
                ligand_element=ligand_atom.element,
                protein_atom=protein_atom.name,
                protein_atom_serial=protein_atom.serial,
                residue_name=protein_atom.res_name,
                residue_seq=protein_atom.res_seq,
                chain=protein_atom.chain,
            )
            key = contact.pair_key
            previous = best.get(key)
            if previous is None or _PRIORITY[kind] < _PRIORITY[previous.interaction_type]:
                best[key] = contact

    contacts = list(best.values())
    if include_rings:
        contacts.extend(_ring_contacts(ligand_heavy, nearby, charges))
    contacts.sort(key=lambda c: (_PRIORITY[c.interaction_type], c.distance_a))
    return contacts


def contacts_by_residue(contacts: list[Contact]) -> dict[str, set[str]]:
    """Остаток → множество типов взаимодействий с ним."""
    mapping: dict[str, set[str]] = {}
    for contact in contacts:
        mapping.setdefault(contact.residue_label, set()).add(contact.interaction_type)
    return mapping


def count_by_type(contacts: list[Contact]) -> dict[str, int]:
    counts = {name: 0 for name in INTERACTION_TYPES}
    for contact in contacts:
        counts[contact.interaction_type] += 1
    return counts
