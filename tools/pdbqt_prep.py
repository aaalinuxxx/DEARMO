"""Подготовка PDBQT для AutoDock Vina: рецептор из PDB и лиганд из SMILES."""
from __future__ import annotations

from pathlib import Path

_METALS = {"ZN", "MG", "MN", "FE", "CU", "CO", "NI"}
_SKIP_RESIDUES = {
    "HOH",
    "WAT",
    "DOD",
    "SO4",
    "GOL",
    "EDO",
    "PEG",
    "DMS",
    "ACT",
    "ACE",
    "NAG",
    "CL",
    "NA",
    "MG",
    "ZN",
    "CA",
    "K",
}

_AD_TYPE = {
    "C": "C",
    "N": "NA",
    "O": "OA",
    "S": "SA",
    "H": "HD",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
}


def ad_type(element: str) -> str:
    return _AD_TYPE.get(element.upper(), "C")


def _pdbqt_atom(serial: int, name: str, resname: str, chain: str, resseq: str, x: float, y: float, z: float, element: str) -> str:
    atom_name = f"{name:>4}"[:4]
    chain_id = (chain or "A")[:1]
    try:
        residue = f"{int(resseq):4d}"
    except ValueError:
        residue = f"{resseq[:4]:>4}"
    core = f"ATOM  {serial:5d} {atom_name} {resname[:3]:>3} {chain_id}{residue}    {x:8.3f}{y:8.3f}{z:8.3f}"
    # Колонки 71–76 — заряд, 78–79 — тип AutoDock. Пробел в заряде Vina отвергает.
    return f"{core}  1.00  0.00    +0.000 {ad_type(element):>2}"


def receptor_pdbqt(pdb_text: str) -> str:
    """Белковые ATOM без кристаллизационных вод и без HETATM лиганда."""
    lines: list[str] = []
    serial = 0
    for raw in pdb_text.splitlines():
        if raw.startswith("ATOM"):
            resname = raw[17:20].strip()
        elif raw.startswith("HETATM") and raw[17:20].strip() in _METALS:
            resname = raw[17:20].strip()
        else:
            continue
        if resname in _SKIP_RESIDUES:
            continue
        try:
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
        except ValueError:
            continue
        element = raw[76:78].strip() if len(raw) >= 78 else raw[12:16].strip()[:1]
        if not element:
            element = raw[12:16].strip()[:1]
        serial += 1
        lines.append(
            _pdbqt_atom(
                serial,
                raw[12:16].strip() or element,
                resname,
                raw[21].strip(),
                raw[22:26].strip() or "1",
                x,
                y,
                z,
                element,
            )
        )
    if not lines:
        raise ValueError("В PDB нет белковых атомов ATOM")
    return "\n".join(lines) + "\n"


def crystal_ligand_center(pdb_text: str) -> tuple[float, float, float] | None:
    """Центр связанного лиганда в кристалле, если он есть. Воды и ионы не учитываются."""
    coords: list[tuple[float, float, float]] = []
    for raw in pdb_text.splitlines():
        if not raw.startswith("HETATM"):
            continue
        resname = raw[17:20].strip()
        if resname in _SKIP_RESIDUES:
            continue
        try:
            coords.append((float(raw[30:38]), float(raw[38:46]), float(raw[46:54])))
        except ValueError:
            continue
    if len(coords) < 3:
        return None
    n = len(coords)
    return (sum(c[0] for c in coords) / n, sum(c[1] for c in coords) / n, sum(c[2] for c in coords) / n)


def protein_center(pdb_text: str) -> tuple[float, float, float]:
    coords: list[tuple[float, float, float]] = []
    for raw in pdb_text.splitlines():
        if not raw.startswith("ATOM") or raw[12:16].strip() not in {"CA", "C"}:
            continue
        try:
            coords.append((float(raw[30:38]), float(raw[38:46]), float(raw[46:54])))
        except ValueError:
            continue
    if not coords:
        raise ValueError("Не удалось найти центр белка")
    n = len(coords)
    return (sum(c[0] for c in coords) / n, sum(c[1] for c in coords) / n, sum(c[2] for c in coords) / n)


def ligand_pdbqt_from_smiles(smiles: str) -> tuple[str, str]:
    """Гибкий лиганд через Meeko, иначе жёсткий PDBQT из 3D-конформера RDKit."""
    try:
        return _meeko_ligand(smiles), "meeko"
    except Exception:
        return _rigid_ligand(smiles), "rdkit-rigid"


def _meeko_ligand(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES не разбирается: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError("RDKit не построил 3D-конформер")
    AllChem.MMFFOptimizeMolecule(mol)
    setups = MoleculePreparation().prepare(mol)
    if not setups:
        raise ValueError("Meeko не подготовил лиганд")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok or not pdbqt:
        raise ValueError(error or "пустой PDBQT лиганда")
    return pdbqt


def _rigid_ligand(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES не разбирается: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError("RDKit не построил 3D-конформер")
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    lines = ["ROOT"]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        lines.append(
            _pdbqt_atom(
                atom.GetIdx() + 1,
                symbol,
                "LIG",
                "A",
                "1",
                pos.x,
                pos.y,
                pos.z,
                symbol,
            )
        )
    lines.extend(["ENDROOT", "TORSDOF 0"])
    return "\n".join(lines) + "\n"


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path
