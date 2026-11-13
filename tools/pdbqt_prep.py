"""Подготовка PDBQT для AutoDock Vina: рецептор из PDB и лиганд из SMILES."""
from __future__ import annotations

import shutil
import subprocess
import sys
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
    # Ионы металлов должны сохранять свой тип. Раньше их здесь не было, и
    # запасной вариант записывал цинк карбоангидразы как углерод: Vina считала
    # бы по нему углеродный потенциал, а сценарии 3 и 4 не отличили бы металл
    # от атома белка.
    "ZN": "Zn",
    "MG": "Mg",
    "MN": "Mn",
    "FE": "Fe",
    "CU": "Cu",
    "CO": "Co",
    "NI": "Ni",
    "CA": "Ca",
    "K": "K",
    "NA_ION": "Na",
}


def ad_type(element: str) -> str:
    """Тип AutoDock по символу элемента.

    Неизвестный элемент раньше молча становился углеродом. Теперь он остаётся
    самим собой: подменять элемент в расчётном файле нельзя, лучше пусть Vina
    отвергнет незнакомый тип явно, чем посчитает не ту химию молча.
    """
    key = element.strip().upper()
    if key in _AD_TYPE:
        return _AD_TYPE[key]
    return key.capitalize() if key else "C"


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


def clean_receptor_pdb(pdb_text: str, chain: str | None = None) -> str:
    """Рецептор как PDB: белок и ионы металлов, без вод и без кристаллического лиганда.

    Нужен отдельно от PDBQT по двум причинам. P2Rank ищет карманы по PDB, и
    оставленный в файле лиганд подсказал бы ему ответ. А в PDBQT теряются имена
    и нумерация остатков, поэтому сценарии 3 и 4 не смогли бы написать «Asp189»
    вместо номера атома.

    Ионы металлов сохраняются: Zn²⁺ карбоангидразы — часть активного центра, и
    без него докинг ацетазоламида считал бы не ту задачу.
    """
    lines: list[str] = []
    for raw in pdb_text.splitlines():
        if raw.startswith("ATOM"):
            resname = raw[17:20].strip()
        elif raw.startswith("HETATM") and raw[17:20].strip().upper() in _METALS:
            resname = raw[17:20].strip()
        else:
            continue
        # Воды и криопротекторы отбрасываются всегда, металлы оставляем, даже
        # если их код совпал со списком пропуска.
        if resname.upper() in _SKIP_RESIDUES and resname.upper() not in _METALS:
            continue
        if raw.startswith("ATOM") and resname.upper() in _SKIP_RESIDUES:
            continue
        if chain and raw[21:22].strip().upper() != chain.strip().upper():
            continue
        lines.append(raw)
    if not lines:
        raise ValueError(
            "В PDB не осталось атомов рецептора" + (f" для цепи {chain}" if chain else "")
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def receptor_pdbqt(pdb_text: str, chain: str | None = None) -> str:
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
        if resname in _SKIP_RESIDUES and resname.upper() not in _METALS:
            continue
        if raw.startswith("ATOM") and resname in _SKIP_RESIDUES:
            continue
        if chain and raw[21:22].strip().upper() != chain.strip().upper():
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


# Заряды ионов металлов для PDBQT. Vina их не использует в скоринге, но формат
# требует число, а неверный заряд сбивает читающих файл людей.
_METAL_CHARGE = {"ZN": 2.0, "MG": 2.0, "MN": 2.0, "FE": 2.0, "CU": 2.0, "CO": 2.0, "NI": 2.0}


def _tool(name: str) -> str | None:
    """Исполняемый файл из того же окружения, где запущен Python."""
    sibling = Path(sys.executable).parent / name
    if sibling.is_file():
        return str(sibling)
    return shutil.which(name)


def _metal_lines(clean_pdb_text: str, start_serial: int) -> list[str]:
    """Ионы металлов из очищенного PDB, оформленные как HETATM для PDBQT."""
    lines: list[str] = []
    serial = start_serial
    for raw in clean_pdb_text.splitlines():
        if not raw.startswith("HETATM"):
            continue
        resname = raw[17:20].strip().upper()
        if resname not in _METALS:
            continue
        try:
            x, y, z = float(raw[30:38]), float(raw[38:46]), float(raw[46:54])
        except ValueError:
            continue
        serial += 1
        chain = (raw[21:22].strip() or "A")[:1]
        resseq = raw[22:26].strip() or "1"
        try:
            residue = f"{int(resseq):4d}"
        except ValueError:
            residue = f"{resseq[:4]:>4}"
        atom_name = f"{raw[12:16].strip() or resname:>4}"[:4]
        lines.append(
            f"HETATM{serial:5d} {atom_name} {resname[:3]:>3} {chain}{residue}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    {_METAL_CHARGE.get(resname, 2.0):+6.3f} {ad_type(resname):>2}"
        )
    return lines


def _free_format_pqr(pqr_text: str) -> tuple[str, int]:
    """Переписать PQR в свободном формате, который читает Meeko.

    pdb2pqr выдаёт строки по фиксированным колонкам, и тогда два случая ломают
    парсер Meeko, делящий строку по пробелам: у остатков с номером от 1000 цепь
    и номер слипаются («A1007», киназа 4WKQ), а вставочный код приклеивается к
    номеру («184A», трипсин 3PTB). Meeko умеет читать цепь, номер и вставочный
    код как отдельные поля, поэтому строки просто разбираются по колонкам и
    собираются заново через пробел. Нумерация остатков при этом не меняется.

    Возвращает (текст, число остатков со вставочным кодом).
    """
    out: list[str] = []
    seen_icode: set[tuple[str, str, str]] = set()
    for line in pqr_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            out.append(line)
            continue
        record = line[:6].strip()
        serial = line[6:11].strip()
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resnum = line[22:26].strip()
        icode = line[26:27].strip()
        tail = line[27:].split()
        if len(tail) < 5 or not (record and serial and name and resname and resnum):
            out.append(line)
            continue
        if icode:
            seen_icode.add((chain, resnum, icode))
        fields = [record, serial, name, resname, chain, resnum, icode, *tail[:5]]
        out.append(" ".join(field for field in fields if field))
    return "\n".join(out) + "\n", len(seen_icode)


def protonated_receptor_pdbqt(clean_pdb_path: Path | str, workdir: Path | str) -> tuple[str, list[str]]:
    """Рецептор с водородами и корректными типами атомов: pdb2pqr + Meeko.

    Зачем это нужно. Прежняя подготовка (`receptor_pdbqt`) записывала координаты
    как есть: все азоты становились акцепторами NA, водородов не было вовсе,
    ароматические углероды шли как алифатические. Для Vina белок без полярных
    водородов не имеет доноров водородных связей, и скоринг теряет половину
    контактов. На 1M17 поза №1 эрлотиниба при этом лежала в 7.9 Å от кристалла,
    после протонирования — в 1.2 Å.

    Ионы металлов pdb2pqr отбрасывает, поэтому они возвращаются из исходного PDB:
    без Zn²⁺ докинг в карбоангидразу считал бы не ту задачу.

    Возвращает (текст PDBQT, заметки). При любой неудаче бросает RuntimeError,
    и вызывающий откатывается на грубую схему с явным предупреждением.
    """
    pdb2pqr = _tool("pdb2pqr")
    meeko = _tool("mk_prepare_receptor.py")
    if not pdb2pqr or not meeko:
        raise RuntimeError("pdb2pqr или mk_prepare_receptor.py не установлены (pip install pdb2pqr meeko)")

    clean_pdb_path = Path(clean_pdb_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    pqr = workdir / "receptor.pqr"
    base = workdir / "receptor_meeko"

    proc = subprocess.run(
        [pdb2pqr, "--ff=AMBER", "--with-ph=7.4", "--drop-water", "--keep-chain", str(clean_pdb_path), str(pqr)],
        capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0 or not pqr.is_file():
        raise RuntimeError(f"pdb2pqr не справился: {(proc.stderr or proc.stdout)[-400:]}")

    free_format, with_icode = _free_format_pqr(pqr.read_text(encoding="utf-8"))
    pqr.write_text(free_format, encoding="utf-8")

    proc = subprocess.run(
        [meeko, "--read_pqr", str(pqr), "-o", str(base), "-p", "--charge_model", "read"],
        capture_output=True, text=True, timeout=600,
    )
    pdbqt = Path(f"{base}.pdbqt")
    if proc.returncode != 0 or not pdbqt.is_file():
        raise RuntimeError(f"Meeko не подготовил рецептор: {(proc.stderr or proc.stdout)[-400:]}")

    kept = [ln for ln in pdbqt.read_text(encoding="utf-8").splitlines() if ln.startswith(("ATOM", "HETATM"))]
    if not kept:
        raise RuntimeError("Meeko вернул пустой рецептор")

    metals = _metal_lines(clean_pdb_path.read_text(encoding="utf-8", errors="replace"), len(kept))
    notes = [
        "Рецептор протонирован pdb2pqr (AMBER, pH 7.4, PROPKA) и типизирован Meeko; "
        "состояния протонирования не оптимизировались под конкретный лиганд.",
    ]
    if with_icode:
        notes.append(f"В структуре есть остатки со вставочными кодами ({with_icode}); нумерация сохранена.")
    if metals:
        notes.append(
            f"Ионы металлов ({len(metals)}) возвращены из исходной структуры: pdb2pqr их отбрасывает. "
            "Протонирование остатков, координирующих металл, не уточнялось."
        )
    return "\n".join(kept + metals) + "\n", notes


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
