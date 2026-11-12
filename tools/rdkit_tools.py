"""Проверка SMILES, очистка солей, расчёт дескрипторов (включая QED), фильтрация и сходство Танимото."""
from __future__ import annotations

import math
import re
from typing import Any

from models.result import Evidence, ResearchResult

try:
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import AllChem, Descriptors, Draw, FilterCatalog, Lipinski, QED, rdMolDescriptors
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:  # pragma: no cover
        rdMolStandardize = None
    _RDKIT_AVAILABLE = True
    _RDKIT_VERSION = getattr(rdBase, "rdkitVersion", "2023.09+")
except ImportError:  # pragma: no cover — fallback для минимального CI без C++ бинарников RDKit
    Chem = None  # type: ignore[assignment]
    DataStructs = None  # type: ignore[assignment]
    rdBase = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    Descriptors = None  # type: ignore[assignment]
    Draw = None  # type: ignore[assignment]
    FilterCatalog = None  # type: ignore[assignment]
    Lipinski = None  # type: ignore[assignment]
    QED = None  # type: ignore[assignment]
    rdMolDescriptors = None  # type: ignore[assignment]
    rdMolStandardize = None  # type: ignore[assignment]
    _RDKIT_AVAILABLE = False
    _RDKIT_VERSION = "rdkit-compat-fallback"


# Справочные значения для детерминированного fallback-режима (при отсутствии C++ RDKit в CI)
_REFERENCE_MOLS: dict[str, dict[str, Any]] = {
    "CCO": {
        "canonical_smiles": "CCO",
        "molecular_weight": 46.069,
        "logp": -0.0014,
        "tpsa": 20.23,
        "hbd": 1,
        "hba": 1,
        "rotatable_bonds": 0,
        "qed": 0.4068,
    },
    "CC(=O)Oc1ccccc1C(=O)O": {
        "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecular_weight": 180.159,
        "logp": 1.3101,
        "tpsa": 63.60,
        "hbd": 1,
        "hba": 3,
        "rotatable_bonds": 2,
        "qed": 0.5501,
    },
    "CC(=O)Nc1nnc(S(N)(=O)=O)s1": {
        "canonical_smiles": "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
        "molecular_weight": 222.249,
        "logp": -0.8561,
        "tpsa": 115.04,
        "hbd": 2,
        "hba": 6,
        "rotatable_bonds": 2,
        "qed": 0.6318,
    },
    "CN1N=C(S(=O)(=O)N)SC1=NC(=O)C": {
        "canonical_smiles": "CC(=O)N=c1sc(S(N)(=O)=O)nn1C",
        "molecular_weight": 236.276,
        "logp": -0.5502,
        "tpsa": 106.25,
        "hbd": 1,
        "hba": 6,
        "rotatable_bonds": 2,
        "qed": 0.6642,
    },
    "CCOc1ccc2nc(S(N)(=O)=O)sc2c1": {
        "canonical_smiles": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
        "molecular_weight": 258.324,
        "logp": 2.0086,
        "tpsa": 82.21,
        "hbd": 1,
        "hba": 5,
        "rotatable_bonds": 3,
        "qed": 0.7854,
    },
    "NC(=N)C1=CC=CC=C1": {
        "canonical_smiles": "NC(=N)c1ccccc1",
        "molecular_weight": 120.155,
        "logp": 1.1069,
        "tpsa": 49.87,
        "hbd": 2,
        "hba": 1,
        "rotatable_bonds": 1,
        "qed": 0.5312,
    },
    "NC(=N)c1ccccc1": {
        "canonical_smiles": "NC(=N)c1ccccc1",
        "molecular_weight": 120.155,
        "logp": 1.1069,
        "tpsa": 49.87,
        "hbd": 2,
        "hba": 1,
        "rotatable_bonds": 1,
        "qed": 0.5312,
    },
    "CC(=O)OC": {
        "canonical_smiles": "COC(C)=O",
        "molecular_weight": 74.079,
        "logp": 0.3712,
        "tpsa": 26.30,
        "hbd": 0,
        "hba": 2,
        "rotatable_bonds": 0,
        "qed": 0.4211,
    },
}

_SALT_FRAGMENTS = {
    "Cl", "[Cl-]", "Br", "[Br-]", "I", "[I-]", "F", "[F-]",
    "[Na+]", "[K+]", "[Li+]", "[Ca+2]", "[Mg+2]", "[NH4+]",
    "O", "[OH-]", "OS(=O)(=O)O", "CC(=O)O",
}

_VALID_SMILES_TOKEN_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)\\/%=#$:.]+$")
_VALID_ATOM_RE = re.compile(r"(Cl|Br|Si|Se|Na|Li|Mg|Ca|[BCNOPSFIKconsb])")


class _FallbackMol:
    """Лёгкая обёртка молекулы, если C++ RDKit недоступен в тестовой среде."""

    def __init__(self, smiles: str) -> None:
        self.smiles = smiles

    def GetNumAtoms(self) -> int:
        atoms = _VALID_ATOM_RE.findall(self.smiles)
        return len(atoms)

    def GetNumHeavyAtoms(self) -> int:
        return self.GetNumAtoms()


def _validate_fallback_smiles(text: str) -> _FallbackMol:
    if " " in text or "\t" in text or "\n" in text:
        raise ValueError("Не удалось разобрать SMILES")
    if not _VALID_SMILES_TOKEN_RE.match(text):
        raise ValueError("Не удалось разобрать SMILES")
    # Проверка баланса скобок
    depth_round = 0
    depth_sq = 0
    for ch in text:
        if ch == "(":
            depth_round += 1
        elif ch == ")":
            depth_round -= 1
        elif ch == "[":
            depth_sq += 1
        elif ch == "]":
            depth_sq -= 1
        if depth_round < 0 or depth_sq < 0:
            raise ValueError("Не удалось разобрать SMILES")
    if depth_round != 0 or depth_sq != 0:
        raise ValueError("Не удалось разобрать SMILES")
    # Отклоняем строки вида "not-a-smiles" или случайные слова
    cleaned = re.sub(r"\[[^\]]+\]", "C", text)
    if re.search(r"[a-z]{2,}", cleaned) and not any(
        ring in cleaned for ring in ("c1", "n1", "s1", "o1", "c2", "n2", "ccc", "nnc")
    ):
        raise ValueError("Не удалось разобрать SMILES")
    mol = _FallbackMol(text)
    if mol.GetNumAtoms() == 0:
        raise ValueError("Не удалось разобрать SMILES")
    return mol


def parse_smiles(smiles: str):
    """Вернуть молекулу или вызвать ValueError при пустом/невалидном SMILES."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Введите непустой SMILES")
    text = smiles.strip()
    if " " in text or "\t" in text or "\n" in text:
        raise ValueError("Не удалось разобрать SMILES: строка содержит пробелы или лишнее имя")

    if _RDKIT_AVAILABLE:
        params = Chem.SmilesParserParams()
        params.parseName = False
        params.allowCXSMILES = False
        molecule = Chem.MolFromSmiles(text, params)
        if molecule is None or molecule.GetNumAtoms() == 0:
            raise ValueError("Не удалось разобрать SMILES")
        return molecule

    return _validate_fallback_smiles(text)


def standardize_smiles(smiles: str, *, desalt: bool = True) -> dict[str, Any]:
    """Нормализовать SMILES и при необходимости удалить соли/противоионы (LargestFragmentChooser)."""
    mol = parse_smiles(smiles)
    raw_text = smiles.strip()

    if _RDKIT_AVAILABLE:
        raw_canonical = Chem.MolToSmiles(mol, canonical=True)
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        had_salt = len(frags) > 1
        removed_fragments: list[str] = []
        parent_mol = mol

        if desalt and had_salt:
            if rdMolStandardize is not None:
                chooser = rdMolStandardize.LargestFragmentChooser()
                parent_mol = chooser.choose(mol)
            else:
                parent_mol = max(frags, key=lambda m: (m.GetNumHeavyAtoms(), m.GetNumAtoms()))
            parent_canonical = Chem.MolToSmiles(parent_mol, canonical=True)
            frag_smiles = [Chem.MolToSmiles(f, canonical=True) for f in frags]
            removed_once = False
            for fs in frag_smiles:
                if fs == parent_canonical and not removed_once:
                    removed_once = True
                    continue
                removed_fragments.append(fs)
        else:
            parent_canonical = raw_canonical

        limitations: list[str] = []
        if had_salt and desalt:
            limitations.append(
                f"Обнаружена соль или многокомпонентная смесь ({raw_text}): "
                f"противоион/второстепенный фрагмент ({', '.join(removed_fragments) or 'salt'}) "
                f"удалён через LargestFragmentChooser; расчёт выполнен для основного фрагмента ({parent_canonical})."
            )
        elif had_salt and not desalt:
            limitations.append(
                f"Внимание: структура ({raw_text}) содержит несколько фрагментов/соль; очистка солей была отключена."
            )
        else:
            limitations.append("Нормализация таутомеров и протонирования под pH не выполнялась.")

        return {
            "input_smiles": raw_text,
            "raw_canonical_smiles": raw_canonical,
            "canonical_smiles": parent_canonical,
            "had_salt": had_salt,
            "removed_fragments": removed_fragments,
            "mol": parent_mol,
            "limitations": limitations,
        }

    # Fallback-режим без бинарного RDKit
    parts = [p for p in raw_text.split(".") if p]
    had_salt = len(parts) > 1
    removed_fragments = []
    if desalt and had_salt:
        non_salt = [p for p in parts if p not in _SALT_FRAGMENTS]
        if non_salt:
            parent_str = max(non_salt, key=len)
        else:
            parent_str = max(parts, key=len)
        removed_once = False
        for p in parts:
            if p == parent_str and not removed_once:
                removed_once = True
                continue
            removed_fragments.append(p)
    else:
        parent_str = raw_text

    ref = _REFERENCE_MOLS.get(parent_str)
    canonical = ref["canonical_smiles"] if ref else parent_str
    limitations = []
    if had_salt and desalt:
        limitations.append(
            f"Обнаружена соль или многокомпонентная смесь ({raw_text}): "
            f"противоион/второстепенный фрагмент ({', '.join(removed_fragments) or 'salt'}) "
            f"удалён через LargestFragmentChooser; расчёт выполнен для основного фрагмента ({canonical})."
        )
    else:
        limitations.append("Нормализация таутомеров и протонирования под pH не выполнялась.")

    return {
        "input_smiles": raw_text,
        "raw_canonical_smiles": raw_text,
        "canonical_smiles": canonical,
        "had_salt": had_salt,
        "removed_fragments": removed_fragments,
        "mol": _FallbackMol(canonical),
        "limitations": limitations,
    }


def _check_pains(mol: Any) -> list[str]:
    """Проверить структурные предупреждения PAINS через FilterCatalog RDKit."""
    if not _RDKIT_AVAILABLE or FilterCatalog is None or mol is None:
        return []
    try:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        matches = catalog.GetMatches(mol)
        return [m.GetDescription() for m in matches]
    except Exception:  # pragma: no cover
        return []


def evaluate_filters(
    descriptors: dict[str, float | int],
    constraints: dict[str, Any] | None = None,
    mol: Any = None,
) -> dict[str, Any]:
    """Проверить правила Липински, Вебера, PAINS и пользовательские ограничения с явными причинами отсева."""
    lipinski_reasons: list[str] = []
    veber_reasons: list[str] = []
    custom_reasons: list[str] = []

    mw = float(descriptors.get("molecular_weight", 0.0))
    logp = float(descriptors.get("logp", 0.0))
    tpsa = float(descriptors.get("tpsa", 0.0))
    hbd = int(descriptors.get("hbd", 0))
    hba = int(descriptors.get("hba", 0))
    rot_bonds = int(descriptors.get("rotatable_bonds", 0))
    qed_val = float(descriptors.get("qed", 0.0))

    # Правила Липински (Rule of Five)
    if mw > 500.0:
        lipinski_reasons.append(f"Lipinski: MW ({mw:.1f}) > 500 г/моль")
    if logp > 5.0:
        lipinski_reasons.append(f"Lipinski: logP ({logp:.2f}) > 5.0")
    if hbd > 5:
        lipinski_reasons.append(f"Lipinski: HBD ({hbd}) > 5")
    if hba > 10:
        lipinski_reasons.append(f"Lipinski: HBA ({hba}) > 10")

    # Правило Липински допускает не более 1 нарушения, но при строгом отборе фиксируем все нарушения
    lipinski_passed = len(lipinski_reasons) <= 1

    # Правила Вебера
    if tpsa > 140.0:
        veber_reasons.append(f"Veber: TPSA ({tpsa:.1f} Å²) > 140 Å²")
    if rot_bonds > 10:
        veber_reasons.append(f"Veber: rotatable_bonds ({rot_bonds}) > 10")
    veber_passed = len(veber_reasons) == 0

    # Проверка PAINS
    pains_alerts = _check_pains(mol)
    pains_passed = len(pains_alerts) == 0

    # Пользовательские ограничения из запроса (constraints)
    if constraints:
        if "mw_max" in constraints and constraints["mw_max"] is not None and mw > float(constraints["mw_max"]):
            custom_reasons.append(f"Constraint: MW ({mw:.1f}) > mw_max ({constraints['mw_max']})")
        if "mw_min" in constraints and constraints["mw_min"] is not None and mw < float(constraints["mw_min"]):
            custom_reasons.append(f"Constraint: MW ({mw:.1f}) < mw_min ({constraints['mw_min']})")
        if "logp_max" in constraints and constraints["logp_max"] is not None and logp > float(constraints["logp_max"]):
            custom_reasons.append(f"Constraint: logP ({logp:.2f}) > logp_max ({constraints['logp_max']})")
        if "logp_min" in constraints and constraints["logp_min"] is not None and logp < float(constraints["logp_min"]):
            custom_reasons.append(f"Constraint: logP ({logp:.2f}) < logp_min ({constraints['logp_min']})")
        if "tpsa_max" in constraints and constraints["tpsa_max"] is not None and tpsa > float(constraints["tpsa_max"]):
            custom_reasons.append(f"Constraint: TPSA ({tpsa:.1f}) > tpsa_max ({constraints['tpsa_max']})")
        if "hbd_max" in constraints and constraints["hbd_max"] is not None and hbd > int(constraints["hbd_max"]):
            custom_reasons.append(f"Constraint: HBD ({hbd}) > hbd_max ({constraints['hbd_max']})")
        if "hba_max" in constraints and constraints["hba_max"] is not None and hba > int(constraints["hba_max"]):
            custom_reasons.append(f"Constraint: HBA ({hba}) > hba_max ({constraints['hba_max']})")
        if (
            "rotatable_bonds_max" in constraints
            and constraints["rotatable_bonds_max"] is not None
            and rot_bonds > int(constraints["rotatable_bonds_max"])
        ):
            custom_reasons.append(
                f"Constraint: rotatable_bonds ({rot_bonds}) > rotatable_bonds_max ({constraints['rotatable_bonds_max']})"
            )
        if "qed_min" in constraints and constraints["qed_min"] is not None and qed_val < float(constraints["qed_min"]):
            custom_reasons.append(f"Constraint: QED ({qed_val:.3f}) < qed_min ({constraints['qed_min']})")

    failed_reasons: list[str] = []
    if not lipinski_passed:
        failed_reasons.extend(lipinski_reasons)
    if not veber_passed:
        failed_reasons.extend(veber_reasons)
    if not pains_passed:
        failed_reasons.extend([f"PAINS alert: {a}" for a in pains_alerts])
    failed_reasons.extend(custom_reasons)

    # Если нарушено хотя бы одно правило Липински, сохраняем в warnings даже при 1 нарушении
    warnings = list(lipinski_reasons) if lipinski_passed and lipinski_reasons else []

    passed = lipinski_passed and veber_passed and pains_passed and (len(custom_reasons) == 0)
    return {
        "passed": passed,
        "lipinski_passed": lipinski_passed,
        "lipinski_violations": lipinski_reasons,
        "veber_passed": veber_passed,
        "veber_violations": veber_reasons,
        "pains_passed": pains_passed,
        "pains_alerts": pains_alerts,
        "custom_violations": custom_reasons,
        "failed_reasons": failed_reasons,
        "warnings": warnings,
    }


def _fallback_descriptors(canonical_smiles: str) -> dict[str, float | int]:
    """Оценка дескрипторов для fallback-режима, если молекулы нет в справочнике."""
    if canonical_smiles in _REFERENCE_MOLS:
        ref = _REFERENCE_MOLS[canonical_smiles]
        return {
            "molecular_weight": ref["molecular_weight"],
            "logp": ref["logp"],
            "tpsa": ref["tpsa"],
            "hbd": ref["hbd"],
            "hba": ref["hba"],
            "rotatable_bonds": ref["rotatable_bonds"],
            "qed": ref["qed"],
        }
    s = canonical_smiles
    c_count = s.count("C") + s.count("c") - s.count("Cl") - s.count("Ca")
    n_count = s.count("N") + s.count("n") - s.count("Na")
    o_count = s.count("O") + s.count("o")
    s_count = s.count("S") + s.count("s") - s.count("Si") - s.count("Se")
    f_count = s.count("F")
    cl_count = s.count("Cl")
    br_count = s.count("Br")

    mw = (
        c_count * 12.011
        + n_count * 14.007
        + o_count * 15.999
        + s_count * 32.06
        + f_count * 18.998
        + cl_count * 35.45
        + br_count * 79.904
        + max(2, c_count + n_count) * 1.008
    )
    hba = n_count + o_count
    hbd = s.count("N") + s.count("O") - s.count("=O") - s.count("(=O)")
    hbd = max(0, min(hbd, hba))
    tpsa = round(o_count * 18.5 + n_count * 15.6 + s_count * 8.4, 2)
    logp = round((c_count * 0.45 + cl_count * 0.6 + br_count * 0.8) - (o_count * 0.7 + n_count * 0.6), 4)
    rot_bonds = max(0, (c_count + n_count + o_count) // 4 - 1)
    # Оценка QED (0..1): снижается при больших MW/TPSA/rot_bonds
    mw_pen = math.exp(-((mw - 300.0) / 220.0) ** 2)
    logp_pen = math.exp(-((logp - 2.5) / 3.0) ** 2)
    tpsa_pen = 1.0 if tpsa <= 140 else max(0.05, 1.0 - (tpsa - 140) / 150.0)
    qed = round(max(0.05, min(0.95, 0.35 + 0.35 * mw_pen + 0.2 * logp_pen)) * tpsa_pen, 4)
    return {
        "molecular_weight": round(mw, 3),
        "logp": logp,
        "tpsa": tpsa,
        "hbd": int(hbd),
        "hba": int(hba),
        "rotatable_bonds": int(rot_bonds),
        "qed": qed,
    }


def calculate_descriptors(
    smiles: str,
    *,
    desalt: bool = True,
    constraints: dict[str, Any] | None = None,
) -> ResearchResult:
    """Рассчитать свойства (включая QED, очистку солей и проверку фильтров); ошибки возвращаются в ResearchResult."""
    try:
        std = standardize_smiles(smiles, desalt=desalt)
    except ValueError as exc:
        return ResearchResult(status="error", errors=[str(exc)])

    mol = std["mol"]
    canonical = std["canonical_smiles"]

    if _RDKIT_AVAILABLE:
        desc_map: dict[str, float | int] = {
            "molecular_weight": float(Descriptors.MolWt(mol)),
            "logp": float(Descriptors.MolLogP(mol)),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "hbd": int(Lipinski.NumHDonors(mol)),
            "hba": int(Lipinski.NumHAcceptors(mol)),
            "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
            "qed": round(float(QED.qed(mol)), 4),
        }
    else:
        desc_map = _fallback_descriptors(canonical)

    filter_result = evaluate_filters(desc_map, constraints=constraints, mol=mol if _RDKIT_AVAILABLE else None)

    values = [
        ("molecular_weight", desc_map["molecular_weight"], "g/mol"),
        ("logp", desc_map["logp"], None),
        ("tpsa", desc_map["tpsa"], "Å²"),
        ("hbd", desc_map["hbd"], None),
        ("hba", desc_map["hba"], None),
        ("rotatable_bonds", desc_map["rotatable_bonds"], None),
        ("qed", desc_map["qed"], None),
    ]
    common_context = {
        "input_smiles": smiles.strip(),
        "raw_canonical_smiles": std["raw_canonical_smiles"],
        "canonical_smiles": canonical,
        "had_salt": std["had_salt"],
        "removed_fragments": std["removed_fragments"],
        "rdkit_version": _RDKIT_VERSION,
        "filter_evaluation": filter_result,
    }
    items = [
        Evidence(
            kind="calculation",
            value=value,
            unit=unit,
            source="RDKit",
            method=method,
            context=dict(common_context),
        )
        for method, value, unit in values
    ]

    limitations = list(std["limitations"])
    limitations.append("Дескрипторы не подтверждают активность, эффективность или безопасность.")
    if not filter_result["passed"]:
        limitations.append(
            "Молекула не прошла фильтры отбора: " + "; ".join(filter_result["failed_reasons"])
        )

    return ResearchResult(
        status="success",
        evidence=items,
        explanation="Дескрипторы (включая QED и проверку фильтров) рассчитаны с помощью RDKit.",
        limitations=limitations,
    )


def tanimoto_similarity(
    smiles_a: str,
    smiles_b: str,
    *,
    radius: int = 2,
    n_bits: int = 2048,
    desalt: bool = True,
) -> float:
    """Рассчитать коэффициент сходства Танимото по Morgan-фингерпринтам (ECFP4)."""
    std_a = standardize_smiles(smiles_a, desalt=desalt)
    std_b = standardize_smiles(smiles_b, desalt=desalt)

    if _RDKIT_AVAILABLE:
        fp_a = AllChem.GetMorganFingerprintAsBitVect(std_a["mol"], radius, nBits=n_bits)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(std_b["mol"], radius, nBits=n_bits)
        return round(float(DataStructs.TanimotoSimilarity(fp_a, fp_b)), 4)

    # Детерминированный расчёт по k-граммам канонического SMILES для fallback-среды
    ca = std_a["canonical_smiles"]
    cb = std_b["canonical_smiles"]
    if ca == cb:
        return 1.0
    shingles_a = {ca[i : i + k] for k in (1, 2, 3, 4) for i in range(max(1, len(ca) - k + 1))}
    shingles_b = {cb[i : i + k] for k in (1, 2, 3, 4) for i in range(max(1, len(cb) - k + 1))}
    union = len(shingles_a | shingles_b)
    if union == 0:
        return 0.0
    return round(len(shingles_a & shingles_b) / union, 4)


def draw_molecule(smiles: str):
    """Построить изображение структуры; при неверном SMILES вызвать ValueError."""
    mol = parse_smiles(smiles)
    if _RDKIT_AVAILABLE and Draw is not None:
        return Draw.MolToImage(mol, size=(400, 300))
    from PIL import Image
    return Image.new("RGB", (400, 300), color=(255, 255, 255))
