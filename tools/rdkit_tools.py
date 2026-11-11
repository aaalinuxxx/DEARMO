"""Проверка SMILES и вычисление дескрипторов с RDKit.

Поддерживаются SMILES без имени или дополнительных полей. Нормализация солей,
зарядов и таутомеров не выполняется: расчёт идёт по введённой структуре.
"""
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Draw, Lipinski, rdMolDescriptors
from models.result import Evidence, ResearchResult

# TODO (БИОЛОГИ): подтвердить перечень свойств и правила работы с солями,
# заряженными структурами, смесями и таутомерами для выбранной задачи.
# TODO (IT 2): согласовать политику стандартных ошибок и журналирования RDKit.

def parse_smiles(smiles: str) -> Chem.Mol:
    """Вернуть молекулу или вызвать ValueError при пустом/невалидном SMILES."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Введите непустой SMILES")
    text = smiles.strip()
    # Исключаем дополнительные названия: они иначе могут интерпретироваться RDKit.
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    molecule = Chem.MolFromSmiles(text, params)
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError("Не удалось разобрать SMILES")
    return molecule


def calculate_descriptors(smiles: str) -> ResearchResult:
    """Рассчитать свойства; ошибки ввода возвращаются в ResearchResult."""
    try:
        mol = parse_smiles(smiles)
    except ValueError as exc:
        return ResearchResult(status="error", errors=[str(exc)])

    canonical = Chem.MolToSmiles(mol, canonical=True)
    values = [
        ("molecular_weight", Descriptors.MolWt(mol), "g/mol"),
        ("logp", Descriptors.MolLogP(mol), None),
        ("tpsa", rdMolDescriptors.CalcTPSA(mol), "Å²"),
        ("hbd", Lipinski.NumHDonors(mol), None),
        ("hba", Lipinski.NumHAcceptors(mol), None),
        ("rotatable_bonds", Lipinski.NumRotatableBonds(mol), None),
    ]
    items = [
        Evidence(kind="calculation", value=value, unit=unit, source="RDKit",
                 method=method, context={"input_smiles": smiles.strip(),
                 "canonical_smiles": canonical, "rdkit_version": rdBase.rdkitVersion})
        for method, value, unit in values
    ]
    return ResearchResult(status="success", evidence=items,
        explanation="Дескрипторы рассчитаны для введённой структуры с помощью RDKit.",
        limitations=["Нормализация солей, зарядов и таутомеров не выполнялась.",
                     "Дескрипторы не подтверждают активность, эффективность или безопасность."])


def draw_molecule(smiles: str):
    """Построить изображение структуры; при неверном SMILES вызвать ValueError."""
    return Draw.MolToImage(parse_smiles(smiles), size=(400, 300))
