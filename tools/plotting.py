"""Графики сравнения докинга по карманам."""
from __future__ import annotations

from pathlib import Path

from models.docking import PoseRow, best_vina_per_ligand_pocket

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None


def pocket_ligand_matrix(rows: list[PoseRow]) -> tuple[list[str], list[str], list[list[float]]]:
    best = best_vina_per_ligand_pocket(rows)
    ligands = sorted({k[0] for k in best})
    pockets = sorted({k[1] for k in best}, key=lambda p: (int(p.replace("pocket", "") or 0) if p.replace("pocket", "").isdigit() else 99, p))
    matrix = []
    for pocket in pockets:
        row_vals = []
        for lig in ligands:
            cell = best.get((lig, pocket))
            row_vals.append(cell.vina_score if cell else float("nan"))
        matrix.append(row_vals)
    return ligands, pockets, matrix


def write_heatmap(rows: list[PoseRow], output_path: Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None:
        output_path.write_text("matplotlib not installed", encoding="utf-8")
        return output_path
    ligands, pockets, matrix = pocket_ligand_matrix(rows)
    fig, ax = plt.subplots(figsize=(max(4, len(ligands)), max(3, len(pockets))))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(ligands)), labels=ligands, rotation=45, ha="right")
    ax.set_yticks(range(len(pockets)), labels=pockets)
    ax.set_title("Лучший Vina score по лиганд × карман")
    fig.colorbar(im, ax=ax, label="Vina score (kcal/mol)")
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if val == val:  # not nan
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=120)
    plt.close(fig)
    return output_path


# --------------------------------------------------------------------------
# Графики сценария 3: сравнение поз
#
# Палитра и правила оформления — валидированный набор по слотам, светлая
# поверхность: PNG уходит в отчёт и на слайд, тёмного режима у него нет.
# Категориальные цвета назначаются по фиксированному порядку слотов и не
# перебираются по кругу: цвет закреплён за сущностью, а не за её номером.
# --------------------------------------------------------------------------

_SURFACE = "#fcfcfb"
_INK_PRIMARY = "#0b0b0b"
_INK_SECONDARY = "#52514e"
_INK_MUTED = "#898781"
_GRID = "#e1e0d9"
_AXIS = "#c3c2b7"

# Категориальные слоты в фиксированном порядке.
_SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
# Расходящаяся пара с нейтральной серединой: «исчез» ↔ «нет изменений» ↔ «появился».
_DIVERGING_LOST = "#2a78d6"
_DIVERGING_NEUTRAL = "#f0efec"
_DIVERGING_GAINED = "#e34948"
# Один синий для величины: оценки Vina — одна серия, не категории.
_MAGNITUDE = "#2a78d6"
_MAGNITUDE_REFERENCE = "#184f95"

# Короткие коды типов взаимодействий. Цвет не несёт тип в одиночку:
# в ячейке стоит код, а подпись расшифровывает его под графиком.
INTERACTION_CODES = {
    "metal": "M",
    "salt_bridge": "S",
    "hydrogen_bond": "H",
    "halogen_bond": "X",
    "pi_stacking": "π",
    "pi_cation": "π+",
    "hydrophobic": "Φ",
}

_CODE_LEGEND = "M — металл · S — солевой мостик · H — водородная связь · X — галогенная · π — стэкинг · π+ — π-катион · Φ — гидрофобный"


def _placeholder(output_path: Path, message: str) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(message, encoding="utf-8")
    return output_path


def _style_axes(ax) -> None:
    ax.set_facecolor(_SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(_AXIS)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=_INK_MUTED, labelsize=9, length=0)
    ax.title.set_color(_INK_PRIMARY)


def write_vina_by_pose(deltas, output_path: Path) -> Path:
    """Оценки Vina по позам и разница с позой отсчёта, ккал/моль.

    Одна серия, поэтому один синий и без легенды: заголовок её называет.
    Поза отсчёта выделена более тёмным шагом того же синего и подписью.
    """
    output_path = Path(output_path)
    if plt is None:
        return _placeholder(output_path, "matplotlib not installed")
    if not deltas:
        return _placeholder(output_path, "нет поз для сравнения")

    labels = [item.key.label for item in deltas]
    scores = [item.vina_score for item in deltas]
    colors = [_MAGNITUDE_REFERENCE if item.is_reference else _MAGNITUDE for item in deltas]

    height = max(2.4, 0.42 * len(deltas) + 1.5)
    fig, ax = plt.subplots(figsize=(8.0, height), facecolor=_SURFACE)
    positions = range(len(deltas))
    ax.barh(list(positions), scores, color=colors, height=0.62)
    ax.set_yticks(list(positions), labels=labels)
    ax.invert_yaxis()
    ax.set_xlabel("Оценка Vina, ккал/моль", color=_INK_SECONDARY, fontsize=10)
    ax.set_title("Оценки Vina по позам и разница с позой отсчёта", fontsize=12, pad=12)
    ax.xaxis.grid(True, color=_GRID, linewidth=1.0)
    ax.set_axisbelow(True)
    _style_axes(ax)

    # Оценки Vina отрицательные, поэтому столбец растёт влево от нуля, а его
    # «внутренняя» сторона — справа от конца столбца. Подпись белым слева от
    # конца легла бы на фон и пропала.
    span = max(abs(min(scores)), 0.001)
    pad = 0.015 * span
    for index, item in enumerate(deltas):
        if item.is_reference:
            note = f"{item.vina_score:.2f} ккал/моль · поза №1"
        else:
            note = f"{item.vina_score:.2f} · {item.delta_to_reference:+.2f}"
        fits_inside = abs(item.vina_score) > 0.45 * span
        ax.text(
            item.vina_score + pad if fits_inside else item.vina_score - pad,
            index,
            note,
            va="center",
            ha="left" if fits_inside else "right",
            fontsize=9,
            color="#ffffff" if fits_inside else _INK_SECONDARY,
        )
    ax.set_xlim(min(scores) * 1.18, 0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    return output_path


def write_contact_map(matrix, output_path: Path, *, max_residues: int = 24) -> Path:
    """Карта «поза × остаток»: где контакты появляются и где исчезают.

    Заливка расходящаяся относительно позы отсчёта: красный — контакт
    появился, синий — исчез, нейтральный серый — есть и там и там.
    Тип взаимодействия несёт код в ячейке, а не цвет.
    """
    output_path = Path(output_path)
    if plt is None:
        return _placeholder(output_path, "matplotlib not installed")
    if not matrix.poses or not matrix.residues:
        return _placeholder(output_path, "нет контактов для карты")

    import matplotlib.patches as mpatches

    pose_labels = [pose.label for pose in matrix.poses]
    reference_label = matrix.reference.label if matrix.reference else None

    def is_baseline(pose) -> bool:
        baseline = matrix.baseline_for(pose)
        return baseline is None or baseline == pose

    # Колонки отбираются по информативности: сначала остатки, по которым
    # контакт меняется, затем остатки позы отсчёта. Остаток, который виден
    # только у одной позы чужого кармана, растягивает карту и ничего не
    # добавляет к сравнению.
    changed = {change.residue_label for change in matrix.changes}
    in_baselines = {
        residue
        for pose in matrix.poses
        if is_baseline(pose)
        for residue in matrix.residues
        if matrix.types_at(pose, residue)
    }
    frequency = {
        residue: sum(1 for pose in matrix.poses if matrix.types_at(pose, residue))
        for residue in matrix.residues
    }
    ranked = sorted(
        matrix.residues,
        key=lambda residue: (
            residue not in changed,
            residue not in in_baselines,
            -frequency[residue],
            residue,
        ),
    )
    residues = sorted(ranked[:max_residues], key=lambda item: (
        "".join(ch for ch in item if not ch.isdigit()),
        int("".join(ch for ch in item if ch.isdigit()) or 0),
    ))
    hidden = len(matrix.residues) - len(residues)

    def baseline_types(pose) -> dict[str, set[str]]:
        """Отсчёт для позы — первая поза её собственного кармана."""
        baseline = matrix.baseline_for(pose)
        if baseline is None:
            return {}
        return {residue: matrix.types_at(baseline, residue) for residue in residues}

    width = max(6.0, 0.46 * len(residues) + 3.2)
    height = max(2.8, 0.36 * len(pose_labels) + 2.6)
    fig, ax = plt.subplots(figsize=(width, height), facecolor=_SURFACE)

    for row, pose in enumerate(matrix.poses):
        pose_baseline = baseline_types(pose)
        for column, residue in enumerate(residues):
            types = matrix.types_at(pose, residue)
            before = pose_baseline.get(residue, set())
            if not types:
                continue
            if is_baseline(pose):
                fill, ink = _DIVERGING_NEUTRAL, _INK_PRIMARY
            elif types - before:
                fill, ink = _DIVERGING_GAINED, "#ffffff"
            elif before - types:
                fill, ink = _DIVERGING_LOST, "#ffffff"
            else:
                fill, ink = _DIVERGING_NEUTRAL, _INK_PRIMARY
            # Зазор 2 px между заливками: соседние ячейки не сливаются.
            ax.add_patch(
                mpatches.FancyBboxPatch(
                    (column - 0.44, row - 0.40),
                    0.88,
                    0.80,
                    boxstyle="round,pad=0,rounding_size=0.08",
                    facecolor=fill,
                    edgecolor=_SURFACE,
                    linewidth=1.4,
                )
            )
            codes = "".join(
                INTERACTION_CODES.get(name, "?")
                for name in sorted(types, key=lambda item: INTERACTION_CODES.get(item, "?"))
            )
            ax.text(column, row, codes, ha="center", va="center", fontsize=8, color=ink)

    # Остаток, контакт с которым исчез целиком, иначе не получил бы ячейки.
    for row, pose in enumerate(matrix.poses):
        if is_baseline(pose):
            continue
        pose_baseline = baseline_types(pose)
        for column, residue in enumerate(residues):
            if matrix.types_at(pose, residue) or not pose_baseline.get(residue):
                continue
            ax.add_patch(
                mpatches.FancyBboxPatch(
                    (column - 0.44, row - 0.40),
                    0.88,
                    0.80,
                    boxstyle="round,pad=0,rounding_size=0.08",
                    facecolor=_DIVERGING_LOST,
                    edgecolor=_SURFACE,
                    linewidth=1.4,
                    alpha=0.35,
                )
            )
            ax.text(column, row, "—", ha="center", va="center", fontsize=9, color=_INK_PRIMARY)

    ax.set_xlim(-0.6, len(residues) - 0.4)
    ax.set_ylim(len(pose_labels) - 0.5, -0.5)
    ax.set_xticks(range(len(residues)), labels=residues, rotation=45, ha="right")
    tick_labels = []
    for pose in matrix.poses:
        if pose.label == reference_label:
            tick_labels.append(f"{pose.label} · поза №1")
        elif is_baseline(pose):
            tick_labels.append(f"{pose.label} · отсчёт кармана")
        else:
            tick_labels.append(pose.label)
    ax.set_yticks(range(len(pose_labels)), labels=tick_labels)
    ax.set_title(
        "Контакты с белком по позам: что появляется и что исчезает\n"
        "относительно первой позы того же кармана",
        fontsize=12,
        pad=12,
    )
    _style_axes(ax)

    handles = [
        mpatches.Patch(facecolor=_DIVERGING_GAINED, label="контакт появился"),
        mpatches.Patch(facecolor=_DIVERGING_NEUTRAL, label="есть и в отсчёте кармана"),
        mpatches.Patch(facecolor=_DIVERGING_LOST, label="контакт исчез"),
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.012, 0.055),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    for text in legend.get_texts():
        text.set_color(_INK_SECONDARY)
    footnote = _CODE_LEGEND
    if hidden > 0:
        footnote += f"\nПоказаны {len(residues)} остатков из {len(residues) + hidden}; остальные есть в pose_contacts.csv."
    fig.text(0.012, 0.012, footnote, fontsize=7.5, color=_INK_MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(output_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    return output_path


def write_distance_series(series, pose_labels: list[str], output_path: Path) -> Path:
    """Одни и те же пары атомов во всех позах, Å.

    Цвет закреплён за парой атомов по порядку слотов. Подписи у конца линии
    дублируют легенду, чтобы ряд опознавался не только по цвету.
    """
    output_path = Path(output_path)
    if plt is None:
        return _placeholder(output_path, "matplotlib not installed")
    if not series or not pose_labels:
        return _placeholder(output_path, "нет рядов расстояний")

    fig, ax = plt.subplots(figsize=(max(7.0, 0.95 * len(pose_labels) + 4.0), 4.8), facecolor=_SURFACE)
    positions = list(range(len(pose_labels)))
    drawn = series[: len(_SERIES)]

    ends: list[tuple[float, str, str]] = []
    for index, entry in enumerate(drawn):
        color = _SERIES[index]
        values = [entry.values.get(label) for label in pose_labels]
        # Разрыв там, где атома нет: пропуск нельзя рисовать нулём.
        xs = [x for x, value in zip(positions, values) if value is not None]
        ys = [value for value in values if value is not None]
        if not xs:
            continue
        ax.plot(xs, ys, color=color, linewidth=2.0, marker="o", markersize=5, label=entry.label)
        ends.append((ys[-1], entry.residue_label, color))

    # Подписи у конца линии: раздвигаются по вертикали, иначе при сходящихся
    # кривых три названия остатка печатаются друг поверх друга.
    if ends:
        low = min(value for value, _, _ in ends)
        high = max(value for value, _, _ in ends)
        gap = max((high - low) * 0.09, 0.4)
        placed: list[float] = []
        for value, label, color in sorted(ends, key=lambda item: item[0]):
            # Один сдвиг вместо цикла «пока ближе gap»: после position = prev + gap
            # разность равна gap лишь с точностью до плавающей запятой и могла
            # оказаться чуть меньше, из-за чего цикл не заканчивался никогда.
            position = max(value, placed[-1] + gap) if placed else value
            placed.append(position)
            if abs(position - value) > 1e-9:
                # Выноска к точке, от которой подпись пришлось сдвинуть.
                ax.plot(
                    [positions[-1], positions[-1] + 0.15],
                    [value, position],
                    color=color,
                    linewidth=0.8,
                    alpha=0.6,
                    clip_on=False,
                )
            ax.text(
                positions[-1] + 0.18,
                position,
                label,
                fontsize=9,
                color=color,
                va="center",
                ha="left",
            )

    ax.set_xlim(-0.4, len(pose_labels) - 1 + max(1.4, 0.1 * len(pose_labels)))
    ax.set_xticks(positions, labels=pose_labels, rotation=30, ha="right")
    ax.set_ylabel("Расстояние, Å", color=_INK_SECONDARY, fontsize=10)
    ax.set_title("Расстояния до ключевых остатков по позам", fontsize=12, pad=12)
    ax.yaxis.grid(True, color=_GRID, linewidth=1.0)
    ax.set_axisbelow(True)
    _style_axes(ax)

    legend = fig.legend(loc="lower left", bbox_to_anchor=(0.012, 0.01), ncol=3, frameon=False, fontsize=8)
    for text in legend.get_texts():
        text.set_color(_INK_SECONDARY)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(output_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)
    return output_path
