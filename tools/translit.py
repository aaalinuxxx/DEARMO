"""Русское название молекулы → латиница, с проверкой по PubChem.

Это орфографика, а не химические знания: никакого списка молекул здесь нет.
Любой вариант написания считается верным, только если PubChem **точно** нашёл
по нему соединение. Нечёткие подсказки сюда не попадают: для «целекоксиб»
fuzzy-поиск PubChem предлагает Trotabresib, то есть другое вещество, и принять
такую подстановку означало бы выдумать результат.
"""
from __future__ import annotations

import re

_BASE = {
    "а": "a", "б": "b", "в": "v", "г": "g", "д": "d", "е": "e", "ё": "e", "ж": "zh",
    "з": "z", "и": "i", "й": "y", "л": "l", "м": "m", "н": "n", "о": "o", "п": "p",
    "р": "r", "с": "s", "т": "t", "у": "u", "ф": "f", "х": "kh", "ч": "ch",
    "ш": "sh", "щ": "sch", "ъ": "", "ы": "y", "ь": "", "э": "e", "ю": "yu", "я": "ya",
}
# Буквы, у которых в химической номенклатуре два равноправных латинских написания.
_AMBIGUOUS = {"ц": ("c", "ts"), "к": ("k", "c")}

# Падежные окончания: «эрлотиниба» -> «эрлотиниб».
_ENDINGS = ("ом", "ами", "ах", "ов", "а", "у", "е", "ы", "и", "ю")

_CYRILLIC = re.compile(r"[а-яё]", re.I)


def looks_russian(word: str) -> bool:
    return bool(_CYRILLIC.search(word))


def _transliterate(word: str, *, ts: str, k: str) -> str:
    out = []
    for char in word.lower():
        if char == "ц":
            out.append(ts)
        elif char == "к":
            out.append(k)
        else:
            out.append(_BASE.get(char, char))
    return "".join(out)


def candidates(word: str, *, limit: int = 14) -> list[str]:
    """Варианты латинского написания в порядке убывания правдоподобия."""
    word = word.strip().lower()
    if not word:
        return []

    stems = [word]
    for ending in _ENDINGS:
        if word.endswith(ending) and len(word) - len(ending) >= 5:
            stems.append(word[: -len(ending)])

    # «кс» в латинских названиях часто пишется как x (целекоксиб -> celecoxib).
    # Заменять нужно до транслитерации: иначе вариант «к -> c» превратит «кс» в «cs».
    spellings: list[str] = []
    for stem in stems:
        spellings.append(stem)
        if "кс" in stem:
            spellings.append(stem.replace("кс", "x"))

    per_spelling: list[list[str]] = []
    for stem in spellings:
        forms: list[str] = []
        for ts in _AMBIGUOUS["ц"]:
            for k in _AMBIGUOUS["к"]:
                base = _transliterate(stem, ts=ts, k=k)
                group = [base]
                if "ks" in base:
                    group.append(base.replace("ks", "x", 1))
                # «ацетазоламид» -> acetazolamide: латинские названия на -amid/-id
                # в английском получают конечную -e.
                for form in list(group):
                    if form.endswith(("amid", "id", "in", "ol")) and not form.endswith("e"):
                        group.append(form + "e")
                forms.extend(form for form in group if form and form not in forms)
        per_spelling.append(forms)

    # Варианты перемешиваются по кругу: иначе лимит съедают написания одного
    # падежного стебля, и до «целекоксиб -> celecoxib» очередь не доходит.
    out: list[str] = []
    for index in range(max((len(group) for group in per_spelling), default=0)):
        for group in per_spelling:
            if index < len(group) and group[index] not in out:
                out.append(group[index])
    return out[:limit]


def resolve_russian_name(word: str, resolver) -> dict | None:
    """Первый вариант написания, который PubChem подтвердил точным совпадением.

    `resolver` — функция вида `resolve_compound(name) -> dict`; передаётся
    параметром, чтобы модуль не зависел от сети в тестах.
    """
    for variant in candidates(word):
        info = resolver(variant)
        if info and info.get("cid") and not info.get("ambiguous"):
            return {**info, "translit_from": word, "translit_variant": variant}
    return None
