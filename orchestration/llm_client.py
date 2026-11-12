"""Локальный Qwen через Ollama. Если сервер недоступен, планировщик берёт ключевые слова."""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any

import httpx

DEFAULT_OLLAMA_MODEL = "qwen2.5:7b-instruct"

PLANNER_JSON_SCHEMA_DESCRIPTION = """
Верни только JSON без markdown:
{
  "scenario_id": "target_for_ligand|dock_pockets|compare_docking|validate_pose",
  "ligands": ["английские названия молекул"],
  "smiles": "SMILES первой молекулы, если достоверно знаешь, иначе null",
  "target_hint": "имя белка или null",
  "pdb_id": "четырёхсимвольный PDB или null",
  "poses_csv_path": "путь или null",
  "pose_file": "путь к позе или null",
  "ligand_id": "идентификатор лиганда или null",
  "explanation_hint": "кратко что хочет пользователь"
}
"""


@dataclass
class PlannerExtraction:
    scenario_id: str | None
    ligands: list[str]
    smiles: str | None
    target_hint: str | None
    pdb_id: str | None
    poses_csv_path: str | None
    pose_file: str | None
    ligand_id: str | None
    explanation_hint: str | None
    source: str


def _extract_json_object(text: str) -> dict[str, Any] | None:
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        data = json.loads(text)
        return data if isinstance(data, dict) else None
    except json.JSONDecodeError:
        match = re.search(r"\{[\s\S]*\}", text)
        if not match:
            return None
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            return None


def _normalize_planner_data(data: dict[str, Any], *, source: str) -> PlannerExtraction:
    ligands = data.get("ligands") or []
    if isinstance(ligands, str):
        ligands = [ligands]
    if not isinstance(ligands, list):
        ligands = []
    ligands = [str(item).strip() for item in ligands if str(item).strip()]
    pdb = data.get("pdb_id")
    if isinstance(pdb, str):
        pdb = pdb.strip().upper()[:4] if pdb.strip() else None
    smiles = data.get("smiles")
    if isinstance(smiles, str):
        smiles = smiles.strip() or None
    else:
        smiles = None
    return PlannerExtraction(
        scenario_id=data.get("scenario_id"),
        ligands=ligands,
        smiles=smiles,
        target_hint=data.get("target_hint"),
        pdb_id=pdb,
        poses_csv_path=data.get("poses_csv_path"),
        pose_file=data.get("pose_file"),
        ligand_id=data.get("ligand_id"),
        explanation_hint=data.get("explanation_hint"),
        source=source,
    )


def provider() -> str:
    return os.environ.get("DEARMO_LLM_PROVIDER", "auto").strip().lower()


def ollama_host() -> str:
    return os.environ.get("OLLAMA_HOST", "http://127.0.0.1:11434").rstrip("/")


def ollama_model() -> str:
    return os.environ.get("OLLAMA_MODEL", DEFAULT_OLLAMA_MODEL)


def _timeout() -> float:
    return float(os.environ.get("OLLAMA_TIMEOUT_S", "120"))


def _debug(message: str) -> None:
    if os.environ.get("DEARMO_LLM_DEBUG", "").lower() in ("1", "true", "yes"):
        print(f"[llm] {message}", flush=True)


def ollama_available() -> bool:
    try:
        response = httpx.get(f"{ollama_host()}/api/tags", timeout=0.4)
        return response.status_code == 200
    except httpx.HTTPError:
        return False


def chat_completion(system: str, user: str, *, temperature: float = 0.1) -> str | None:
    payload: dict[str, Any] = {
        "model": ollama_model(),
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "temperature": temperature,
        "stream": False,
    }
    try:
        response = httpx.post(
            f"{ollama_host()}/v1/chat/completions",
            json=payload,
            timeout=_timeout(),
        )
        response.raise_for_status()
        data = response.json()
        return data["choices"][0]["message"]["content"]
    except (httpx.HTTPError, KeyError, IndexError, json.JSONDecodeError) as exc:
        _debug(f"ollama chat failed: {exc}")
        return None


def ask_json(system: str, user: str) -> dict[str, Any] | None:
    """Короткий JSON от локальной Qwen. При keyword-режиме и недоступном Ollama возвращает None."""
    mode = provider()
    if mode == "keyword":
        return None
    if mode in ("auto", "ollama") and not ollama_available():
        return None
    content = chat_completion(system, user)
    if not content:
        return None
    data = _extract_json_object(content)
    if data:
        return data
    repair = chat_completion("Верни только исправленный JSON-объект без markdown.", content)
    return _extract_json_object(repair or "")


@dataclass
class Decision:
    """Выбор LLM из списка, который собрали инструменты. Числа сюда не попадают."""

    index: int
    rationale: str
    source: str
    option: dict[str, Any]


def llm_available() -> bool:
    mode = provider()
    if mode == "keyword":
        return False
    return ollama_available()


def decide(
    task: str,
    options: list[dict[str, Any]],
    *,
    criteria: str,
    context: str = "",
) -> Decision | None:
    """LLM выбирает ОДИН вариант из готового списка и обосновывает выбор.

    Модель не может ни придумать вариант, ни вернуть число: допустим только индекс
    в пределах списка. Всё, что вне диапазона, отбрасывается, и агент работает
    своей детерминированной эвристикой.
    """
    if not options or not llm_available():
        return None
    lines = []
    for index, option in enumerate(options):
        summary = "; ".join(
            f"{key}={value}" for key, value in option.items() if value not in (None, "", [], {})
        )
        lines.append(f"[{index}] {summary}")
    system = (
        "Ты научный агент drug discovery. Выбери ОДИН вариант из списка.\n"
        "Запрещено: придумывать варианты, менять числа, называть источники, которых нет в списке.\n"
        'Ответ строго JSON: {"index": <целое>, "rationale": "<1-2 предложения по-русски>"}'
    )
    user = (
        f"Задача: {task}\n"
        f"Критерии выбора: {criteria}\n"
        + (f"Контекст запроса пользователя: {context}\n" if context else "")
        + f"Варианты:\n" + "\n".join(lines)
    )
    data = ask_json(system, user)
    if not data:
        return None
    raw_index = data.get("index")
    try:
        index = int(raw_index)
    except (TypeError, ValueError):
        _debug(f"decide: нечисловой index {raw_index!r}")
        return None
    if not 0 <= index < len(options):
        _debug(f"decide: index {index} вне диапазона 0..{len(options) - 1}")
        return None
    rationale = str(data.get("rationale") or "").strip()
    return Decision(
        index=index,
        rationale=rationale or "Модель не дала обоснования.",
        source=f"ollama:{ollama_model()}",
        option=options[index],
    )


def critique(question: str, facts: str) -> dict[str, Any] | None:
    """LLM-критик: хватает ли собранного, что доработать. Возвращает решение, а не числа."""
    if not llm_available():
        return None
    system = (
        "Ты научный рецензент. По собранным фактам реши, отвечают ли они на вопрос пользователя.\n"
        "Не добавляй новых чисел и источников. Оценивай только полноту уже собранного.\n"
        'Ответ строго JSON: {"sufficient": true|false, "missing": "<чего не хватает или пусто>", '
        '"next_scenario": "target_for_ligand|dock_pockets|compare_docking|validate_pose|null"}'
    )
    data = ask_json(system, f"Вопрос пользователя: {question}\n\nСобранные факты:\n{facts}")
    if not data:
        return None
    return {
        "sufficient": bool(data.get("sufficient", True)),
        "missing": str(data.get("missing") or "").strip(),
        "next_scenario": data.get("next_scenario") if data.get("next_scenario") not in ("null", "") else None,
        "source": f"ollama:{ollama_model()}",
    }


def plan_with_ollama(user_text: str, system_prompt: str) -> PlannerExtraction | None:
    mode = provider()
    if mode == "keyword":
        return None
    if mode == "auto" and not ollama_available():
        return None
    if mode == "ollama" and not ollama_available():
        _debug("ollama requested but not reachable")
        return None
    content = chat_completion(
        system_prompt + "\n\n" + PLANNER_JSON_SCHEMA_DESCRIPTION,
        user_text,
    )
    if not content:
        return None
    data = _extract_json_object(content)
    if not data:
        repair = chat_completion("Верни только исправленный JSON-объект без markdown.", content)
        data = _extract_json_object(repair or "")
    if not data:
        return None
    return _normalize_planner_data(data, source=f"ollama:{ollama_model()}")


def plan_with_llm(user_text: str, system_prompt: str) -> PlannerExtraction | None:
    return plan_with_ollama(user_text, system_prompt)


_GLOSSARY = """
- molecular_weight — молекулярная масса, г/моль
- logp — липофильность logP
- tpsa — топологическая площадь полярной поверхности, Å²
- hbd — число доноров водородной связи
- hba — число акцепторов водородной связи
- rotatable_bonds — число вращаемых связей
- qed — QED, оценка drug-likeness от 0 до 1
- pubchem_cid — идентификатор соединения в PubChem
- vina_score — оценка AutoDock Vina, ккал/моль (расчёт, не измеренная аффинность)
- rmsd_to_crystal — RMSD до кристаллического лиганда, Å
- checks_passed_percent — доля пройденных геометрических проверок позы
- pdb_structure — код структуры в RCSB PDB
- resolution — разрешение структуры, Å
- Ki, Kd, IC50 — экспериментальные константы связывания или ингибирования, нМ
"""


def narrate_with_llm(evidence_summary: str, limitations: list[str]) -> str | None:
    mode = provider()
    if mode == "keyword":
        return None
    if mode == "auto" and not ollama_available():
        return None
    text = chat_completion(
        "Сформируй краткое объяснение на русском по фактам (3-5 предложений).\n"
        "Не добавляй чисел, ссылок и утверждений, которых нет в фактах.\n"
        "Используй принятые русские термины и не переводи их заново:\n"
        f"{_GLOSSARY}\n"
        "Если термина нет в списке, оставь английское название как есть.",
        f"Факты:\n{evidence_summary}\n\nОграничения:\n" + "\n".join(limitations),
        temperature=0.2,
    )
    return text.strip() if text else None
