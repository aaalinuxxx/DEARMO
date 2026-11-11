"""Клиент GigaChat: при отсутствии ключа планировщик использует keyword fallback."""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any

PLANNER_JSON_SCHEMA_DESCRIPTION = """
Верни только JSON без markdown:
{
  "scenario_id": "target_for_ligand|dock_pockets|compare_docking|validate_pose",
  "ligands": ["список имён или SMILES"],
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
    target_hint: str | None
    pdb_id: str | None
    poses_csv_path: str | None
    pose_file: str | None
    ligand_id: str | None
    explanation_hint: str | None
    source: str  # gigachat | keyword


def _extract_json_object(text: str) -> dict[str, Any] | None:
    text = text.strip()
  # strip code fences
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        data = json.loads(text)
        return data if isinstance(data, dict) else None
    except json.JSONDecodeError:
        match = re.search(r"\{[\s\S]*\}", text)
        if match:
            try:
                return json.loads(match.group(0))
            except json.JSONDecodeError:
                return None
    return None


def _normalize_planner_data(data: dict[str, Any]) -> PlannerExtraction:
    ligands = data.get("ligands") or []
    if isinstance(ligands, str):
        ligands = [ligands]
    if not isinstance(ligands, list):
        ligands = []
    ligands = [str(x).strip() for x in ligands if str(x).strip()]
    pdb = data.get("pdb_id")
    if isinstance(pdb, str):
        pdb = pdb.strip().upper()[:4] if pdb.strip() else None
    return PlannerExtraction(
        scenario_id=data.get("scenario_id"),
        ligands=ligands,
        target_hint=data.get("target_hint"),
        pdb_id=pdb,
        poses_csv_path=data.get("poses_csv_path"),
        pose_file=data.get("pose_file"),
        ligand_id=data.get("ligand_id"),
        explanation_hint=data.get("explanation_hint"),
        source="gigachat",
    )


def gigachat_available() -> bool:
    return bool(os.environ.get("GIGACHAT_CREDENTIALS") or os.environ.get("GIGACHAT_API_KEY"))


def plan_with_gigachat(user_text: str, system_prompt: str) -> PlannerExtraction | None:
    if not gigachat_available():
        return None
    try:
        from gigachat import GigaChat
    except ImportError:
        return None
    credentials = os.environ.get("GIGACHAT_CREDENTIALS") or os.environ.get("GIGACHAT_API_KEY")
    model = os.environ.get("GIGACHAT_MODEL", "GigaChat")
    verify = os.environ.get("GIGACHAT_VERIFY_SSL_CERTS", "false").lower() in ("1", "true", "yes")
    payload = f"{system_prompt}\n\n{PLANNER_JSON_SCHEMA_DESCRIPTION}\n\nЗапрос пользователя:\n{user_text}"
    try:
        with GigaChat(credentials=credentials, model=model, verify_ssl_certs=verify) as client:
            response = client.chat(payload)
            content = response.choices[0].message.content
    except Exception:
        return None
    data = _extract_json_object(content or "")
    if not data:
        return None
    return _normalize_planner_data(data)


def narrate_with_gigachat(evidence_summary: str, limitations: list[str]) -> str | None:
    if not gigachat_available():
        return None
    try:
        from gigachat import GigaChat
    except ImportError:
        return None
    credentials = os.environ.get("GIGACHAT_CREDENTIALS") or os.environ.get("GIGACHAT_API_KEY")
    model = os.environ.get("GIGACHAT_MODEL", "GigaChat")
    verify = os.environ.get("GIGACHAT_VERIFY_SSL_CERTS", "false").lower() in ("1", "true", "yes")
    prompt = (
        "Сформируй краткое объяснение на русском по фактам ниже. "
        "Не добавляй чисел, ссылок и утверждений, которых нет в фактах.\n\n"
        f"Факты:\n{evidence_summary}\n\nОграничения:\n" + "\n".join(limitations)
    )
    try:
        with GigaChat(credentials=credentials, model=model, verify_ssl_certs=verify) as client:
            response = client.chat(prompt)
            return (response.choices[0].message.content or "").strip()
    except Exception:
        return None
