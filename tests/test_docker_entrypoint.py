"""Контракт контейнера: громко сообщать о недоступной модели и не плавать по версиям.

Тесты не требуют Docker: entrypoint — обычный sh-скрипт, а зафиксированные версии
проверяются по самому requirements.txt.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
_ENTRYPOINT = _ROOT / "docker" / "entrypoint.sh"


def _run(extra_env: dict[str, str], *args: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    # entrypoint зовёт `python`: подсовываем интерпретатор, под которым идут тесты.
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    env.update(extra_env)
    return subprocess.run(
        ["sh", str(_ENTRYPOINT), *args],
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )


@pytest.mark.skipif(not _ENTRYPOINT.is_file(), reason="нет docker/entrypoint.sh")
def test_entrypoint_is_valid_shell():
    assert subprocess.run(["sh", "-n", str(_ENTRYPOINT)], capture_output=True).returncode == 0


def test_unreachable_llm_stops_the_container_instead_of_silently_falling_back():
    """Пользователь выбрал Qwen. Если её нет, ответ «как будто от Qwen» был бы обманом."""
    proc = _run(
        {
            "DRUGDIS_LLM_PROVIDER": "ollama",
            "OLLAMA_HOST": "http://127.0.0.1:9",
            "OLLAMA_WAIT_S": "0",
            "DRUGDIS_ALLOW_FALLBACK": "0",
        },
        "--help",
    )
    assert proc.returncode == 3, proc.stdout + proc.stderr
    assert "языковая модель недоступна" in proc.stderr
    assert "docker compose" in proc.stderr, "ошибка должна подсказывать, что делать"


def test_fallback_is_allowed_only_when_asked_and_is_announced():
    proc = _run(
        {
            "DRUGDIS_LLM_PROVIDER": "ollama",
            "OLLAMA_HOST": "http://127.0.0.1:9",
            "OLLAMA_WAIT_S": "0",
            "DRUGDIS_ALLOW_FALLBACK": "1",
        },
        "--help",
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "БЕЗ LLM" in proc.stderr, "переход на ключевые слова обязан быть объявлен"


def test_keyword_provider_does_not_wait_for_ollama():
    proc = _run({"DRUGDIS_LLM_PROVIDER": "keyword", "OLLAMA_HOST": "http://127.0.0.1:9"}, "--help")
    assert proc.returncode == 0
    assert "жду Ollama" not in proc.stdout + proc.stderr


def test_requirements_are_pinned():
    """Без фиксации `pip install` через неделю соберёт другое окружение, и тесты станут красными."""
    lines = [
        line.strip()
        for line in (_ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert lines
    unpinned = [line for line in lines if not re.fullmatch(r"[A-Za-z0-9_.\-]+==[0-9][A-Za-z0-9_.\-]*", line)]
    assert not unpinned, f"версии не зафиксированы: {unpinned}"


def test_docker_image_does_not_ship_host_binaries():
    """tools/bin собран под macOS: в Linux-образе он бесполезен и вытеснил бы настоящий бинарник."""
    ignore = (_ROOT / ".dockerignore").read_text(encoding="utf-8")
    assert "tools/bin/" in ignore
    dockerfile = (_ROOT / "Dockerfile").read_text(encoding="utf-8")
    assert "VINA_BIN=/usr/local/bin/vina" in dockerfile
    assert "P2RANK_BIN=/opt/p2rank/prank" in dockerfile
