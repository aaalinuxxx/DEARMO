"""Самопроверка окружения: всё ли на месте, чтобы сценарии реально считались.

    python -m tools.doctor           # быстрая проверка (секунды)
    python -m tools.doctor --deep    # плюс настоящий прогон P2Rank и Vina на встроенной структуре
    docker compose run --rm dearmo doctor

Зачем она нужна. Без неё сломанное окружение проявляется как странный ответ посреди
сценария: нет Vina — «докинг не выполнен», нет Java — «P2Rank не нашёл кармана», нет модели —
тихий переход на ключевые слова. Здесь каждая причина названа сразу и по имени.

Код выхода 1, если не работает что-то обязательное. Недоступность внешних баз — предупреждение,
а не ошибка: ChEMBL, например, иногда отвечает HTTP 500, и система это переносит.
"""
from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

OK, WARN, FAIL = "ok", "warn", "fail"
_MARK = {OK: "[ OK ]", WARN: "[WARN]", FAIL: "[FAIL]"}
_results: list[tuple[str, str, str]] = []


def _report(level: str, name: str, detail: str = "") -> None:
    _results.append((level, name, detail))
    print(f"{_MARK[level]} {name}" + (f": {detail}" if detail else ""), flush=True)


def _run(cmd: list[str], timeout: float = 30) -> tuple[int, str]:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 1, str(exc)
    return proc.returncode, ((proc.stdout or "") + (proc.stderr or "")).strip()


def check_python_stack() -> None:
    _report(OK, "Python", sys.version.split()[0])
    for module in ("rdkit", "meeko", "gemmi", "numpy", "scipy", "pandas", "matplotlib", "httpx"):
        try:
            mod = importlib.import_module(module)
        except Exception as exc:  # noqa: BLE001 — любая причина отказа важна
            _report(FAIL, module, f"не импортируется ({exc}). pip install -r requirements.txt")
            continue
        _report(OK, module, str(getattr(mod, "__version__", "")))

    pdb2pqr = Path(sys.executable).parent / "pdb2pqr"
    found = str(pdb2pqr) if pdb2pqr.is_file() else shutil.which("pdb2pqr")
    if found:
        _report(OK, "pdb2pqr", found)
    else:
        _report(
            WARN,
            "pdb2pqr",
            "не найден: рецептор будет подготовлен упрощённо, без водородов, качество докинга ниже",
        )


def check_vina() -> None:
    from tools.vina_runner import vina_binary

    binary = vina_binary()
    if not binary:
        _report(FAIL, "AutoDock Vina", "бинарник не найден (VINA_BIN или tools/bin/vina)")
        return
    code, out = _run([binary, "--version"])
    if code == 0 and "Vina" in out:
        _report(OK, "AutoDock Vina", f"{out.splitlines()[0]} ({binary})")
    else:
        _report(FAIL, "AutoDock Vina", f"{binary} не запускается: {out[:160]}")


def check_p2rank() -> None:
    from tools.p2rank_runner import bundled_java_home, p2rank_binary

    binary = p2rank_binary()
    java_home = os.environ.get("JAVA_HOME") or bundled_java_home()
    java = str(Path(java_home) / "bin" / "java") if java_home else shutil.which("java")
    if not java:
        _report(FAIL, "Java", "не найдена: P2Rank без неё не работает (нужна версия 21)")
    else:
        code, out = _run([java, "-version"])
        first = out.splitlines()[0] if out else ""
        _report(OK if code == 0 else FAIL, "Java", first or "не запускается")
        if code == 0 and '"21' not in first and '"17' not in first:
            _report(
                WARN,
                "Версия Java",
                "P2Rank 2.5 собран на Groovy 4 и на JDK новее 21 падает; рекомендуется 21",
            )
    if binary:
        _report(OK, "P2Rank", binary)
    else:
        _report(FAIL, "P2Rank", "не найден (P2RANK_BIN или tools/bin/p2rank/prank)")


def check_data() -> None:
    required = _ROOT / "data" / "external" / "rcsb" / "3PTB.pdb"
    if required.is_file():
        _report(OK, "Опорные данные", "data/external/rcsb")
    else:
        _report(WARN, "Опорные данные", "data/external/rcsb пуст: структуры будут скачиваться из RCSB")


def check_llm() -> None:
    from orchestration import llm_client

    provider = llm_client.provider()
    host = llm_client.ollama_host()
    model = llm_client.ollama_model()
    if provider == "keyword":
        _report(WARN, "Языковая модель", "провайдер keyword: Qwen не используется, планировщик работает по ключевым словам")
        return
    required = provider == "ollama"
    level_missing = FAIL if required else WARN

    try:
        with urllib.request.urlopen(host.rstrip("/") + "/api/tags", timeout=4) as resp:
            names = {m.get("name") for m in json.load(resp).get("models", [])}
    except Exception as exc:  # noqa: BLE001
        _report(level_missing, "Ollama", f"{host} не отвечает ({exc.__class__.__name__}). docker compose up -d ollama")
        return
    _report(OK, "Ollama", host)
    if model not in names:
        _report(level_missing, f"Модель {model}", "не скачана. docker compose run --rm model-pull")
        return
    _report(OK, f"Модель {model}", "скачана")

    # Модель может быть на диске, но не загружаться: не хватает памяти. Узнать это можно только запросом.
    started = time.perf_counter()
    answer = llm_client.chat_completion("Отвечай одним словом.", "Скажи «готово».", temperature=0.0)
    elapsed = time.perf_counter() - started
    if answer:
        _report(OK, "Генерация Qwen", f"ответила за {elapsed:.0f} с: «{answer.strip()[:40]}»")
    else:
        _report(
            level_missing,
            "Генерация Qwen",
            "модель не ответила: вероятно, не хватает памяти (для qwen2.5:7b нужно около 6 ГБ свободной RAM)",
        )


def check_network() -> None:
    probes = {
        "PubChem": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularWeight/JSON",
        "RCSB PDB": "https://data.rcsb.org/rest/v1/core/entry/3PTB",
        "UniProt": "https://rest.uniprot.org/uniprotkb/P00918?fields=accession&format=json",
        "BindingDB": "https://www.bindingdb.org/rwd/bind/index.jsp",
        "ChEMBL": "https://www.ebi.ac.uk/chembl/api/data/status.json",
    }
    for name, url in probes.items():
        request = urllib.request.Request(url, headers={"User-Agent": "DEARMO-doctor/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=12) as resp:
                _report(OK, f"Сеть: {name}", f"HTTP {resp.status}")
        except Exception as exc:  # noqa: BLE001
            note = " (источник нестабилен, система это учитывает)" if name == "ChEMBL" else ""
            _report(WARN, f"Сеть: {name}", f"недоступен: {str(exc)[:80]}{note}")


def _first_chain(pdb_path: Path) -> str | None:
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") and len(line) > 21 and line[21].strip():
            return line[21]
    return None


def check_deep() -> None:
    """Настоящий расчёт: P2Rank и Vina на встроенной структуре, без сети."""
    from tools.live_docking import dock_pair

    structure = _ROOT / "data" / "external" / "rcsb" / "3PTB.pdb"
    if not structure.is_file():
        _report(WARN, "Глубокая проверка", "нет встроенной структуры 3PTB, пропущена")
        return
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="dearmo-doctor-") as tmp:
        try:
            # SMILES бензамидина передаётся прямо: проверка не должна зависеть от PubChem.
            # Цепь берётся из самой структуры (первая цепь с белковыми атомами), а не из памяти.
            run = dock_pair("3PTB", "NC(=N)c1ccccc1", "benzamidine", Path(tmp), chain=_first_chain(structure), crystal_comp_id="BEN")
        except Exception as exc:  # noqa: BLE001
            _report(FAIL, "Глубокая проверка", f"докинг не выполнен: {str(exc)[:200]}")
            return
        elapsed = time.perf_counter() - started
        if not run.pockets:
            _report(FAIL, "P2Rank", "карманов не найдено")
        else:
            _report(OK, "P2Rank (расчёт)", f"найдено карманов: {len(run.pockets)}")
        score = run.best.vina_score
        if -20.0 < score < 0.0:
            _report(OK, "AutoDock Vina (расчёт)", f"лучшая поза {score} ккал/моль, {len(run.poses)} поз, {elapsed:.0f} с")
        else:
            _report(FAIL, "AutoDock Vina (расчёт)", f"нефизичная оценка {score}")
        rmsd = (run.control or {}).get("best_pose_rmsd_to_crystal_a")
        if rmsd is not None:
            _report(
                OK if rmsd <= 2.0 else WARN,
                "Контроль редокинга 3PTB",
                f"RMSD позы №1 до кристалла {rmsd:.2f} Å (порог 2 Å)",
            )


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    deep = "--deep" in args
    offline = "--offline" in args

    print("DEARMO doctor: проверка окружения\n")
    check_python_stack()
    check_vina()
    check_p2rank()
    check_data()
    check_llm()
    if not offline:
        check_network()
    if deep:
        check_deep()

    failed = [name for level, name, _ in _results if level == FAIL]
    warned = [name for level, name, _ in _results if level == WARN]
    print()
    if failed:
        print(f"Итог: НЕ ГОТОВО. Не работает: {', '.join(failed)}.")
        return 1
    print("Итог: готово" + (f", есть предупреждения ({len(warned)})." if warned else "."))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
