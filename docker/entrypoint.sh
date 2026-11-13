#!/bin/sh
# Точка входа контейнера DRUGDIS.
#
#   doctor            самопроверка окружения (Python, RDKit, Vina, P2Rank, Ollama, сеть)
#   doctor --deep     то же плюс настоящий прогон P2Rank и Vina на встроенной структуре
#   demo | ask ...    команды cli.py; перед ними ждём Qwen
#   pytest | sh ...   запускаются как есть
#
# Если выбран провайдер ollama, а модель недоступна, контейнер СТОИТ и объясняет почему.
# Молча уйти на ключевые слова нельзя: пользователь решил бы, что ответила Qwen.
set -eu

case "${1:-}" in
    doctor) shift; exec python -m tools.doctor "$@" ;;
    pytest|python|sh|bash) exec "$@" ;;
esac

wait_for_llm() {
    host="${OLLAMA_HOST:-http://ollama:11434}"
    model="${OLLAMA_MODEL:-qwen2.5:7b-instruct}"
    limit="${OLLAMA_WAIT_S:-600}"
    waited=0
    echo "[drugdis] жду Ollama ($host) и модель $model, не дольше ${limit} с..."
    while :; do
        status="$(python - "$host" "$model" <<'PY' 2>/dev/null || echo down
import json, sys, urllib.request
host, model = sys.argv[1], sys.argv[2]
try:
    with urllib.request.urlopen(host.rstrip("/") + "/api/tags", timeout=3) as resp:
        names = {m.get("name") for m in json.load(resp).get("models", [])}
except Exception:
    print("down"); sys.exit(0)
print("ready" if model in names else "no-model")
PY
)"
        case "$status" in
            ready) echo "[drugdis] Qwen готова: $model"; return 0 ;;
            no-model) hint="Ollama отвечает, но модели $model нет. Скачайте её: docker compose run --rm model-pull" ;;
            *) hint="Ollama ($host) не отвечает. Запущен ли контейнер: docker compose up -d ollama" ;;
        esac
        if [ "$waited" -ge "$limit" ]; then
            echo "[drugdis] ОШИБКА: языковая модель недоступна. $hint" >&2
            if [ "${DRUGDIS_ALLOW_FALLBACK:-0}" = "1" ]; then
                echo "[drugdis] DRUGDIS_ALLOW_FALLBACK=1: продолжаю БЕЗ LLM (провайдер keyword). Ответ модели в этом прогоне нет." >&2
                export DRUGDIS_LLM_PROVIDER=keyword
                return 0
            fi
            echo "[drugdis] Чтобы работать без модели явно: DRUGDIS_LLM_PROVIDER=keyword или DRUGDIS_ALLOW_FALLBACK=1" >&2
            exit 3
        fi
        [ $((waited % 30)) -eq 0 ] && echo "[drugdis] ...$hint"
        sleep 5
        waited=$((waited + 5))
    done
}

if [ "${DRUGDIS_LLM_PROVIDER:-auto}" = "ollama" ]; then
    wait_for_llm
fi

exec python cli.py "$@"
