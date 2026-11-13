#!/bin/sh
# Скачать Qwen после старта контейнера Ollama.
set -eu
ollama pull "${OLLAMA_MODEL:-qwen2.5:7b-instruct}"
