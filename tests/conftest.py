"""Тесты по умолчанию не ходят в Ollama."""
import os

os.environ.setdefault("DEARMO_LLM_PROVIDER", "keyword")
