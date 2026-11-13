"""Тесты по умолчанию не ходят в Ollama."""
import os

os.environ.setdefault("DRUGDIS_LLM_PROVIDER", "keyword")
