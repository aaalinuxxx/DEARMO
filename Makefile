.DEFAULT_GOAL := help
PY ?= .venv/bin/python

.PHONY: help venv install test test-fast demo demo-offline fetch-data lint clean check

help: ## Показать список команд
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

venv: ## Создать виртуальное окружение
	python3 -m venv .venv
	$(PY) -m pip install --upgrade pip

install: ## Поставить зависимости
	$(PY) -m pip install -r requirements.txt

test: ## Все тесты
	$(PY) -m pytest -q

test-fast: ## Тесты без интеграционных
	$(PY) -m pytest -q -m "not integration"

demo: ## Прогнать сценарии и обновить docs/scenarios.md
	DEARMO_LLM_PROVIDER=ollama $(PY) cli.py demo

demo-offline: ## Прогнать сценарии без LLM, не трогая документацию
	DEARMO_LLM_PROVIDER=keyword $(PY) cli.py demo --no-docs

fetch-data: ## Скачать опорный набор данных в data/external/
	$(PY) cli.py fetch-data

check: test demo-offline ## Полная проверка перед коммитом

clean: ## Убрать кеши и результаты прогонов
	rm -rf outputs .pytest_cache
	find . -name __pycache__ -type d -not -path "./.venv/*" -exec rm -rf {} +
