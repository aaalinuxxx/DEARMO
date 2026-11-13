.DEFAULT_GOAL := help
PY ?= .venv/bin/python

.PHONY: help venv install test test-fast demo demo-offline fetch-data lint clean check redock scenario1 scenario2

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

scenario1: ## Сценарий 1: мишень для ацетазоламида (нужна сеть и Ollama)
	DEARMO_LLM_PROVIDER=ollama $(PY) cli.py ask "Подбери белок-мишень для ацетазоламида и обоснуй выбор экспериментальными данными."

scenario2: ## Сценарий 2: P2Rank + докинг в топ-3 кармана для эрлотиниба в 1M17
	DEARMO_LLM_PROVIDER=keyword $(PY) cli.py ask "Найди три наиболее перспективных кармана, проведи докинг эрлотиниба в 1M17 в каждом и покажи лучшую позу по Vina."

redock: ## Редокинг-бенчмарк: доля поз №1 в пределах 2 Å от кристалла, пишет docs/redocking.md
	$(PY) -m tools.redock_benchmark

fetch-data: ## Скачать опорный набор данных в data/external/
	$(PY) cli.py fetch-data

check: test demo-offline ## Полная проверка перед коммитом

clean: ## Убрать кеши и результаты прогонов
	rm -rf outputs .pytest_cache
	find . -name __pycache__ -type d -not -path "./.venv/*" -exec rm -rf {} +
