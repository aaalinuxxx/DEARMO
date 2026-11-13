# syntax=docker/dockerfile:1
#
# Образ приложения DEARMO. Внутри всё, что нужно расчётам: Python-стек, AutoDock Vina,
# P2Rank и Java. Языковая модель живёт в отдельном контейнере `ollama` (см. docker-compose.yml).
#
# Образ собирается и на x86_64, и на arm64: бинарники Vina статические и берутся под
# архитектуру сборки. Бинарники из tools/bin/ (они собраны под macOS) в образ не попадают.

FROM python:3.14-slim-trixie

ARG VINA_VERSION=1.2.7
ARG P2RANK_VERSION=2.5.1

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app \
    MPLBACKEND=Agg \
    VINA_BIN=/usr/local/bin/vina \
    P2RANK_BIN=/opt/p2rank/prank \
    DEARMO_LLM_PROVIDER=ollama \
    OLLAMA_HOST=http://ollama:11434 \
    OLLAMA_MODEL=qwen2.5:7b-instruct

# Java 21 нужна P2Rank: сборка на Groovy 4 не читает классы более новых JDK.
RUN apt-get update && apt-get install -y --no-install-recommends \
        openjdk-21-jre-headless \
        curl \
        ca-certificates \
        libxrender1 \
        libxext6 \
        libsm6 \
    && rm -rf /var/lib/apt/lists/*

# AutoDock Vina: готовый статический бинарник под архитектуру сборки.
RUN set -eux; \
    case "$(dpkg --print-architecture)" in \
        amd64) suffix=x86_64 ;; \
        arm64) suffix=aarch64 ;; \
        *) echo "Архитектура не поддерживается: $(dpkg --print-architecture)" >&2; exit 1 ;; \
    esac; \
    curl -fsSL -o "${VINA_BIN}" \
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_${suffix}"; \
    chmod +x "${VINA_BIN}"; \
    "${VINA_BIN}" --version

# P2Rank: платформенно-независимый релиз (Java). Документацию и тестовые данные не берём.
RUN set -eux; \
    curl -fsSL "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/p2rank_${P2RANK_VERSION}.tar.gz" \
        | tar -xz -C /opt; \
    mv "/opt/p2rank_${P2RANK_VERSION}" /opt/p2rank; \
    rm -rf /opt/p2rank/test_data /opt/p2rank/doc; \
    chmod +x "${P2RANK_BIN}"

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY . /app
RUN chmod +x /app/docker/entrypoint.sh \
    && mkdir -p /app/outputs /app/data/cache

# Работает от root намеренно: каталоги outputs/ и data/cache/ монтируются с хоста, и пользователь
# с другим uid не смог бы в них писать. Контейнер ходит только к публичным научным API.
ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["demo", "--export-dir", "outputs/demo", "--no-docs"]
