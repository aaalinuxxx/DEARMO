"""Один устойчивый HTTP-клиент на все научные API.

Раньше `_get_json_with_retry` был скопирован в три клиента с разными таймаутами
и разными списками перехваченных исключений, поэтому «источник недоступен»
означало в каждом файле своё. Здесь это одно поведение: ретраи только на
временных кодах, остальные ошибки возвращают None сразу.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from typing import Any

# Коды, при которых повтор имеет смысл: сервер перегружен или временно упал.
_RETRY_CODES = frozenset({429, 500, 502, 503, 504})

_DEFAULT_HEADERS = {"Accept": "application/json"}


class Transport:
    """Исход последнего запроса — чтобы агент отличал «нет данных» от «нет связи».

    Это разные ответы пользователю: в первом случае измерений действительно нет,
    во втором мы просто ничего не знаем и обязаны это сказать.
    """

    def __init__(self) -> None:
        self.reachable: bool | None = None
        self.detail: str | None = None

    def reset(self) -> None:
        self.reachable = None
        self.detail = None

    def ok(self) -> None:
        self.reachable = True
        self.detail = None

    def fail(self, detail: str) -> None:
        self.reachable = False
        self.detail = detail

    def as_dict(self) -> dict[str, Any]:
        return {"reachable": self.reachable, "detail": self.detail}


def request_json(
    url: str,
    *,
    body: Any = None,
    timeout: float = 10.0,
    max_retries: int = 3,
    headers: dict[str, str] | None = None,
    transport: Transport | None = None,
) -> Any | None:
    """GET, либо POST с JSON-телом, если передан `body`. None — запрос не удался."""
    merged = dict(_DEFAULT_HEADERS)
    if body is not None:
        merged["Content-Type"] = "application/json"
    merged.update(headers or {})
    payload = json.dumps(body).encode("utf-8") if body is not None else None

    last_detail = "неизвестная ошибка"
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(
                url,
                data=payload,
                headers=merged,
                method="POST" if payload is not None else "GET",
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            if transport is not None:
                transport.ok()
            return data
        except urllib.error.HTTPError as exc:
            last_detail = f"HTTP {exc.code}"
            if exc.code in _RETRY_CODES and attempt + 1 < max_retries:
                time.sleep(0.2 * (attempt + 1))
                continue
            break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_detail = f"сеть недоступна: {exc}"
            if attempt + 1 < max_retries:
                time.sleep(0.1 * (attempt + 1))
                continue
            break
        except (json.JSONDecodeError, ValueError) as exc:
            # Ответ пришёл, но это не JSON: повтор не поможет.
            last_detail = f"ответ не является JSON: {exc}"
            break

    if transport is not None:
        transport.fail(last_detail)
    return None
