import httpx

from orchestration.llm_client import plan_with_ollama


class _Response:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise httpx.HTTPStatusError("fail", request=None, response=None)


def test_plan_with_ollama_parses_json(monkeypatch):
    monkeypatch.setenv("DRUGDIS_LLM_PROVIDER", "ollama")
    monkeypatch.setenv("OLLAMA_HOST", "http://ollama.test")

    def fake_get(url, timeout=0):
        assert url.endswith("/api/tags")
        return _Response({"models": []})

    def fake_post(url, json=None, timeout=0):
        assert url.endswith("/v1/chat/completions")
        return _Response(
            {
                "choices": [
                    {
                        "message": {
                            "content": '{"scenario_id":"target_for_ligand","ligands":["acetazolamide"],"pdb_id":"3HS4"}'
                        }
                    }
                ]
            }
        )

    monkeypatch.setattr(httpx, "get", fake_get)
    monkeypatch.setattr(httpx, "post", fake_post)
    extracted = plan_with_ollama("подбери белок", "system")
    assert extracted is not None
    assert extracted.scenario_id == "target_for_ligand"
    assert extracted.pdb_id == "3HS4"
    assert extracted.source.startswith("ollama:")
