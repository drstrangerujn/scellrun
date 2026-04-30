"""
Schema-only tests for the v1.0 agent-demo deliverables under docs/.

These check that the dogfooded run on real OA cartilage data produced the
artifacts the README + agent-demo.md reference, and that the demo dialogue
is grounded in the actual run-dir (not template prose).

The test suite is deliberately structural — it does not re-run the
analyze pipeline (that is a 5-7 minute hospital-server-only process).
The artifacts are committed under docs/v1demo/.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEMO_DIR = REPO_ROOT / "docs" / "v1demo"
INDEX_HTML = DEMO_DIR / "index.html"
DECISIONS_JSONL = DEMO_DIR / "decisions.jsonl"
RUN_JSON = DEMO_DIR / "run.json"
AGENT_DEMO_MD = REPO_ROOT / "docs" / "agent-demo.md"


def test_index_html_present_and_nonempty() -> None:
    """The dogfooded index.html must exist and have meaningful content."""
    assert INDEX_HTML.exists(), f"missing {INDEX_HTML}"
    body = INDEX_HTML.read_text(encoding="utf-8")
    assert len(body) > 1000, f"index.html suspiciously small ({len(body)} bytes)"
    # The "At a glance" block landed in v0.9.1 and is required for v1.0
    assert "At a glance" in body, "index.html should embed the At-a-glance block"


def test_decisions_jsonl_has_rows_with_required_fields() -> None:
    """decisions.jsonl must have >=10 rows; every row must be valid JSON
    with schema_version and attempt_id present (v0.9.1 schema)."""
    assert DECISIONS_JSONL.exists(), f"missing {DECISIONS_JSONL}"
    rows = [
        json.loads(line)
        for line in DECISIONS_JSONL.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(rows) >= 10, f"decisions.jsonl has only {len(rows)} rows; expected >=10"
    for i, row in enumerate(rows):
        assert "schema_version" in row, f"row {i} missing schema_version"
        assert "attempt_id" in row, f"row {i} missing attempt_id"


def test_run_json_present_and_parseable() -> None:
    """run.json must exist and parse as JSON with a stages list."""
    assert RUN_JSON.exists(), f"missing {RUN_JSON}"
    data = json.loads(RUN_JSON.read_text(encoding="utf-8"))
    assert "stages" in data, "run.json must have a 'stages' key"
    assert isinstance(data["stages"], list), "run.json stages must be a list"
    assert len(data["stages"]) >= 1, "run.json should have at least one stage entry"


def test_agent_demo_md_exists_and_has_turn_markers() -> None:
    """agent-demo.md must contain the four conversational turn markers."""
    assert AGENT_DEMO_MD.exists(), f"missing {AGENT_DEMO_MD}"
    text = AGENT_DEMO_MD.read_text(encoding="utf-8")
    user_turns = re.findall(r"^> \*\*User:\*\*", text, flags=re.MULTILINE)
    agent_turns = re.findall(r"^> \*\*Agent:\*\*", text, flags=re.MULTILINE)
    assert len(user_turns) >= 4, f"agent-demo.md has {len(user_turns)} User turns; expected >=4"
    assert len(agent_turns) >= 4, f"agent-demo.md has {len(agent_turns)} Agent turns; expected >=4"


def test_agent_demo_md_references_actual_run_dir() -> None:
    """agent-demo.md must cite the run-dir from run.json, proving the
    dialogue is grounded in a real run, not a template."""
    text = AGENT_DEMO_MD.read_text(encoding="utf-8")
    data = json.loads(RUN_JSON.read_text(encoding="utf-8"))

    run_dir_path = data.get("run_dir")
    if not run_dir_path:
        # fallback: scrape the first stage's out_dir for a run-* segment
        stages = data.get("stages", [])
        for s in stages:
            params = s.get("params") or {}
            out = params.get("out_dir") or ""
            m = re.search(r"run-\d{8}-\d{6}", str(out))
            if m:
                run_dir_path = m.group(0)
                break

    assert run_dir_path, "could not extract run-dir id from run.json"
    run_id = Path(run_dir_path).name if "/" in str(run_dir_path) else str(run_dir_path)
    assert run_id in text, (
        f"agent-demo.md should reference the actual run-dir id {run_id!r}; "
        "the demo must be specific to the dogfooded run, not generic prose"
    )
