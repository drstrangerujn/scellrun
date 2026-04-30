"""Tests for the v0.7 decision log."""
from __future__ import annotations

import json
import threading

import anndata as ad
import numpy as np
import pytest

from scellrun.decisions import (
    DECISIONS_FILENAME,
    SCHEMA_VERSION,
    Decision,
    decisions_path,
    group_by_stage,
    read_decisions,
    record,
    record_many,
    truncate_stage,
)


def test_decision_roundtrip(tmp_path):
    """Write one Decision, read it back as a dict with the right fields."""
    d = Decision(
        stage="qc",
        key="max_pct_mt",
        value=20.0,
        default=10.0,
        source="auto",
        rationale="joint tissue is stress-prone",
        attempt_id="abc123",
    )
    record(tmp_path, d)
    p = decisions_path(tmp_path)
    assert p.exists()
    assert p.name == DECISIONS_FILENAME

    rows = read_decisions(tmp_path)
    assert len(rows) == 1
    assert rows[0]["stage"] == "qc"
    assert rows[0]["key"] == "max_pct_mt"
    assert rows[0]["value"] == 20.0
    assert rows[0]["default"] == 10.0
    assert rows[0]["source"] == "auto"
    assert "stress-prone" in rows[0]["rationale"]
    assert "ts" in rows[0]
    # v0.9.1: schema_version + attempt_id + fix_payload round-trip.
    assert rows[0]["schema_version"] == SCHEMA_VERSION
    assert rows[0]["attempt_id"] == "abc123"
    assert rows[0]["fix_payload"] is None


def test_decision_from_choice_user_override():
    """from_choice with is_user_override=True yields source='user'."""
    d = Decision.from_choice(
        stage="qc",
        key="max_pct_mt",
        value=10.0,
        default=20.0,
        is_user_override=True,
        rationale="user tightened it",
    )
    assert d.source == "user"


def test_decision_from_choice_auto_when_value_differs_from_default():
    """
    A2: value != default does NOT imply user override. The orchestrator's
    chosen_resolution_for_annotate is `auto` regardless of whether the
    chosen value happens to differ from a hardcoded default.
    """
    d = Decision.from_choice(
        stage="analyze",
        key="chosen_resolution_for_annotate",
        value=0.3,
        default=None,
        is_user_override=False,
        rationale="picked from quality table",
    )
    assert d.source == "auto"


def test_decision_fix_payload_round_trip(tmp_path):
    """fix_payload survives the JSONL round-trip and stays a dict."""
    d = Decision(
        stage="qc",
        key="self_check.qc_low_pass_rate.suggest",
        value="raise --max-pct-mt to 25",
        source="auto",
        rationale="...",
        fix_payload={"max_pct_mt": 25},
    )
    record(tmp_path, d)
    rows = read_decisions(tmp_path)
    assert rows[0]["fix_payload"] == {"max_pct_mt": 25}


def test_decision_appends(tmp_path):
    """Multiple decisions append; one JSONL line each, in order."""
    record(tmp_path, Decision(stage="qc", key="a", value=1, source="auto"))
    record(tmp_path, Decision(stage="qc", key="b", value=2, source="user"))
    record(tmp_path, Decision(stage="integrate", key="c", value=3, source="ai"))

    p = decisions_path(tmp_path)
    lines = p.read_text().strip().splitlines()
    assert len(lines) == 3
    for line in lines:
        json.loads(line)  # each line must parse on its own

    rows = read_decisions(tmp_path)
    assert [r["key"] for r in rows] == ["a", "b", "c"]


def test_record_many(tmp_path):
    """record_many writes a list as one-line-each."""
    record_many(
        tmp_path,
        [
            Decision(stage="qc", key="x", value="x_val", source="auto"),
            Decision(stage="qc", key="y", value="y_val", source="auto"),
        ],
    )
    rows = read_decisions(tmp_path)
    assert len(rows) == 2


def test_record_none_run_dir_is_noop(tmp_path):
    """Passing run_dir=None must not blow up; no file is written."""
    record(None, Decision(stage="qc", key="k", value="v", source="auto"))
    record_many(None, [Decision(stage="qc", key="k", value="v", source="auto")])
    assert not (tmp_path / DECISIONS_FILENAME).exists()


def test_invalid_source_raises():
    with pytest.raises(ValueError):
        Decision(stage="qc", key="k", value=1, source="bogus")


def test_read_decisions_missing_returns_empty(tmp_path):
    assert read_decisions(tmp_path) == []


def test_group_by_stage_preserves_order(tmp_path):
    decisions = [
        {"stage": "qc", "key": "a"},
        {"stage": "integrate", "key": "b"},
        {"stage": "qc", "key": "c"},
        {"stage": "annotate", "key": "d"},
    ]
    grouped = group_by_stage(decisions)
    assert list(grouped.keys()) == ["qc", "integrate", "annotate"]
    assert [d["key"] for d in grouped["qc"]] == ["a", "c"]


def test_jsonable_handles_paths_dataclasses_and_collections(tmp_path):
    """Path / dataclass / set / tuple values should serialise without error."""
    from dataclasses import dataclass

    @dataclass
    class Inner:
        a: int
        b: str

    d = Decision(
        stage="qc",
        key="mixed",
        value={"path": tmp_path, "inner": Inner(1, "x"), "tup": (1, 2), "set": {3, 4}},
        source="auto",
    )
    record(tmp_path, d)
    rows = read_decisions(tmp_path)
    assert rows[0]["value"]["inner"]["a"] == 1
    assert rows[0]["value"]["tup"] == [1, 2]


# ---------- end-to-end: pipeline writes a non-empty decisions log -----------


@pytest.fixture
def planted_h5ad(tmp_path):
    """Same as test_analyze.planted_h5ad, copied here so the tests are independent."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 500
    gene_means = rng.gamma(shape=2.0, scale=1.5, size=n_genes).astype(np.float32) + 0.3
    cell_scaling = rng.gamma(shape=4.0, scale=0.25, size=n_cells).astype(np.float32) + 0.5
    lam = np.outer(cell_scaling, gene_means)
    counts = rng.poisson(lam=lam).astype(np.int32)
    counts[:100, 50:70] = rng.poisson(lam=25.0, size=(100, 20))
    counts[100:, 70:90] = rng.poisson(lam=25.0, size=(100, 20))

    a = ad.AnnData(X=counts.astype(np.float32))
    a.var_names = (
        [f"MT-CO{i}" for i in range(1, 6)]
        + [f"HBB{i}" for i in range(1, 4)]
        + [f"RPS{i}" for i in range(1, 11)]
        + [f"GENE{i}" for i in range(n_genes - 18)]
    )
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["sample"] = ["A"] * 100 + ["B"] * 100
    out = tmp_path / "planted.h5ad"
    a.write_h5ad(out)
    return out


def test_analyze_pipeline_writes_non_empty_decisions(planted_h5ad, tmp_path, monkeypatch):
    """End-to-end: a clean analyze run produces a non-empty 00_decisions.jsonl
    with decisions from every stage including the orchestrator's own pick."""
    from scellrun.analyze import run_analyze

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"

    run_analyze(
        planted_h5ad,
        profile="joint-disease",
        species="human",
        tissue="synthetic",
        resolutions=(0.3, 0.5),
        use_ai=False,
        lang="en",
        run_dir=run_dir,
        force=False,
        method="none",
        regress_cell_cycle=False,
        use_pubmed=False,
    )

    log = decisions_path(run_dir)
    assert log.exists()
    rows = read_decisions(run_dir)
    # Spec calls for ~10-15 decisions in a no-override run; sanity-check it's
    # at least 8 to avoid pinning to an exact count that future tweaks would break.
    assert len(rows) >= 8

    stages_seen = {r["stage"] for r in rows}
    # Every stage logged at least one decision.
    for s in ("qc", "integrate", "markers", "annotate", "analyze"):
        assert s in stages_seen, f"no decisions logged for stage {s!r}"

    # The orchestrator's chosen-resolution decision is present and points at one
    # of the requested resolutions.
    chosen = [
        r for r in rows
        if r["stage"] == "analyze" and r["key"] == "chosen_resolution_for_annotate"
    ]
    assert len(chosen) == 1
    assert float(chosen[0]["value"]) in (0.3, 0.5)


def test_truncate_stage_drops_only_matching_rows(tmp_path):
    """truncate_stage(run_dir, 'qc') removes every qc.* row, keeps others."""
    record_many(
        tmp_path,
        [
            Decision(stage="qc", key="profile", value="default", source="auto"),
            Decision(stage="qc", key="max_pct_mt", value=20.0, source="auto"),
            Decision(stage="integrate", key="method", value="harmony", source="auto"),
            Decision(stage="analyze", key="chosen", value=0.5, source="auto"),
        ],
    )
    n = truncate_stage(tmp_path, "qc")
    assert n == 2
    rows = read_decisions(tmp_path)
    assert {(r["stage"], r["key"]) for r in rows} == {
        ("integrate", "method"),
        ("analyze", "chosen"),
    }


def test_truncate_stage_no_op_when_file_missing(tmp_path):
    assert truncate_stage(tmp_path / "nope", "qc") == 0
    assert truncate_stage(tmp_path, "qc") == 0  # file doesn't exist yet


def test_concurrent_appends_no_truncated_lines(tmp_path):
    """
    A4: two threads writing record_many concurrently produce a JSONL file
    where line_count == thread_a_writes + thread_b_writes and every line
    parses as JSON. Without flock, interleaving from competing writers can
    leave half-flushed lines.
    """
    barrier = threading.Barrier(2)
    n_per_thread = 50
    decisions_a = [
        Decision(
            stage="qc",
            key=f"a_{i}",
            value=i,
            source="auto",
            rationale="x" * 200,  # big-ish payload to make interleaving plausible
        )
        for i in range(n_per_thread)
    ]
    decisions_b = [
        Decision(
            stage="integrate",
            key=f"b_{i}",
            value=i,
            source="auto",
            rationale="y" * 200,
        )
        for i in range(n_per_thread)
    ]

    def writer_a() -> None:
        barrier.wait()
        record_many(tmp_path, decisions_a)

    def writer_b() -> None:
        barrier.wait()
        record_many(tmp_path, decisions_b)

    t_a = threading.Thread(target=writer_a)
    t_b = threading.Thread(target=writer_b)
    t_a.start()
    t_b.start()
    t_a.join()
    t_b.join()

    p = decisions_path(tmp_path)
    raw_lines = p.read_text(encoding="utf-8").splitlines()
    assert len(raw_lines) == n_per_thread * 2
    for line in raw_lines:
        json.loads(line)  # must parse — no half-flushed half-line
    rows = read_decisions(tmp_path)
    keys = {r["key"] for r in rows}
    assert len(keys) == n_per_thread * 2


def test_multistage_report_renders_decision_summary(planted_h5ad, tmp_path, monkeypatch):
    """The aggregated index.html mentions 'Decision summary' / '决策摘要'."""
    from scellrun.analyze import run_analyze
    from scellrun.report import build_report

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"
    run_analyze(
        planted_h5ad,
        profile="joint-disease",
        run_dir=run_dir,
        resolutions=(0.5,),
        method="none",
        use_ai=False,
    )

    # English report (run_analyze already produced 05_report/index.html)
    en_html = (run_dir / "05_report" / "index.html").read_text(encoding="utf-8")
    assert "Decision summary" in en_html
    # Some specific decision keys should surface
    assert "max_pct_mt" in en_html
    assert "sample_key" in en_html

    # ZH render goes through build_report explicitly
    zh_out = tmp_path / "zh_report"
    artifacts = build_report(run_dir, zh_out, lang="zh")
    zh_html = artifacts["index"].read_text(encoding="utf-8")
    assert "决策摘要" in zh_html
