"""
Tests for the v1.3.0 review module + analyze --apply-overrides flow.

Covers:
- ``make_review_overrides_json`` writes a schema-conformant file.
- The Flask app (test client; no real browser) renders the cluster
  table from a fixture annotations.csv.
- POST /save round-trips through the JSON shape.
- ``analyze --apply-overrides`` lands the right `source="user"` rows in
  the decision log.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scellrun.review import (
    OVERRIDES_FILENAME,
    OVERRIDES_SCHEMA_VERSION,
    create_app,
    make_review_overrides_json,
    overrides_path,
)


def _seed_run_dir(tmp_path: Path) -> Path:
    """Build a minimal run-dir with the artifacts the review UI reads."""
    run = tmp_path / "run-review-test"
    (run / "01_qc").mkdir(parents=True)
    (run / "04_annotate").mkdir(parents=True)
    (run / "06_views").mkdir(parents=True)

    # 04_annotate/annotations.csv — three clusters
    (run / "04_annotate" / "annotations.csv").write_text(
        "cluster,panel_label,panel_score,panel_margin,panel_rationale,"
        "ai_label,ai_rationale,top_markers\n"
        "0,Chondrocytes,0.42,0.18,COL2A1,,,COL2A1\n"
        "1,Macrophages,0.50,0.30,CD68,,,CD68\n"
        "2,InfC,0.20,0.05,IL6,,,IL6\n"
    )

    # 00_decisions.jsonl — with a few QC threshold rows so the slider
    # default-fill logic exercises both the "found in log" and the
    # "fall back to constant" paths.
    rows = [
        {"schema_version": 1, "stage": "qc", "key": "max_pct_mt", "value": 20,
         "default": 20, "source": "auto", "rationale": "default",
         "fix_payload": None, "attempt_id": "abc",
         "ts": "2026-04-30T09:45:30+00:00"},
        {"schema_version": 1, "stage": "qc", "key": "max_genes", "value": 4500,
         "default": 4000, "source": "user", "rationale": "user override",
         "fix_payload": None, "attempt_id": "abc",
         "ts": "2026-04-30T09:45:30+00:00"},
    ]
    (run / "00_decisions.jsonl").write_text(
        "\n".join(json.dumps(r) for r in rows) + "\n"
    )

    # 00_run.json — pin lang to en
    (run / "00_run.json").write_text(json.dumps({
        "created_at": "2026-04-30T09:45:30+00:00",
        "stages": [{"command": "analyze:qc", "ran_at": "...", "params": {"lang": "en"}}],
    }))

    return run


def test_make_review_overrides_json_writes_valid_schema(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    out = make_review_overrides_json(
        run,
        cluster_label_overrides={"5": "Pericyte"},
        cell_exclusions=["AAACCTG-1-1", "AAACCTG-2-1"],
        threshold_overrides={"max_pct_mt": 18, "max_genes": 4500},
        notes="cluster 5 looks like pericyte not chondrocyte",
        reviewer="dr_test",
    )
    assert out.exists()
    assert out.name == OVERRIDES_FILENAME
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["schema_version"] == OVERRIDES_SCHEMA_VERSION
    assert data["reviewer"] == "dr_test"
    assert data["cluster_label_overrides"] == {"5": "Pericyte"}
    assert data["cell_exclusions"] == ["AAACCTG-1-1", "AAACCTG-2-1"]
    assert data["threshold_overrides"]["max_pct_mt"] == 18
    assert "saved_at" in data and "run_dir" in data
    assert data["notes"].startswith("cluster 5")


def test_review_index_renders_cluster_table_from_annotations(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    app = create_app(run, read_only=False, lang="en")
    client = app.test_client()
    resp = client.get("/")
    assert resp.status_code == 200
    body = resp.get_data(as_text=True)
    # Cluster IDs and current labels from annotations.csv must appear.
    assert "Chondrocytes" in body
    assert "Macrophages" in body
    assert "InfC" in body
    # Threshold sliders pre-fill from the decision log (max_genes=4500).
    assert "value=\"4500\"" in body or "value='4500'" in body
    # Save button is editable (not read-only) so the disabled attribute
    # is absent.
    assert "<button id=\"save_btn\">" in body or "<button id=\"save_btn\" >" in body


def test_review_save_round_trips(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    app = create_app(run, read_only=False, lang="en")
    client = app.test_client()
    payload = {
        "cluster_label_overrides": {"5": "Pericyte", "9": "Osteoclast"},
        "cell_exclusions": "AAACCTG-1-1,AAACCTG-2-1\nAAACCTG-3-1",
        "threshold_overrides": {"max_pct_mt": 18, "max_genes": 4500, "min_counts": 600},
        "notes": "edited via test client",
    }
    resp = client.post("/save", json=payload)
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["ok"] is True
    saved = json.loads(overrides_path(run).read_text(encoding="utf-8"))
    assert saved["schema_version"] == 1
    assert saved["cluster_label_overrides"] == {"5": "Pericyte", "9": "Osteoclast"}
    assert saved["cell_exclusions"] == ["AAACCTG-1-1", "AAACCTG-2-1", "AAACCTG-3-1"]
    assert saved["threshold_overrides"]["max_pct_mt"] == 18
    assert "edited via test client" in saved["notes"]


def test_review_read_only_rejects_save(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    app = create_app(run, read_only=True, lang="en")
    client = app.test_client()
    resp = client.post("/save", json={"notes": "shouldn't save"})
    assert resp.status_code == 403
    assert not overrides_path(run).exists()
    body = client.get("/").get_data(as_text=True)
    # Read-only label surfaces on the page.
    assert "read-only" in body or "read-only" in body.lower()
    # Save button has disabled attribute.
    assert "disabled" in body


def test_review_zh_template_uses_zh_strings(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    app = create_app(run, read_only=False, lang="zh")
    body = app.test_client().get("/").get_data(as_text=True)
    assert "复审" in body
    assert "聚类标签" in body


def test_apply_overrides_writes_user_decisions(tmp_path: Path, monkeypatch) -> None:
    """
    `analyze --apply-overrides` must apply the overrides and record
    `source="user"` rows in the decision log. We test the orchestrator
    machinery directly via a fake annotated run-dir + a manual
    `_apply_post_annotate_overrides` call rather than the full pipeline
    (the full pipeline is exercised by test_analyze.py).
    """
    from scellrun.analyze import _apply_post_annotate_overrides, _load_overrides
    from scellrun.decisions import read_decisions

    run = _seed_run_dir(tmp_path)
    annot_out = run / "04_annotate"
    integrated = run / "02_integrate" / "integrated.h5ad"
    integrated.parent.mkdir(parents=True, exist_ok=True)
    integrated.touch()

    overrides_p = make_review_overrides_json(
        run,
        cluster_label_overrides={"0": "ProC", "2": "Pericyte"},
        cell_exclusions=["GHOST-CELL-A", "GHOST-CELL-B"],
        notes="manual review",
    )

    data = _load_overrides(overrides_p)
    assert data["cluster_label_overrides"] == {"0": "ProC", "2": "Pericyte"}

    _apply_post_annotate_overrides(
        run_dir=run,
        annot_out=annot_out,
        integrated_h5ad_path=integrated,
        cluster_label_overrides=data["cluster_label_overrides"],
        cell_exclusions=data["cell_exclusions"],
        attempt_id="testattempt",
    )

    # annotations.csv now has final_label / label_source columns.
    rows = list(csv.DictReader((annot_out / "annotations.csv").open()))
    by_cluster = {r["cluster"]: r for r in rows}
    assert by_cluster["0"]["final_label"] == "ProC"
    assert by_cluster["0"]["label_source"] == "user"
    assert by_cluster["1"]["final_label"] == "Macrophages"  # untouched → auto
    assert by_cluster["1"]["label_source"] == "auto"
    assert by_cluster["2"]["final_label"] == "Pericyte"
    assert by_cluster["2"]["label_source"] == "user"

    # Decision log carries the user-source rows.
    decisions = read_decisions(run)
    user_rows = [d for d in decisions if d.get("source") == "user"]
    label_keys = {
        d["key"] for d in user_rows
        if d["key"].startswith("label_override.cluster_")
    }
    assert "label_override.cluster_0" in label_keys
    assert "label_override.cluster_2" in label_keys
    excl_rows = [d for d in user_rows if d["key"] == "cell_exclusions"]
    assert excl_rows and excl_rows[0]["value"] == 2


def test_load_overrides_handles_missing_path() -> None:
    from scellrun.analyze import _load_overrides

    assert _load_overrides(None) == {}


def test_load_overrides_rejects_non_dict(tmp_path: Path) -> None:
    from scellrun.analyze import _load_overrides

    bad = tmp_path / "bad.json"
    bad.write_text("[1, 2, 3]")
    with pytest.raises(ValueError):
        _load_overrides(bad)
