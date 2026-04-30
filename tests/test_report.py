"""Tests for the multi-stage report aggregator."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scellrun.report import _scan_stage, build_report
from scellrun.runlayout import write_run_meta


def _make_run_dir(tmp_path: Path) -> Path:
    """Build a minimal run-dir layout with stub stage outputs."""
    run = tmp_path / "run-test"
    (run / "01_qc").mkdir(parents=True)
    (run / "01_qc" / "report.html").write_text("<html>qc report</html>")
    (run / "01_qc" / "per_cell_metrics.csv").write_text("a,b\n1,2\n")

    (run / "02_integrate").mkdir()
    (run / "02_integrate" / "report.html").write_text("<html>integrate report</html>")
    (run / "02_integrate" / "umap_grid.png").write_bytes(b"\x89PNG fake")
    return run


def test_scan_stage_present(tmp_path):
    run = _make_run_dir(tmp_path)
    info = _scan_stage(run, "qc")
    assert info.present is True
    assert info.report_path is not None
    assert "report.html" in info.artifacts
    assert "per_cell_metrics.csv" in info.artifacts


def test_scan_stage_absent(tmp_path):
    run = _make_run_dir(tmp_path)
    info = _scan_stage(run, "annotate")
    assert info.present is False
    assert info.report_path is None
    assert info.artifacts == {}


def test_build_report_emits_index(tmp_path):
    run = _make_run_dir(tmp_path)
    write_run_meta(run, command="scrna qc", params={"foo": 1})
    write_run_meta(run, command="scrna integrate", params={"bar": 2})

    out = run / "05_report"
    artifacts = build_report(run, out)

    assert artifacts["index"].exists()
    html = artifacts["index"].read_text()
    # mentions present stages
    assert "qc" in html
    assert "integrate" in html
    # mentions absent ones too (greyed out)
    assert "markers" in html
    assert "annotate" in html
    # manifest entries surface
    assert "scrna qc" in html
    assert "scrna integrate" in html


def test_build_report_zh(tmp_path):
    run = _make_run_dir(tmp_path)
    write_run_meta(run, command="scrna qc", params={"foo": 1})
    out = run / "05_report"
    artifacts = build_report(run, out, lang="zh")
    html = artifacts["index"].read_text(encoding="utf-8")
    assert "完整流程报告" in html


def test_build_report_no_manifest_ok(tmp_path):
    """Aggregator should still work if 00_run.json doesn't exist."""
    run = _make_run_dir(tmp_path)
    out = run / "05_report"
    artifacts = build_report(run, out)
    assert artifacts["index"].exists()
