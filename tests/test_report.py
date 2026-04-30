"""Tests for the multi-stage report aggregator."""
from __future__ import annotations

from pathlib import Path

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


def test_build_report_at_a_glance_block_renders(tmp_path):
    """
    v0.9.1 (#007): the index page now carries an "At a glance" section
    that surfaces QC pass-rate, cluster count, and top labels from the
    stage artifacts (not just stage links).
    """
    run = _make_run_dir(tmp_path)
    # Add a per_cell_metrics.csv with the QC pass column populated.
    (run / "01_qc" / "per_cell_metrics.csv").write_text(
        "cell,n_genes_by_counts,total_counts,pct_counts_mt,pct_counts_ribo,pct_counts_hb,scellrun_qc_pass\n"
        "c0,500,1000,5.0,10.0,1.0,True\n"
        "c1,300,800,8.0,15.0,2.0,True\n"
        "c2,100,200,40.0,5.0,1.0,False\n"
    )
    # Add an annotations.csv with three clusters, two labeled the same.
    (run / "04_annotate").mkdir()
    (run / "04_annotate" / "annotations.csv").write_text(
        "cluster,panel_label,panel_score,panel_margin,ai_label,ai_rationale,top_markers\n"
        "0,Macrophages,0.4,0.2,,,CD68\n"
        "1,Macrophages,0.3,0.1,,,CD163\n"
        "2,T cells,0.5,0.3,,,CD3D\n"
    )
    out = run / "05_report"
    artifacts = build_report(run, out)
    html = artifacts["index"].read_text(encoding="utf-8")
    assert "At a glance" in html
    assert "QC pass rate" in html
    # 2/3 = 66.7%
    assert "66.7" in html
    # cluster total of 3
    assert "Macrophages" in html or "T cells" in html


def test_build_report_at_a_glance_zh(tmp_path):
    """ZH template carries the same glance block in Chinese."""
    run = _make_run_dir(tmp_path)
    out = run / "05_report"
    artifacts = build_report(run, out, lang="zh")
    html = artifacts["index"].read_text(encoding="utf-8")
    assert "一览" in html
