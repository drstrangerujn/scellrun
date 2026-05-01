"""
Tests for the v1.2.0 views layer (06_views/).

Build a fake run-dir with the upstream stage artifacts the views layer
references, call build_views, then assert each cross-section page
exists and pulls in the expected data.

The views layer is intentionally tolerant of missing inputs (best-effort
aggregator), so the bulk of these tests build a happy-path layout and
spot-check missing-input branches separately.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scellrun.views import build_views


def _make_run_dir(tmp_path: Path) -> Path:
    run = tmp_path / "run-views-test"

    # 01_qc — minimal so build_views' callers don't trip up
    (run / "01_qc").mkdir(parents=True)
    (run / "01_qc" / "report.html").write_text("<html>qc</html>")

    # 02_integrate — UMAP grid + cluster sizes + cluster_by_sample.csv
    (run / "02_integrate").mkdir()
    (run / "02_integrate" / "umap_grid.png").write_bytes(b"\x89PNG fake")
    (run / "02_integrate" / "cluster_sizes.png").write_bytes(b"\x89PNG fake")
    (run / "02_integrate" / "cluster_by_sample.csv").write_text(
        "leiden_res_0_5,sampleA,sampleB\n"
        "0,40,30\n"
        "1,20,18\n"
        "2,10,3\n"
    )

    # 03_markers — one CSV at res 0.5
    (run / "03_markers").mkdir()
    (run / "03_markers" / "markers_res_0.5.csv").write_text(
        "cluster,gene,log2fc,pct_in,pct_out\n"
        "0,COL2A1,3.2,0.95,0.05\n"
        "0,ACAN,2.8,0.90,0.10\n"
        "1,CD68,4.1,0.88,0.07\n"
        "1,CD163,3.5,0.80,0.10\n"
        "2,IL6,2.1,0.70,0.20\n"
    )

    # 04_annotate — annotations.csv at the chosen resolution
    (run / "04_annotate").mkdir()
    (run / "04_annotate" / "annotations.csv").write_text(
        "cluster,panel_label,panel_score,panel_margin,panel_rationale,ai_label,ai_rationale,top_markers\n"
        "0,Chondrocytes,0.42,0.18,\"Chondrocytes (0.42): COL2A1, ACAN\",,,COL2A1\n"
        "1,Macrophages,0.50,0.30,\"Macrophages (0.50): CD68, CD163\",,,CD68\n"
        "2,InfC,0.20,0.05,\"InfC (0.20): IL6\",,,IL6\n"
    )

    # 00_decisions.jsonl — a few rows including the chosen-resolution row
    decisions = [
        {
            "schema_version": 1,
            "stage": "qc",
            "key": "profile",
            "value": "joint-disease",
            "default": "default",
            "source": "user",
            "rationale": "user passed --profile joint-disease",
            "fix_payload": None,
            "attempt_id": "abc",
            "ts": "2026-04-30T09:45:30+00:00",
        },
        {
            "schema_version": 1,
            "stage": "analyze",
            "key": "chosen_resolution_for_annotate",
            "value": 0.5,
            "default": None,
            "source": "auto",
            "rationale": "largest n_clusters among non-fragmented",
            "fix_payload": None,
            "attempt_id": "abc",
            "ts": "2026-04-30T09:48:01+00:00",
        },
        {
            "schema_version": 1,
            "stage": "integrate",
            "key": "resolution_recommended",
            "value": 0.5,
            "default": None,
            "source": "ai",
            "rationale": "LLM picked 0.5 from the sweep",
            "fix_payload": None,
            "attempt_id": "abc",
            "ts": "2026-04-30T09:48:01+00:00",
        },
        {
            "schema_version": 1,
            "stage": "qc",
            "key": "min_genes",
            "value": 200,
            "default": 200,
            "source": "auto",
            "rationale": "AIO NRMI=200",
            "fix_payload": None,
            "attempt_id": "abc",
            "ts": "2026-04-30T09:45:30+00:00",
        },
    ]
    (run / "00_decisions.jsonl").write_text(
        "\n".join(json.dumps(d) for d in decisions) + "\n"
    )

    return run


def test_build_views_emits_expected_files(tmp_path: Path) -> None:
    run = _make_run_dir(tmp_path)
    out = build_views(run)

    # Each top-level view must exist on disk.
    assert "index" in out and out["index"].exists()
    assert "by_decision_source" in out and out["by_decision_source"].exists()


def test_by_resolution_index_contains_marker_genes(tmp_path: Path) -> None:
    run = _make_run_dir(tmp_path)
    build_views(run)

    page = run / "06_views" / "by_resolution" / "0.5" / "index.html"
    assert page.exists()
    html = page.read_text(encoding="utf-8")

    # Marker genes from markers_res_0.5.csv must appear inline.
    assert "COL2A1" in html
    assert "CD68" in html
    # The footer link to the full markers CSV is present (relative path).
    assert "../../../03_markers/markers_res_0.5.csv" in html


def test_by_cluster_page_contains_label_from_annotations(tmp_path: Path) -> None:
    run = _make_run_dir(tmp_path)
    build_views(run)

    page = run / "06_views" / "by_cluster" / "cluster_0.html"
    assert page.exists()
    html = page.read_text(encoding="utf-8")

    # Label from annotations.csv ("Chondrocytes" for cluster 0).
    assert "Chondrocytes" in html
    # Top markers from the markers CSV for cluster 0.
    assert "COL2A1" in html
    assert "ACAN" in html
    # Sample-distribution table from cluster_by_sample.csv.
    assert "sampleA" in html
    assert "sampleB" in html

    # Per-cluster cell barcode CSV exists next to the page (footer link).
    cells_csv = run / "06_views" / "by_cluster" / "cluster_0_cells.csv"
    assert cells_csv.exists()


def test_by_decision_source_groups_three_sources(tmp_path: Path) -> None:
    run = _make_run_dir(tmp_path)
    build_views(run)

    page = run / "06_views" / "by_decision_source" / "index.html"
    assert page.exists()
    html = page.read_text(encoding="utf-8")

    # The three sections must be rendered.
    assert "AI-driven" in html
    assert "User overrides" in html
    assert "Auto picks" in html

    # An ai-row's rationale must surface.
    assert "LLM picked 0.5 from the sweep" in html
    # A user-row.
    assert "user passed --profile joint-disease" in html
    # The boring auto row (min_genes default) should be filtered out.
    assert "AIO NRMI=200" not in html


def test_views_index_links_resolutions_and_clusters(tmp_path: Path) -> None:
    run = _make_run_dir(tmp_path)
    build_views(run)

    landing = run / "06_views" / "index.html"
    assert landing.exists()
    html = landing.read_text(encoding="utf-8")

    # Resolution links.
    assert "by_resolution/0.5/index.html" in html
    # Cluster links.
    assert "by_cluster/cluster_0.html" in html
    assert "by_cluster/cluster_1.html" in html
    # Back-link to the main report.
    assert "../05_report/index.html" in html


def test_no_symlinks_or_copies_used(tmp_path: Path) -> None:
    """
    Hard rule (v1.2.0): the views layer is pure HTML + relative paths.
    Verify build_views doesn't create any symlinks under 06_views/.
    """
    run = _make_run_dir(tmp_path)
    build_views(run)

    views_root = run / "06_views"
    for p in views_root.rglob("*"):
        assert not p.is_symlink(), f"unexpected symlink at {p}"


def test_build_views_is_idempotent(tmp_path: Path) -> None:
    """Re-running build_views over the same run-dir should be a no-op overwrite."""
    run = _make_run_dir(tmp_path)
    build_views(run)
    # Mutate the resolution page; second call must restore it.
    page = run / "06_views" / "by_resolution" / "0.5" / "index.html"
    page.write_text("STALE", encoding="utf-8")
    build_views(run)
    assert "STALE" not in page.read_text(encoding="utf-8")


def test_build_views_handles_missing_artifacts(tmp_path: Path) -> None:
    """
    With an empty run-dir, build_views should still emit a landing page
    and the by_decision_source page; the other views degrade to
    placeholders rather than crash.
    """
    run = tmp_path / "empty-run"
    run.mkdir()
    out = build_views(run)
    assert out["index"].exists()
    assert out["by_decision_source"].exists()


@pytest.mark.parametrize("name", ["index.html", "by_decision_source/index.html"])
def test_emitted_pages_use_relative_paths_only(tmp_path: Path, name: str) -> None:
    """
    The view pages must reference upstream artifacts via relative paths
    (../*) — not by absolute paths or file:// URLs that would tie them
    to one machine.
    """
    run = _make_run_dir(tmp_path)
    build_views(run)
    html = (run / "06_views" / name).read_text(encoding="utf-8")
    assert "file://" not in html
    # No absolute filesystem paths leaked in (heuristic: no /tmp/ anywhere
    # the test run-dir actually lives on tmp_path).
    assert str(tmp_path) not in html
