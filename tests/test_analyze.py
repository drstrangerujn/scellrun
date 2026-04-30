"""Tests for the v0.6 one-shot analyze orchestrator."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from scellrun.analyze import _pick_best_resolution, run_analyze


@pytest.fixture
def planted_h5ad(tmp_path):
    """
    Synthetic h5ad with planted two-population structure that survives QC and
    yields >= 2 Leiden clusters at res=0.5 after integrate.

    Two cell groups (rows 0-99 vs 100-199) with non-overlapping high-expression
    blocks; profile-default thresholds let all 200 cells through QC. Includes
    MT-/RPS-/HBB- genes so QC sub-tests have signal, and a 'sample' obs column
    so harmony / sample-key auto-detection has something to look at.
    """
    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 500
    # Raw integer counts so the QC raw-counts heuristic is happy.
    counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.int32)
    # Plant a strong group-A signal in cols 50..70, group-B signal in cols 70..90.
    counts[:100, 50:70] = rng.poisson(lam=20.0, size=(100, 20))
    counts[100:, 70:90] = rng.poisson(lam=20.0, size=(100, 20))
    # A handful of high-expressed "background" genes in the first 10 cols
    counts[:, :10] = rng.poisson(lam=15.0, size=(n_cells, 10))

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


def test_pick_best_resolution_prefers_more_clusters_when_clean():
    quality = {
        0.3: {"n_clusters": 3, "n_singletons": 0, "largest_pct": 50.0, "smallest_pct": 10.0,
              "mixing_entropy": 0.9},
        0.5: {"n_clusters": 5, "n_singletons": 0, "largest_pct": 40.0, "smallest_pct": 8.0,
              "mixing_entropy": 0.9},
        0.8: {"n_clusters": 7, "n_singletons": 3, "largest_pct": 35.0, "smallest_pct": 1.0,
              "mixing_entropy": 0.85},
    }
    assert _pick_best_resolution(quality) == 0.5


def test_pick_best_resolution_falls_back_when_empty():
    assert _pick_best_resolution({}) == 0.5
    assert _pick_best_resolution(None) == 0.5
    assert _pick_best_resolution({}, fallback=0.3) == 0.3


def test_pick_best_resolution_handles_only_fragmented():
    """If every resolution has many singletons, still pick the max n_clusters one."""
    quality = {
        0.3: {"n_clusters": 4, "n_singletons": 3, "largest_pct": 60.0, "smallest_pct": 0.5,
              "mixing_entropy": 0.5},
        0.5: {"n_clusters": 6, "n_singletons": 4, "largest_pct": 55.0, "smallest_pct": 0.4,
              "mixing_entropy": 0.5},
    }
    assert _pick_best_resolution(quality) == 0.5


def test_analyze_end_to_end_creates_all_stage_dirs(planted_h5ad, tmp_path, monkeypatch):
    """
    Synthetic h5ad → full pipeline runs end to end → all five stage subdirs
    exist → 05_report/index.html mentions every stage.

    --ai is OFF (CI has no API key).
    """
    # Prevent harmony import / batch correction failure from polluting the
    # test by using --method none. The orchestrator still exercises every
    # other code path.
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"

    result = run_analyze(
        planted_h5ad,
        profile="default",
        species="human",
        tissue="synthetic",
        resolutions=(0.3, 0.5),
        use_ai=False,
        lang="en",
        run_dir=run_dir,
        force=False,
        method="none",  # single-batch synthetic, no harmony
        regress_cell_cycle=False,
        use_pubmed=False,
    )

    # All five stage subdirs exist
    for sub in ("01_qc", "02_integrate", "03_markers", "04_annotate", "05_report"):
        assert (run_dir / sub).is_dir(), f"missing stage dir: {sub}"

    # Each stage produced its report.html (or for stage 5, index.html)
    assert (run_dir / "01_qc" / "report.html").exists()
    assert (run_dir / "01_qc" / "qc.h5ad").exists()
    assert (run_dir / "02_integrate" / "report.html").exists()
    assert (run_dir / "02_integrate" / "integrated.h5ad").exists()
    assert (run_dir / "03_markers" / "report.html").exists()
    assert (run_dir / "04_annotate" / "report.html").exists()
    assert (run_dir / "04_annotate" / "annotations.csv").exists()
    assert (run_dir / "05_report" / "index.html").exists()

    # 00_run.json captures every orchestrated stage
    manifest = (run_dir / "00_run.json").read_text()
    for cmd in (
        "analyze:qc",
        "analyze:integrate",
        "analyze:markers",
        "analyze:annotate",
        "analyze:report",
    ):
        assert cmd in manifest, f"manifest missing {cmd!r}"

    # Top-level index.html mentions every stage by name
    index_html = (run_dir / "05_report" / "index.html").read_text()
    for stage in ("qc", "integrate", "markers", "annotate"):
        assert stage in index_html

    # Result object exposes the chosen resolution + stage list
    assert set(result.stages_completed) == {"qc", "integrate", "markers", "annotate", "report"}
    assert result.report_index == run_dir / "05_report" / "index.html"
    assert result.chosen_resolution in (0.3, 0.5)


def test_analyze_force_allows_rerun(planted_h5ad, tmp_path, monkeypatch):
    """A repeat run with --force overwrites existing stage subdirs."""
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"
    run_analyze(
        planted_h5ad,
        run_dir=run_dir,
        resolutions=(0.5,),
        method="none",
        use_ai=False,
    )
    # second run on same run-dir without force should fail
    from scellrun.analyze import StageFailure

    with pytest.raises(StageFailure):
        run_analyze(
            planted_h5ad,
            run_dir=run_dir,
            resolutions=(0.5,),
            method="none",
            use_ai=False,
            force=False,
        )
    # with force it works
    run_analyze(
        planted_h5ad,
        run_dir=run_dir,
        resolutions=(0.5,),
        method="none",
        use_ai=False,
        force=True,
    )
    assert (run_dir / "05_report" / "index.html").exists()
