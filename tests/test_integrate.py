"""Tests for scellrun.scrna.integrate."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scellrun.scrna.integrate import IntegrationError, run_integrate, write_artifacts


@pytest.fixture
def post_qc_synthetic(tmp_path):
    """Synthetic h5ad with scellrun_qc_pass + a 'sample' batch key."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 500
    counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.int32)
    a = ad.AnnData(X=counts.astype(np.float32))
    a.var_names = [f"GENE{i}" for i in range(n_genes)]
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["scellrun_qc_pass"] = True
    a.obs["sample"] = ["A"] * 100 + ["B"] * 100
    return a


def test_integrate_runs_and_clusters(post_qc_synthetic):
    result, integrated = run_integrate(
        post_qc_synthetic,
        method="none",  # skip harmony to keep test deps light
        n_pcs=10,
        resolutions=(0.3, 0.5),
    )
    assert result.n_cells_used == 200
    assert "X_pca" in integrated.obsm
    assert "X_umap" in integrated.obsm
    assert "leiden_res_0_3" in integrated.obs.columns
    assert "leiden_res_0_5" in integrated.obs.columns


def test_integrate_drops_qc_fail(post_qc_synthetic):
    post_qc_synthetic.obs["scellrun_qc_pass"] = [True] * 150 + [False] * 50
    result, _ = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.5,),
        drop_qc_fail=True,
    )
    assert result.n_cells_used == 150


def test_integrate_keeps_qc_fail_when_asked(post_qc_synthetic):
    post_qc_synthetic.obs["scellrun_qc_pass"] = [True] * 150 + [False] * 50
    result, _ = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.5,),
        drop_qc_fail=False,
    )
    assert result.n_cells_used == 200


def test_integrate_detects_sample_key(post_qc_synthetic):
    result, _ = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.5,),
    )
    assert result.sample_key == "sample"


def test_integrate_stores_raw_before_scale(post_qc_synthetic):
    """v0.3+: .raw must be set so v0.3 markers reads log-normalized values, not scaled."""
    _, integrated = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.5,),
    )
    assert integrated.raw is not None
    assert integrated.raw.X.max() < 30  # log-normalized values, not raw counts


def test_integrate_rpca_not_implemented(post_qc_synthetic):
    """rpca and cca aren't no-ops anymore — they raise."""
    with pytest.raises(NotImplementedError, match="rpca"):
        run_integrate(post_qc_synthetic, method="rpca", n_pcs=10, resolutions=(0.5,))
    with pytest.raises(NotImplementedError, match="cca"):
        run_integrate(post_qc_synthetic, method="cca", n_pcs=10, resolutions=(0.5,))


def test_integrate_unknown_method_raises(post_qc_synthetic):
    with pytest.raises(ValueError, match="unknown method"):
        run_integrate(post_qc_synthetic, method="bogus", n_pcs=10, resolutions=(0.5,))


def test_integrate_harmony_no_sample_key_raises():
    """Explicit --method harmony without a detectable sample key must hard-fail."""
    rng = np.random.default_rng(0)
    a = ad.AnnData(X=rng.poisson(lam=2.0, size=(50, 100)).astype(np.float32))
    a.var_names = [f"GENE{i}" for i in range(100)]
    a.obs_names = [f"cell{i}" for i in range(50)]
    with pytest.raises(IntegrationError, match="sample"):
        run_integrate(a, method="harmony", n_pcs=10, resolutions=(0.5,), drop_qc_fail=False)


def test_integrate_harmony_single_level_raises():
    """harmony with only one batch level is meaningless; must hard-fail.

    v1.3.1+ note: a single-level sample column is now treated identically
    to a missing sample column by ``select_effective_sample_key``, so
    harmony falls into the "no sample/batch key" branch instead of a
    dedicated "multiple levels" branch. Behavior is the same — a hard
    fail before integration runs — only the error wording differs.
    """
    rng = np.random.default_rng(0)
    a = ad.AnnData(X=rng.poisson(lam=2.0, size=(50, 100)).astype(np.float32))
    a.var_names = [f"GENE{i}" for i in range(100)]
    a.obs_names = [f"cell{i}" for i in range(50)]
    a.obs["sample"] = ["only_one"] * 50
    with pytest.raises(IntegrationError, match="sample/batch key"):
        run_integrate(a, method="harmony", n_pcs=10, resolutions=(0.5,), drop_qc_fail=False)


def test_select_effective_sample_key_skips_single_level():
    """select_effective_sample_key must skip 1-level columns (v1.3.1)."""
    from scellrun.scrna.integrate import select_effective_sample_key

    df = pd.DataFrame({
        "orig.ident": ["A"] * 10,           # only 1 level
        "sample": ["P1", "P2"] * 5,         # 2 levels — picked
        "donor": ["X", "Y", "Z"] * 3 + ["X"],  # 3 levels but later in order
    })
    key, present = select_effective_sample_key(df)
    assert key == "sample"
    assert set(present) == {"orig.ident", "sample", "donor"}


def test_select_effective_sample_key_all_single_level_returns_none():
    """All-single-level candidates → key is None, present is non-empty."""
    from scellrun.scrna.integrate import select_effective_sample_key

    df = pd.DataFrame({"orig.ident": ["A"] * 5, "sample": ["B"] * 5})
    key, present = select_effective_sample_key(df)
    assert key is None
    assert set(present) == {"orig.ident", "sample"}


def test_select_effective_sample_key_no_candidates():
    """No candidate column at all → both None and empty list."""
    from scellrun.scrna.integrate import select_effective_sample_key

    df = pd.DataFrame({"unrelated_col": ["X"] * 5})
    key, present = select_effective_sample_key(df)
    assert key is None
    assert present == []


def test_integrate_subsets_to_hvg_for_scale(post_qc_synthetic):
    """v1.3.1: scale + PCA must run on HVG-subset .X (preserving .raw).

    Closes the v0.7 dogfood follow-up where ``sc.pp.scale`` was paying the
    full-matrix memory cost for genes PCA was about to ignore via the
    ``mask_var="highly_variable"`` parameter.
    """
    result, integrated = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.3,),
    )
    n_top_genes_default = 2000
    n_full_vars = 500  # post_qc_synthetic fixture uses 500 genes
    expected_n_hvg = min(n_top_genes_default, n_full_vars)
    assert integrated.n_vars == expected_n_hvg, (
        f"expected .X subset to {expected_n_hvg} HVGs, got {integrated.n_vars}"
    )
    # .raw must retain the full pre-subset gene set so downstream markers /
    # annotate (use_raw=True) still see every gene.
    assert integrated.raw is not None
    assert integrated.raw.n_vars == n_full_vars


def test_integrate_records_scale_subset_decision(post_qc_synthetic, tmp_path):
    """The scale_subset_to_hvg decision row must land in 00_decisions.jsonl."""
    import json

    run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.3,),
        run_dir=tmp_path,
    )
    decisions_path = tmp_path / "00_decisions.jsonl"
    assert decisions_path.exists()
    rows = [json.loads(line) for line in decisions_path.read_text().splitlines() if line.strip()]
    scale_rows = [r for r in rows if r.get("key") == "scale_subset_to_hvg"]
    assert len(scale_rows) >= 1, "scale_subset_to_hvg decision row missing"
    row = scale_rows[-1]
    assert row["value"] is True
    assert row["source"] == "auto"
    assert row["stage"] == "integrate"


def test_integrate_cc_regress_no_match_warns(capsys, post_qc_synthetic):
    """If CC genes don't match adata.var_names, warn loudly and skip CC regression."""
    # synthetic h5ad has GENE0..GENE499, no Tirosh genes match
    run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.5,),
        regress_cell_cycle=True,
        species="human",
    )
    captured = capsys.readouterr()
    assert "WARNING" in captured.out
    assert "Cell-cycle regression SKIPPED" in captured.out


def test_write_artifacts_produces_report_and_umap(post_qc_synthetic, tmp_path):
    result, integrated = run_integrate(
        post_qc_synthetic,
        method="none",
        n_pcs=10,
        resolutions=(0.3, 0.5),
    )
    out = tmp_path / "out"
    artifacts = write_artifacts(result, integrated, out)
    assert artifacts["report"].exists()
    assert artifacts["umap_grid"].exists()
    assert artifacts["integrated_h5ad"].exists()
    # Roundtrip: integrated h5ad reloads with cluster cols intact
    reloaded = ad.read_h5ad(artifacts["integrated_h5ad"])
    assert "leiden_res_0_3" in reloaded.obs.columns
