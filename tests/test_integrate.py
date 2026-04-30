"""Tests for scellrun.scrna.integrate."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from scellrun.scrna.integrate import run_integrate, write_artifacts


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
