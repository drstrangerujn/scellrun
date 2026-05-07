"""Tests for scellrun.scrna.markers."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from scellrun.scrna.markers import run_markers, write_artifacts


@pytest.fixture
def integrated_synthetic(tmp_path):
    """
    Synthetic AnnData mimicking what `scrna integrate` writes:
    - already log-normalized X
    - leiden_res_0_3 + leiden_res_0_5 obs columns with 2-3 clusters each
    """
    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 100

    # Build two cell groups with one differentially expressed gene each
    base = rng.normal(loc=0.5, scale=0.3, size=(n_cells, n_genes)).clip(0, None)
    base[:100, 0] += 3.0  # gene 0 high in group A
    base[100:, 1] += 3.0  # gene 1 high in group B

    a = ad.AnnData(X=base.astype(np.float32))
    a.var_names = [f"GENE{i}" for i in range(n_genes)]
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["leiden_res_0_3"] = (["0"] * 100 + ["1"] * 100)
    a.obs["leiden_res_0_5"] = (["0"] * 50 + ["1"] * 50 + ["2"] * 100)
    a.obs["leiden_res_0_3"] = a.obs["leiden_res_0_3"].astype("category")
    a.obs["leiden_res_0_5"] = a.obs["leiden_res_0_5"].astype("category")
    return a


def test_markers_runs_on_all_resolutions(integrated_synthetic):
    result, per_res = run_markers(integrated_synthetic, only_positive=True)
    assert set(result.resolutions) == {0.3, 0.5}
    assert 0.3 in per_res
    assert 0.5 in per_res
    assert len(per_res[0.3]) > 0


def test_markers_picks_signal_gene(integrated_synthetic):
    """Cluster 0 at res=0.3 should have GENE0 as a top marker (we baked it in)."""
    result, per_res = run_markers(integrated_synthetic, only_positive=True, logfc_threshold=0.5)
    df_03 = per_res[0.3]
    cluster_0_top = df_03[df_03["cluster"] == "0"].sort_values("log2fc", ascending=False).head(3)
    assert "GENE0" in cluster_0_top["gene"].values


def test_markers_filter_by_resolution(integrated_synthetic):
    result, per_res = run_markers(integrated_synthetic, resolutions=(0.5,))
    assert tuple(result.resolutions) == (0.5,)
    assert 0.3 not in per_res


def test_markers_unknown_resolution_raises(integrated_synthetic):
    with pytest.raises(ValueError):
        run_markers(integrated_synthetic, resolutions=(2.5,))


def test_markers_no_leiden_columns_raises():
    a = ad.AnnData(X=np.zeros((10, 5), dtype=np.float32))
    a.obs["something_else"] = ["x"] * 10
    with pytest.raises(ValueError, match="leiden_res"):
        run_markers(a)


def test_write_artifacts_produces_csv_and_report(integrated_synthetic, tmp_path):
    result, per_res = run_markers(integrated_synthetic, only_positive=True)
    out = tmp_path / "out"
    artifacts = write_artifacts(result, per_res, out)
    assert artifacts["report"].exists()
    assert (out / "markers_res_0.3.csv").exists()
    assert (out / "markers_res_0.5.csv").exists()


def test_markers_caps_at_default_n_per_cluster():
    """v1.3.1: rank_genes_groups must cap per-cluster output to bound memory.

    With a 1500-gene fixture, each cluster's pre-filter row count must be
    <= DEFAULT_RANK_GENES_N_PER_CLUSTER. Using ``only_positive=False`` and
    very loose ``logfc_threshold`` / ``pct_min`` keeps the markers DataFrame
    in proportion with the rank cap so we can assert directly on it.
    """
    from scellrun.scrna.markers import DEFAULT_RANK_GENES_N_PER_CLUSTER

    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 1500
    X = rng.normal(loc=0.5, scale=0.3, size=(n_cells, n_genes)).clip(0, None).astype(np.float32)
    a = ad.AnnData(X=X)
    a.var_names = [f"G{i:05d}" for i in range(n_genes)]
    a.obs_names = [f"c{i:04d}" for i in range(n_cells)]
    a.obs["leiden_res_0_5"] = (["0"] * 100 + ["1"] * 100)
    a.obs["leiden_res_0_5"] = a.obs["leiden_res_0_5"].astype("category")

    _, per_res = run_markers(
        a,
        only_positive=False,
        logfc_threshold=0.0,
        pct_min=0.0,
    )
    df = per_res[0.5]
    # 2 clusters × cap → max possible rows
    assert len(df) <= 2 * DEFAULT_RANK_GENES_N_PER_CLUSTER
    # Per-cluster row count must respect the cap
    for cluster in df["cluster"].unique():
        n_rows_cluster = (df["cluster"] == cluster).sum()
        assert n_rows_cluster <= DEFAULT_RANK_GENES_N_PER_CLUSTER, (
            f"cluster {cluster} has {n_rows_cluster} rows, "
            f"exceeds cap {DEFAULT_RANK_GENES_N_PER_CLUSTER}"
        )
