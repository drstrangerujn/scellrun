"""
Test fixtures.

Tests use a small synthetic AnnData rather than scanpy.datasets.pbmc3k
to keep CI under 1 minute and avoid the 5MB download on every run.
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture
def synthetic_h5ad(tmp_path):
    """
    600 cells × 500 genes synthetic raw-count AnnData.
    Some genes named MT-* and HBB so qc_vars annotation has signal.
    """
    import anndata as ad

    rng = np.random.default_rng(42)
    n_cells, n_genes = 600, 500
    counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.int32)
    counts[:, :10] = rng.poisson(lam=20.0, size=(n_cells, 10))  # high-expressed genes
    a = ad.AnnData(X=counts.astype(np.float32))
    a.var_names = (
        [f"MT-CO{i}" for i in range(1, 6)]
        + [f"HBB{i}" for i in range(1, 4)]
        + [f"RPS{i}" for i in range(1, 11)]
        + [f"GENE{i}" for i in range(n_genes - 18)]
    )
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    out = tmp_path / "synthetic.h5ad"
    a.write_h5ad(out)
    return out
