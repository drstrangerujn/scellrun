"""Tests for `_autopick_panel_for_data` panel-selection heuristic.

The v1.1.0 cold-validation surfaced a trap where the chondrocyte panel
picked up 1-2 weak cluster hits on immune-rich subchondral-bone data and
was kept by default — even though the broad panel had 7+ cluster hits.
The 1.5× margin requirement closes that gap. These tests guard the
margin so future tweaks don't silently regress to "any chondrocyte hit
keeps chondrocyte".
"""
from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scellrun.analyze import (
    PANEL_AUTOPICK_CHONDRO_MARGIN,
    _autopick_panel_for_data,
)


def _make_profile(*, chondro: dict | None = None, broad: dict | None = None) -> SimpleNamespace:
    """Mock a profile module with the marker_panel attributes the picker reads."""
    attrs: dict[str, dict] = {}
    if chondro is not None:
        attrs["chondrocyte_markers"] = chondro
    if broad is not None:
        attrs["celltype_broad"] = broad
    return SimpleNamespace(**attrs)


def _adata_with_clusters(n_clusters: int, top_markers_per_cluster: list[list[str]]) -> ad.AnnData:
    """Synthesize a minimal AnnData carrying the post-integrate state the
    panel picker reads:

    - obs[``leiden_res_0_5``] (=== `chosen_res=0.5`)
    - var_names spanning every gene the markers refer to
    - rank_genes_groups already populated under the picker's expected key
    """
    chosen_res = 0.5
    res_key = f"leiden_res_{chosen_res:g}".replace(".", "_")
    rank_key = f"rank_panel_pick_{res_key}"

    n_cells_per_cluster = 10
    n_cells = n_clusters * n_cells_per_cluster

    all_genes = sorted({g for cluster_top in top_markers_per_cluster for g in cluster_top})
    # add filler so var_names is non-trivial
    filler = [f"FILLER{i:04d}" for i in range(50)]
    all_genes = all_genes + filler
    n_genes = len(all_genes)

    rng = np.random.default_rng(0)
    X = rng.normal(loc=0.5, scale=0.3, size=(n_cells, n_genes)).clip(0, None).astype(np.float32)
    a = ad.AnnData(X=X)
    a.var_names = all_genes
    a.obs_names = [f"c{i:04d}" for i in range(n_cells)]
    cluster_labels = []
    for ci in range(n_clusters):
        cluster_labels.extend([str(ci)] * n_cells_per_cluster)
    a.obs[res_key] = pd.Categorical(cluster_labels)

    # Pre-populate rank_genes_groups so the picker doesn't try to recompute
    # it (and so we can guarantee the top markers per cluster).
    cluster_names = [str(ci) for ci in range(n_clusters)]
    max_top = max(len(m) for m in top_markers_per_cluster) if top_markers_per_cluster else 0
    # 30 is the picker's window; pad with filler genes that aren't in any panel.
    top_n = max(30, max_top)
    names_arr = np.zeros(top_n, dtype=[(c, "U64") for c in cluster_names])
    for ci, top in enumerate(top_markers_per_cluster):
        padded = list(top) + [f"PAD{ci}_{j}" for j in range(top_n - len(top))]
        for i, g in enumerate(padded[:top_n]):
            names_arr[i][cluster_names[ci]] = g
    a.uns[rank_key] = {"names": names_arr}
    return a


def test_autopick_only_broad_present_returns_broad():
    profile = _make_profile(broad={"Macrophage": ["CD68", "CD14"]})
    adata = _adata_with_clusters(2, [["CD68"], ["CD14"]])
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "celltype_broad"
    assert "only celltype_broad" in reason


def test_autopick_only_chondro_present_returns_chondro():
    profile = _make_profile(chondro={"ProC": ["COL2A1", "ACAN"]})
    adata = _adata_with_clusters(2, [["COL2A1"], ["ACAN"]])
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "chondrocyte_markers"
    assert "only chondrocyte_markers" in reason


def test_autopick_chondro_clears_margin_kept():
    """3 chondro hits vs 2 broad hits → 1.5× exact → kept (>= margin)."""
    profile = _make_profile(
        chondro={"ProC": ["COL2A1", "ACAN", "SOX9"]},
        broad={"Macrophage": ["CD68"], "Endothelial": ["PECAM1"]},
    )
    # 5 clusters: 3 chondro-marker hits, 2 broad-marker hits, exclusive
    adata = _adata_with_clusters(
        5, [["COL2A1"], ["ACAN"], ["SOX9"], ["CD68"], ["PECAM1"]]
    )
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "chondrocyte_markers"
    assert "kept" in reason
    assert "chondrocyte_hits=3" in reason
    assert "broad_hits=2" in reason


def test_autopick_chondro_below_margin_swaps_to_broad():
    """v1.1.0 cold-validation gap: chondro 7 vs broad 5 = 1.4× < 1.5× → swap."""
    profile = _make_profile(
        chondro={"ProC": ["COL2A1", "ACAN"]},
        broad={"Mac": ["CD68"], "Endo": ["PECAM1"]},
    )
    # 7 chondro hits, 5 broad hits → ratio 1.4 < 1.5 → swap to broad
    adata = _adata_with_clusters(
        12,
        [
            # 7 clusters with chondro markers
            ["COL2A1"], ["COL2A1"], ["COL2A1"], ["ACAN"], ["ACAN"], ["ACAN"], ["ACAN"],
            # 5 clusters with broad markers
            ["CD68"], ["CD68"], ["PECAM1"], ["PECAM1"], ["PECAM1"],
        ],
    )
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "celltype_broad"
    assert "swapped" in reason
    assert f">={PANEL_AUTOPICK_CHONDRO_MARGIN:g}x margin" in reason


def test_autopick_immune_rich_data_swaps_to_broad():
    """v1.1.0 BML_1 reproducer: chondro panel 1 hit, broad panel 9 hits → swap."""
    profile = _make_profile(
        chondro={"ProC": ["COL2A1"]},
        broad={"Mac": ["CD68"], "Endo": ["PECAM1"], "T": ["CD3E"]},
    )
    adata = _adata_with_clusters(
        10,
        [
            ["COL2A1"],  # 1 chondro hit
            ["CD68"], ["CD68"], ["CD68"],
            ["PECAM1"], ["PECAM1"], ["PECAM1"],
            ["CD3E"], ["CD3E"], ["CD3E"],
        ],
    )
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "celltype_broad"
    assert "chondrocyte_hits=1" in reason
    assert "broad_hits=9" in reason


def test_autopick_zero_broad_hits_keeps_chondro_by_inf_ratio():
    """Edge: broad panel matches no cluster → chondro / broad ratio is +inf
    → the chondrocyte panel is kept regardless of how few hits it has."""
    profile = _make_profile(
        chondro={"ProC": ["COL2A1"]},
        broad={"NonMatching": ["GENE_NEVER_PRESENT"]},
    )
    adata = _adata_with_clusters(2, [["COL2A1"], ["FILLER0001"]])
    panel, reason = _autopick_panel_for_data(profile, adata, chosen_res=0.5)
    assert panel == "chondrocyte_markers"
    assert "kept" in reason


def test_autopick_missing_leiden_col_falls_back_to_chondro():
    """If the chosen resolution's leiden column isn't in obs, default to
    chondrocyte_markers (current pre-1.1.0 behavior preserved)."""
    profile = _make_profile(
        chondro={"ProC": ["COL2A1"]},
        broad={"Mac": ["CD68"]},
    )
    a = ad.AnnData(X=np.zeros((10, 5), dtype=np.float32))
    a.var_names = ["A", "B", "C", "D", "E"]
    a.obs_names = [f"c{i}" for i in range(10)]
    # No leiden_res_0_5 column!
    panel, reason = _autopick_panel_for_data(profile, a, chosen_res=0.5)
    assert panel == "chondrocyte_markers"
    assert "no leiden col" in reason


def test_autopick_no_panels_at_all_returns_empty():
    """Profile with neither panel → empty string panel name."""
    profile = _make_profile()  # no attrs
    a = ad.AnnData(X=np.zeros((4, 3), dtype=np.float32))
    a.var_names = ["A", "B", "C"]
    a.obs_names = [f"c{i}" for i in range(4)]
    panel, _ = _autopick_panel_for_data(profile, a, chosen_res=0.5)
    assert panel == ""
