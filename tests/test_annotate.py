"""Tests for scellrun.scrna.annotate."""
from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest

from scellrun.scrna.annotate import (
    _best_panel_match,
    _score_cluster_against_panel,
    _select_panel,
    run_annotate,
    write_artifacts,
)


def test_score_cluster_against_panel_overlap_fraction():
    panel = {"A": ["G1", "G2", "G3"], "B": ["G4", "G5"]}
    cluster_top = ["G1", "G2", "G99"]
    s = _score_cluster_against_panel(cluster_top, panel)
    assert s["A"] == pytest.approx(2 / 3)
    assert s["B"] == 0.0


def test_best_panel_match_returns_label_and_margin():
    label, score, margin = _best_panel_match({"X": 0.8, "Y": 0.5, "Z": 0.1})
    assert label == "X"
    assert score == 0.8
    assert margin == pytest.approx(0.3)


def test_best_panel_match_empty_returns_unassigned():
    label, score, margin = _best_panel_match({})
    assert label == "Unassigned"


def test_select_panel_prefers_chondrocyte_then_broad():
    mod = SimpleNamespace(
        __name__="scellrun.profiles.test",
        celltype_broad={"X": ["G1"]},
        chondrocyte_markers={"Y": ["G2"]},
    )
    name, panel = _select_panel(mod, None)
    assert name == "chondrocyte_markers"


def test_select_panel_falls_back_to_broad():
    mod = SimpleNamespace(
        __name__="scellrun.profiles.test",
        celltype_broad={"X": ["G1"]},
    )
    name, panel = _select_panel(mod, None)
    assert name == "celltype_broad"


def test_select_panel_explicit_override():
    mod = SimpleNamespace(
        __name__="scellrun.profiles.test",
        celltype_broad={"X": ["G1"]},
        my_panel={"Y": ["G2"]},
    )
    name, panel = _select_panel(mod, "my_panel")
    assert name == "my_panel"


def test_select_panel_no_panel_raises():
    mod = SimpleNamespace(__name__="scellrun.profiles.test")
    with pytest.raises(ValueError, match="no marker panels"):
        _select_panel(mod, None)


@pytest.fixture
def adata_with_panel_signal(tmp_path):
    """
    Two-cluster synthetic AnnData where cluster 0 expresses GENE_A_1/2/3
    (matching panel A) and cluster 1 expresses GENE_B_1/2/3 (matching panel B).
    """
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 50
    base = rng.normal(loc=0.5, scale=0.2, size=(n_cells, n_genes)).clip(0, None)
    # cluster 0 (rows 0-99): high GENE_A_1, _2, _3 (cols 0-2)
    base[:100, 0:3] += 5.0
    # cluster 1 (rows 100-199): high GENE_B_1, _2, _3 (cols 3-5)
    base[100:, 3:6] += 5.0

    a = ad.AnnData(X=base.astype(np.float32))
    a.var_names = (
        ["GENE_A_1", "GENE_A_2", "GENE_A_3", "GENE_B_1", "GENE_B_2", "GENE_B_3"]
        + [f"GENE_{i}" for i in range(n_genes - 6)]
    )
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["leiden_res_0_5"] = (["0"] * 100 + ["1"] * 100)
    a.obs["leiden_res_0_5"] = a.obs["leiden_res_0_5"].astype("category")
    return a


def test_run_annotate_picks_correct_label(adata_with_panel_signal):
    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={
            "TypeA": ["GENE_A_1", "GENE_A_2", "GENE_A_3"],
            "TypeB": ["GENE_B_1", "GENE_B_2", "GENE_B_3"],
        },
    )
    result = run_annotate(
        adata_with_panel_signal,
        profile_module,
        resolution=0.5,
        use_ai=False,
        use_pubmed=False,
    )
    by_cluster = {a.cluster: a.panel_label for a in result.annotations}
    assert by_cluster["0"] == "TypeA"
    assert by_cluster["1"] == "TypeB"


def test_run_annotate_unknown_resolution_raises(adata_with_panel_signal):
    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={"X": ["G1"]},
    )
    with pytest.raises(ValueError, match="not present"):
        run_annotate(adata_with_panel_signal, profile_module, resolution=99.0)


def test_run_annotate_writes_obs_column(adata_with_panel_signal):
    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={
            "TypeA": ["GENE_A_1", "GENE_A_2", "GENE_A_3"],
            "TypeB": ["GENE_B_1", "GENE_B_2", "GENE_B_3"],
        },
    )
    run_annotate(adata_with_panel_signal, profile_module, resolution=0.5)
    assert "scellrun_celltype_panel" in adata_with_panel_signal.obs.columns


def test_write_artifacts_produces_csv_and_report(adata_with_panel_signal, tmp_path):
    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={
            "TypeA": ["GENE_A_1", "GENE_A_2", "GENE_A_3"],
            "TypeB": ["GENE_B_1", "GENE_B_2", "GENE_B_3"],
        },
    )
    result = run_annotate(adata_with_panel_signal, profile_module, resolution=0.5)
    out = tmp_path / "out"
    artifacts = write_artifacts(result, adata_with_panel_signal, out)
    assert artifacts["annotations"].exists()
    assert artifacts["report"].exists()
    assert artifacts["annotated_h5ad"].exists()
