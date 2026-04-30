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


@pytest.fixture
def adata_immune_rich_BML():
    """
    v1.1.0: synthetic BML_1-like fixture for the panel-autopick swap test.

    13 clusters, with top markers chosen to mimic the cold-validation BML_1
    profile: predominantly immune (T/NK/B/plasma/mast/macrophage/dendritic)
    plus a couple of stromal clusters and one fibroblast / chondrocyte signal.
    The autopick should count chondrocyte hits << 1.5x broad hits and swap
    to celltype_broad.

    Built without scanpy: we plant a fake rank_genes_groups uns dict so
    `_autopick_panel_for_data` can read top markers without rerunning ranking.
    """
    rng = np.random.default_rng(0)
    n_cells, n_genes = 130, 50
    a = ad.AnnData(X=rng.normal(loc=0.5, scale=0.2, size=(n_cells, n_genes)).clip(0, None).astype(np.float32))
    a.var_names = [f"GENE_{i}" for i in range(n_genes)]
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    # 13 clusters x 10 cells each.
    cluster_labels: list[str] = []
    for i in range(13):
        cluster_labels.extend([str(i)] * 10)
    a.obs["leiden_res_0_3"] = cluster_labels
    a.obs["leiden_res_0_3"] = a.obs["leiden_res_0_3"].astype("category")
    return a


def _plant_fake_rank_genes(adata, res_key, top_per_cluster):
    """
    Plant a minimal fake rank_genes_groups uns dict so the autopick step
    can read top markers without invoking scanpy's wilcoxon path.
    """
    cluster_names = list(top_per_cluster.keys())
    max_top = max(len(v) for v in top_per_cluster.values())
    # rank_genes_groups stores names as a structured array with one column per
    # cluster. Build that shape so `rg["names"][c][:30]` works downstream.
    dtype = np.dtype([(c, "U64") for c in cluster_names])
    names = np.zeros(max_top, dtype=dtype)
    for c, genes in top_per_cluster.items():
        padded = list(genes) + [""] * (max_top - len(genes))
        names[c] = padded
    rank_key = f"rank_panel_pick_{res_key}"
    adata.uns[rank_key] = {"names": names}


def test_autopick_swaps_to_broad_on_immune_rich_data(adata_immune_rich_BML, tmp_path):
    """
    v1.1.0 gap 1: on a BML_1-shaped immune-rich dataset (13 clusters, mostly
    T/NK/B/plasma/mast/macrophage), `_autopick_panel_for_data` swaps from
    chondrocyte_markers to celltype_broad and the rationale carries the
    margin explanation with hit counts.
    """
    from scellrun.analyze import _autopick_panel_for_data
    from scellrun.profiles import load as load_profile

    prof = load_profile("joint-disease")

    # Plant top markers per cluster mimicking BML_1's actual labels.
    # Most are immune (broad-panel hits, no chondrocyte_markers genes).
    # One cluster is fibroblast-like which hits FC of chondrocyte_markers
    # via COL1A1/COL1A2 — that's the one chondrocyte-marker hit in the data.
    top_per_cluster = {
        "0":  ["RAMP2", "ENG", "A2M", "PECAM1", "VWF", "CDH5"],            # endothelial (broad)
        "1":  ["LYZ", "AIF1", "HLA-DRA", "TYROBP", "CD68", "CD163"],       # macrophage (broad)
        "2":  ["CD3D", "CD3E", "IL32", "CCL5", "B2M"],                     # T cell (broad)
        "3":  ["NKG7", "KLRD1", "GZMA", "GNLY", "CD8A"],                   # NK / CD8 (broad)
        "4":  ["COL1A2", "DCN", "COL1A1", "MMP2"],                         # fibroblast — chondrocyte FC also fits
        "5":  ["RGS5", "MCAM", "ACTA2", "MYL9", "TAGLN"],                  # pericyte / smooth muscle (broad)
        "6":  ["CD37", "CD79A", "MS4A1", "CD19"],                          # B cell (broad)
        "7":  ["IRF7", "IRF8", "IL3RA", "PLD4"],                           # plasmacytoid DC (broad)
        "8":  ["HLA-DPA1", "CD74", "HLA-DRB1"],                            # generic activated immune
        "9":  ["ACP5", "CTSK", "MMP9", "CD68"],                            # osteoclast (broad)
        "10": ["KLRB1", "IL7R", "CD3D"],                                   # memory T (broad)
        "11": ["MZB1", "JCHAIN", "XBP1", "IGHG1", "CD79A"],                # plasma cell (broad)
        "12": ["CPA3", "TPSAB1", "TPSB2", "MS4A2", "CTSG"],                # mast cell (broad)
    }
    _plant_fake_rank_genes(adata_immune_rich_BML, "leiden_res_0_3", top_per_cluster)

    panel, rationale = _autopick_panel_for_data(prof, adata_immune_rich_BML, 0.3)

    assert panel == "celltype_broad", (
        f"expected swap to celltype_broad on immune-rich fixture, got {panel!r} "
        f"(rationale: {rationale})"
    )
    # Rationale must carry the new tie-break explanation with hit counts.
    assert "chondrocyte_hits=" in rationale
    assert "broad_hits=" in rationale
    assert "1.5x" in rationale
    assert "celltype_broad" in rationale.lower()


def test_autopick_keeps_chondrocyte_with_clear_dominance(tmp_path):
    """
    v1.1.0 gap 1, opposite case: when chondrocyte hits clearly clear the
    1.5x margin (fixture is mostly chondrocyte), the panel is kept and
    the rationale says so.
    """
    from scellrun.analyze import _autopick_panel_for_data
    from scellrun.profiles import load as load_profile

    prof = load_profile("joint-disease")

    rng = np.random.default_rng(1)
    n_cells, n_genes = 50, 30
    a = ad.AnnData(X=rng.normal(loc=0.5, scale=0.2, size=(n_cells, n_genes)).clip(0, None).astype(np.float32))
    a.var_names = [f"GENE_{i}" for i in range(n_genes)]
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["leiden_res_0_3"] = (["0"] * 10 + ["1"] * 10 + ["2"] * 10 + ["3"] * 10 + ["4"] * 10)
    a.obs["leiden_res_0_3"] = a.obs["leiden_res_0_3"].astype("category")

    # 4/5 clusters hit chondrocyte panel; 1/5 hits broad. 4 vs 1 = 4x margin.
    top_per_cluster = {
        "0": ["BMP2", "C11orf96", "HMGA1"],            # ProC chondrocyte
        "1": ["FRZB", "CHRDL2", "CYTL1"],              # EC chondrocyte
        "2": ["HSPA1A", "HSPA1B", "DDIT3"],            # HomC chondrocyte
        "3": ["COL10A1", "SPP1", "IBSP"],              # HTC chondrocyte
        "4": ["CD3D", "CD3E", "IL32"],                 # T cell (broad only)
    }
    _plant_fake_rank_genes(a, "leiden_res_0_3", top_per_cluster)

    panel, rationale = _autopick_panel_for_data(prof, a, 0.3)

    assert panel == "chondrocyte_markers"
    assert "kept chondrocyte_markers" in rationale
    assert "chondrocyte_hits=" in rationale
    assert "broad_hits=" in rationale


def test_autopick_decision_row_carries_margin_rationale(adata_immune_rich_BML, tmp_path):
    """
    v1.1.0 gap 1: end-to-end check that the orchestrator's
    `analyze.annotate.auto_panel` decision row's rationale string contains
    the new margin explanation, not the legacy "fine-subtype panel preferred"
    fixed string.
    """
    from scellrun.analyze import _autopick_panel_for_data
    from scellrun.profiles import load as load_profile

    prof = load_profile("joint-disease")

    # Recycle the fixture's planted markers and force the swap branch.
    top_per_cluster = {
        "0":  ["RAMP2", "ENG", "A2M"],
        "1":  ["LYZ", "AIF1", "HLA-DRA"],
        "2":  ["CD3D", "CD3E", "IL32"],
        "3":  ["NKG7", "KLRD1", "GZMA"],
        "4":  ["COL1A2", "DCN"],
        "5":  ["RGS5", "MCAM"],
        "6":  ["CD37", "CD79A"],
        "7":  ["IRF7", "IL3RA"],
        "8":  ["HLA-DPA1", "CD74"],
        "9":  ["ACP5", "CTSK"],
        "10": ["KLRB1", "IL7R"],
        "11": ["MZB1", "JCHAIN"],
        "12": ["CPA3", "TPSAB1"],
    }
    _plant_fake_rank_genes(adata_immune_rich_BML, "leiden_res_0_3", top_per_cluster)
    panel, rationale = _autopick_panel_for_data(prof, adata_immune_rich_BML, 0.3)

    assert panel == "celltype_broad"
    # The exact margin explanation pattern that the orchestrator records
    # via its `analyze.annotate.auto_panel` decision row.
    assert "swapped to celltype_broad" in rationale
    assert "required >=1.5x margin to keep chondrocyte panel" in rationale


def test_run_annotate_orchestrator_panel_is_auto_not_user(adata_with_panel_signal, tmp_path):
    """
    v1.1.0 gap 3: when the orchestrator passes a panel name with
    panel_name_user_supplied=False, the recorded decision row reads
    source="auto", not source="user".
    """
    from scellrun.decisions import read_decisions

    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={
            "TypeA": ["GENE_A_1", "GENE_A_2", "GENE_A_3"],
            "TypeB": ["GENE_B_1", "GENE_B_2", "GENE_B_3"],
        },
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    run_annotate(
        adata_with_panel_signal,
        profile_module,
        resolution=0.5,
        panel_name="chondrocyte_markers",  # orchestrator-injected
        panel_name_user_supplied=False,    # but NOT a user override
        run_dir=run_dir,
    )

    rows = read_decisions(run_dir)
    panel_rows = [r for r in rows if r["stage"] == "annotate" and r["key"] == "panel"]
    assert len(panel_rows) == 1
    assert panel_rows[0]["source"] == "auto", (
        f"expected source='auto' for orchestrator-injected panel, "
        f"got {panel_rows[0]['source']!r}"
    )


def test_run_annotate_user_panel_is_user(adata_with_panel_signal, tmp_path):
    """
    v1.1.0 gap 3: when the user explicitly passes --panel via the per-stage
    CLI, the per-stage CLI sets panel_name_user_supplied=True and the row
    reads source="user".
    """
    from scellrun.decisions import read_decisions

    profile_module = SimpleNamespace(
        __name__="scellrun.profiles.test",
        chondrocyte_markers={
            "TypeA": ["GENE_A_1", "GENE_A_2", "GENE_A_3"],
            "TypeB": ["GENE_B_1", "GENE_B_2", "GENE_B_3"],
        },
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    run_annotate(
        adata_with_panel_signal,
        profile_module,
        resolution=0.5,
        panel_name="chondrocyte_markers",  # user-typed
        panel_name_user_supplied=True,     # actually a user override
        run_dir=run_dir,
    )

    rows = read_decisions(run_dir)
    panel_rows = [r for r in rows if r["stage"] == "annotate" and r["key"] == "panel"]
    assert len(panel_rows) == 1
    assert panel_rows[0]["source"] == "user"


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
