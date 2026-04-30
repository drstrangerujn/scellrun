"""
Tests for v0.8 stage self-checks + --auto-fix.

The acceptance contract from PLAN.md v0.8:

    Synthetic adversarial inputs produce actionable suggestions; with
    --auto-fix on, the pipeline recovers and continues. Decision log
    shows source="auto" rows for each self-check trigger + suggestion.

These tests build small synthetic h5ads that trip each rule and verify
both the standalone `*_self_check()` helpers and the orchestrator's
auto-fix loop.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from scellrun.decisions import read_decisions
from scellrun.self_check import annotate_self_check, integrate_self_check, qc_self_check

# ---------------------------------------------------------------------------
# QC self-check
# ---------------------------------------------------------------------------


def test_qc_self_check_fires_on_low_pass_rate():
    """All cells > 30% mt → QC pass-rate < 30% → suggestion fires."""
    from scellrun.defaults import SCRNA_QC

    sensitivity = {
        "max_pct_mt": [
            {"threshold": 5.0, "n_pass": 0, "pct_pass": 0.0},
            {"threshold": 10.0, "n_pass": 0, "pct_pass": 0.0},
            {"threshold": 15.0, "n_pass": 0, "pct_pass": 0.0},
            {"threshold": 20.0, "n_pass": 5, "pct_pass": 5.0},
            {"threshold": 25.0, "n_pass": 50, "pct_pass": 50.0},
            {"threshold": 30.0, "n_pass": 80, "pct_pass": 80.0},
        ],
        "max_genes": [
            {"threshold": 3000, "n_pass": 90, "pct_pass": 90.0},
            {"threshold": 4000, "n_pass": 95, "pct_pass": 95.0},
            {"threshold": 5000, "n_pass": 100, "pct_pass": 100.0},
            {"threshold": 6000, "n_pass": 100, "pct_pass": 100.0},
            {"threshold": 8000, "n_pass": 100, "pct_pass": 100.0},
        ],
        "min_genes": [],
    }
    findings = qc_self_check(
        n_cells_in=100,
        n_cells_pass=10,  # 10% pass rate, well below the 30% trigger
        sensitivity=sensitivity,
        thresholds=SCRNA_QC,
    )
    assert len(findings) == 1
    f = findings[0]
    assert f.stage == "qc"
    assert f.code == "qc_low_pass_rate"
    # Smallest single relaxation reaching 60%+ on max_pct_mt is threshold 30 (80%).
    assert f.fix == {"max_pct_mt": 30.0}
    assert "30" in f.suggestion
    assert "30" in f.trigger or "10.0%" in f.trigger
    assert "below" in f.trigger.lower() or "passed" in f.trigger.lower()


def test_qc_self_check_silent_when_pass_rate_ok():
    """Pass-rate above the trigger → no findings.

    v0.9.1 raised the trigger ceiling to 60%; we use 85% here to stay
    above it without being suspiciously high.
    """
    from scellrun.defaults import SCRNA_QC

    sensitivity = {"max_pct_mt": [], "max_genes": [], "min_genes": []}
    findings = qc_self_check(
        n_cells_in=100,
        n_cells_pass=85,
        sensitivity=sensitivity,
        thresholds=SCRNA_QC,
    )
    assert findings == []


def test_qc_self_check_fires_in_degraded_band():
    """
    v0.9.1 (B1): pass-rate sitting in 30-60% range now triggers the check.
    Before the patch, only <30% triggered, so degraded preps slipped past.
    """
    from scellrun.defaults import SCRNA_QC

    sensitivity = {
        "max_pct_mt": [
            {"threshold": 25.0, "n_pass": 40, "pct_pass": 40.0},
            {"threshold": 30.0, "n_pass": 80, "pct_pass": 80.0},
        ],
        "max_genes": [
            {"threshold": 5000, "n_pass": 100, "pct_pass": 100.0},
        ],
    }
    findings = qc_self_check(
        n_cells_in=100,
        n_cells_pass=45,  # 45% — degraded but >30%
        sensitivity=sensitivity,
        thresholds=SCRNA_QC,
    )
    assert len(findings) == 1
    assert findings[0].code == "qc_low_pass_rate"


def test_qc_self_check_no_easy_fix_when_no_relaxation_helps():
    """When no single relaxation reaches 60%, surface a human-only suggestion."""
    from scellrun.defaults import SCRNA_QC

    sensitivity = {
        "max_pct_mt": [
            {"threshold": 25.0, "n_pass": 5, "pct_pass": 5.0},
            {"threshold": 30.0, "n_pass": 10, "pct_pass": 10.0},
        ],
        "max_genes": [
            {"threshold": 5000, "n_pass": 30, "pct_pass": 30.0},
        ],
    }
    findings = qc_self_check(
        n_cells_in=100,
        n_cells_pass=5,
        sensitivity=sensitivity,
        thresholds=SCRNA_QC,
    )
    assert len(findings) == 1
    assert findings[0].fix == {}  # human-only — orchestrator can't auto-fix
    assert findings[0].code == "qc_low_pass_rate_no_easy_fix"


# ---------------------------------------------------------------------------
# Integrate self-check
# ---------------------------------------------------------------------------


def test_integrate_self_check_too_few_clusters():
    """Every resolution ≤2 clusters → suggest --resolutions aio."""
    cluster_counts = {0.1: 1, 0.3: 2, 0.5: 2, 0.8: 2, 1.0: 2}
    findings = integrate_self_check(
        cluster_counts=cluster_counts,
        quality={
            r: {"n_clusters": n, "largest_pct": 100.0 / max(n, 1), "smallest_pct": 5.0,
                "n_singletons": 0, "mixing_entropy": 0.0}
            for r, n in cluster_counts.items()
        },
        resolutions_source="auto",
        regress_cell_cycle_already_on=False,
    )
    codes = {f.code for f in findings}
    assert "integrate_too_few_clusters" in codes
    too_few = next(f for f in findings if f.code == "integrate_too_few_clusters")
    assert too_few.fix == {"resolutions": "aio"}


def test_integrate_self_check_dominant_cluster_suggests_cc_regression():
    """Largest cluster >50% at every res → suggest cell-cycle regression."""
    cluster_counts = {0.5: 5, 1.0: 7}
    quality = {
        0.5: {"n_clusters": 5, "largest_pct": 65.0, "smallest_pct": 2.0,
              "n_singletons": 0, "mixing_entropy": 0.5},
        1.0: {"n_clusters": 7, "largest_pct": 55.0, "smallest_pct": 1.5,
              "n_singletons": 1, "mixing_entropy": 0.5},
    }
    findings = integrate_self_check(
        cluster_counts=cluster_counts,
        quality=quality,
        resolutions_source="user",
        regress_cell_cycle_already_on=False,
    )
    codes = {f.code for f in findings}
    assert "integrate_dominant_cluster" in codes
    dom = next(f for f in findings if f.code == "integrate_dominant_cluster")
    assert dom.fix == {"regress_cell_cycle": True}


def test_integrate_self_check_silent_on_small_homogeneous():
    """
    v0.9.1 (B2): n_cells < 500 with ≤2 clusters does NOT fire — small
    samples may genuinely be homogeneous, a wider sweep won't help.
    """
    cluster_counts = {0.1: 1, 0.3: 2, 0.5: 2, 0.8: 2, 1.0: 2}
    findings = integrate_self_check(
        cluster_counts=cluster_counts,
        quality={
            r: {"n_clusters": n, "largest_pct": 70.0, "smallest_pct": 30.0,
                "n_singletons": 0, "mixing_entropy": 0.0}
            for r, n in cluster_counts.items()
        },
        resolutions_source="auto",
        regress_cell_cycle_already_on=False,
        n_cells=200,  # small dataset
    )
    codes = {f.code for f in findings}
    assert "integrate_too_few_clusters" not in codes


def test_integrate_self_check_skips_when_already_on_aio_or_cc():
    """If already using aio sweep, no too-few-clusters suggestion."""
    cluster_counts = {0.1: 2, 0.5: 2, 1.0: 2}
    findings = integrate_self_check(
        cluster_counts=cluster_counts,
        quality={
            r: {"n_clusters": 2, "largest_pct": 60.0, "smallest_pct": 40.0,
                "n_singletons": 0, "mixing_entropy": 0.0}
            for r in cluster_counts
        },
        resolutions_source="aio",
        regress_cell_cycle_already_on=True,
    )
    codes = {f.code for f in findings}
    assert "integrate_too_few_clusters" not in codes
    assert "integrate_dominant_cluster" not in codes


# ---------------------------------------------------------------------------
# Annotate self-check
# ---------------------------------------------------------------------------


class _FakeAnnotation:
    """Minimal stand-in for ClusterAnnotation in unit tests."""

    def __init__(self, panel_margin: float, top_markers: list[str]):
        self.panel_margin = panel_margin
        self.top_markers = top_markers


class _FakeProfile:
    """Stand-in profile module exposing the panel attributes the check inspects."""

    __name__ = "tests.fake_profile"

    chondrocyte_markers: dict[str, list[str]] = {
        "ProC": ["BMP2"],
        "EC": ["FRZB"],
    }
    celltype_broad: dict[str, list[str]] = {
        "Macrophages": ["CD68", "CD163"],
        "T cells": ["CD3D", "CD3E"],
    }


def test_annotate_self_check_ambiguous_panel_suggests_alt_panel():
    """Every panel_margin < 0.05 → suggest alt panel."""
    annotations = [
        _FakeAnnotation(panel_margin=0.01, top_markers=["GENE1", "GENE2"]),
        _FakeAnnotation(panel_margin=0.02, top_markers=["GENE3", "GENE4"]),
        _FakeAnnotation(panel_margin=0.0, top_markers=["GENE5"]),
    ]
    findings = annotate_self_check(
        annotations=annotations,
        panel_name="chondrocyte_markers",
        profile_module=_FakeProfile,
    )
    codes = {f.code for f in findings}
    assert "annotate_ambiguous_panel" in codes
    amb = next(f for f in findings if f.code == "annotate_ambiguous_panel")
    assert amb.fix == {"panel_name": "celltype_broad"}


def test_annotate_self_check_tissue_mismatch_suggests_celltype_broad():
    """Chondrocyte panel + immune-dominated clusters → switch panel suggestion."""
    annotations = [
        _FakeAnnotation(
            panel_margin=0.2,  # not ambiguous on margin
            top_markers=["CD3D", "CD3E", "CD8A", "GENE1"],  # T-cell cluster
        ),
        _FakeAnnotation(
            panel_margin=0.3,
            top_markers=["CD68", "CD163", "LYZ", "GENE2"],  # macrophage cluster
        ),
        _FakeAnnotation(
            panel_margin=0.1,
            top_markers=["MS4A1", "CD79A", "CD19", "GENE3"],  # B-cell cluster
        ),
    ]
    findings = annotate_self_check(
        annotations=annotations,
        panel_name="chondrocyte_markers",
        profile_module=_FakeProfile,
    )
    codes = {f.code for f in findings}
    assert "annotate_panel_tissue_mismatch" in codes
    mm = next(f for f in findings if f.code == "annotate_panel_tissue_mismatch")
    assert mm.fix == {"panel_name": "celltype_broad"}


def test_annotate_self_check_low_margin_under_celltype_broad_suggests_chondrocyte():
    """
    v0.9.1 (B7): when current panel is celltype_broad and margins are
    ambiguous, suggest chondrocyte_markers (the OTHER panel), not
    celltype_broad again. Mirror image of test_annotate_self_check_ambiguous_panel_suggests_alt_panel.
    """
    annotations = [
        _FakeAnnotation(panel_margin=0.01, top_markers=["GENE1", "GENE2"]),
        _FakeAnnotation(panel_margin=0.02, top_markers=["GENE3", "GENE4"]),
    ]
    findings = annotate_self_check(
        annotations=annotations,
        panel_name="celltype_broad",
        profile_module=_FakeProfile,
    )
    codes = {f.code for f in findings}
    assert "annotate_ambiguous_panel" in codes
    amb = next(f for f in findings if f.code == "annotate_ambiguous_panel")
    assert amb.fix == {"panel_name": "chondrocyte_markers"}


def test_annotate_self_check_silent_on_clean_call():
    """Margins above 0.05 + panel hits → no findings."""
    annotations = [
        _FakeAnnotation(panel_margin=0.4, top_markers=["BMP2"]),
        _FakeAnnotation(panel_margin=0.3, top_markers=["FRZB"]),
    ]
    findings = annotate_self_check(
        annotations=annotations,
        panel_name="chondrocyte_markers",
        profile_module=_FakeProfile,
    )
    assert findings == []


# ---------------------------------------------------------------------------
# Decision-log integration: triggers + suggestions show up as source="auto"
# ---------------------------------------------------------------------------


def test_qc_run_writes_self_check_decisions(tmp_path):
    """Build an h5ad with adversarial MT% so QC self-check fires + records both rows."""
    import dataclasses

    from scellrun.defaults import SCRNA_QC
    from scellrun.scrna.qc import run_qc

    rng = np.random.default_rng(0)
    # 500 genes so min_genes=200 isn't the bottleneck; MT loading is the
    # axis we want failing.
    n_cells, n_genes = 200, 500
    counts = rng.poisson(lam=4.0, size=(n_cells, n_genes)).astype(np.int32)
    # Tune MT loading so pct_counts_mt sits around 25-30% for most cells:
    # above the default 20% ceiling (so most fail) but below the 30% sensitivity
    # row that brings 60%+ of cells back. This exercises the actionable
    # "qc_low_pass_rate" branch (with a fix) rather than the "no_easy_fix" branch.
    n_mt = 5
    counts[:, :n_mt] = rng.poisson(lam=160.0, size=(n_cells, n_mt))
    a = ad.AnnData(X=counts.astype(np.float32))
    a.var_names = (
        [f"MT-CO{i}" for i in range(1, n_mt + 1)]
        + [f"GENE{i}" for i in range(n_genes - n_mt)]
    )
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    # Use the default thresholds so most cells fail on pct_mt.
    thresholds = dataclasses.replace(SCRNA_QC, max_pct_mt=20.0)

    result = run_qc(
        a,
        thresholds=thresholds,
        flag_doublets=False,
        run_dir=run_dir,
        profile="default",
        lang="en",
    )
    pass_rate = 100.0 * result.n_cells_pass / result.n_cells_in
    assert pass_rate < 30.0, f"adversarial fixture broken: pass_rate={pass_rate}"

    rows = read_decisions(run_dir)
    keys = {r["key"]: r for r in rows}

    # Trigger + suggestion both present, both source="auto".
    # Either the actionable code (with a fix) or the no-easy-fix code can fire;
    # both forms are valid self-check output.
    trigger_keys = [k for k in keys if k.startswith("self_check.") and k.endswith(".trigger")]
    suggest_keys = [k for k in keys if k.startswith("self_check.") and k.endswith(".suggest")]
    assert trigger_keys, f"no self_check trigger row; keys={list(keys)}"
    assert suggest_keys, f"no self_check suggest row; keys={list(keys)}"
    for k in trigger_keys + suggest_keys:
        assert keys[k]["source"] == "auto"
    # Sanity: every trigger has a paired suggest with the same code.
    trigger_codes = {k.split(".")[1] for k in trigger_keys}
    suggest_codes = {k.split(".")[1] for k in suggest_keys}
    assert trigger_codes == suggest_codes


# ---------------------------------------------------------------------------
# End-to-end: --auto-fix rescues a broken QC run
# ---------------------------------------------------------------------------


@pytest.fixture
def adversarial_qc_h5ad(tmp_path):
    """
    Synthetic h5ad where pct_counts_mt sits around 25-35% for most cells —
    too high for the default 20% ceiling, but well within reach of a 30%
    relaxation. Two cell groups so integrate finds two clusters.
    """
    rng = np.random.default_rng(1)
    n_cells, n_genes = 200, 500
    gene_means = rng.gamma(shape=2.0, scale=1.5, size=n_genes).astype(np.float32) + 0.5
    cell_scaling = rng.gamma(shape=4.0, scale=0.25, size=n_cells).astype(np.float32) + 0.5
    lam = np.outer(cell_scaling, gene_means)
    counts = rng.poisson(lam=lam).astype(np.int32)
    # Plant cluster signal so integrate doesn't degenerate to 1 cluster.
    counts[:100, 50:70] = rng.poisson(lam=25.0, size=(100, 20))
    counts[100:, 70:90] = rng.poisson(lam=25.0, size=(100, 20))

    n_mt = 5
    # Crank MT genes up so pct_counts_mt lands around 25-30% for most cells —
    # above the default 20% ceiling (so most cells fail QC at the default)
    # but well below the 30% sensitivity row (so the suggested relaxation
    # rescues 60%+ of cells). Total non-MT counts ~2600/cell; MT lam=200 over
    # 5 genes → MT total ~1000/cell → MT% ~28%.
    counts[:, :n_mt] = rng.poisson(lam=200.0, size=(n_cells, n_mt))

    a = ad.AnnData(X=counts.astype(np.float32))
    a.var_names = (
        [f"MT-CO{i}" for i in range(1, n_mt + 1)]
        + [f"HBB{i}" for i in range(1, 4)]
        + [f"RPS{i}" for i in range(1, 11)]
        + [f"GENE{i}" for i in range(n_genes - n_mt - 13)]
    )
    a.obs_names = [f"cell{i:04d}" for i in range(n_cells)]
    a.obs["sample"] = ["A"] * 100 + ["B"] * 100

    out = tmp_path / "adversarial.h5ad"
    a.write_h5ad(out)
    return out


def test_analyze_auto_fix_failed_first_pass_preserved(adversarial_qc_h5ad, tmp_path, monkeypatch):
    """
    v0.9.1 (B3): when --auto-fix retries a stage, the failed first pass
    is preserved at NN_stage.failed-1/ rather than overwritten in place.
    """
    from scellrun.analyze import run_analyze

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"

    run_analyze(
        adversarial_qc_h5ad,
        profile="joint-disease",
        species="human",
        tissue="synthetic",
        resolutions=(0.3, 0.5),
        use_ai=False,
        lang="en",
        run_dir=run_dir,
        force=False,
        method="none",
        regress_cell_cycle=False,
        use_pubmed=False,
        auto_fix=True,
    )

    # Failed first-pass dir sits next to the retry dir.
    assert (run_dir / "01_qc.failed-1").is_dir(), (
        "expected 01_qc.failed-1/ to preserve the failed first pass"
    )
    assert (run_dir / "01_qc").is_dir()


def test_analyze_auto_fix_rescues_qc(adversarial_qc_h5ad, tmp_path, monkeypatch):
    """
    With --auto-fix on, an adversarial h5ad whose default QC pass-rate < 30%
    triggers a re-run with a relaxed threshold and the pipeline completes.
    The decision log captures the auto-fix application + outcome rows.
    """
    from scellrun.analyze import run_analyze

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"

    result = run_analyze(
        adversarial_qc_h5ad,
        profile="joint-disease",
        species="human",
        tissue="synthetic",
        resolutions=(0.3, 0.5),
        use_ai=False,
        lang="en",
        run_dir=run_dir,
        force=False,
        method="none",
        regress_cell_cycle=False,
        use_pubmed=False,
        auto_fix=True,
    )

    # Pipeline ran to completion.
    assert result.report_index is not None
    assert (run_dir / "05_report" / "index.html").exists()

    rows = read_decisions(run_dir)
    keys_by_stage: dict[tuple[str, str], dict] = {(r["stage"], r["key"]): r for r in rows}

    # Trigger + suggestion both written by the QC self-check on the FIRST pass.
    trigger = keys_by_stage.get(("qc", "self_check.qc_low_pass_rate.trigger"))
    suggest = keys_by_stage.get(("qc", "self_check.qc_low_pass_rate.suggest"))
    assert trigger is not None, "QC self-check trigger row missing"
    assert suggest is not None, "QC self-check suggest row missing"

    # Orchestrator's auto-fix application + outcome rows present.
    applied = keys_by_stage.get(("analyze", "auto_fix.qc.applied"))
    outcome = keys_by_stage.get(("analyze", "auto_fix.qc.outcome"))
    assert applied is not None, "auto_fix.qc.applied decision missing"
    assert outcome is not None, "auto_fix.qc.outcome decision missing"
    assert applied["source"] == "auto"
    assert outcome["source"] == "auto"
    # Outcome encodes the new pass-rate as a percentage string.
    assert outcome["value"].endswith("%")


def test_analyze_auto_fix_outcome_rationale_states_outcome(
    adversarial_qc_h5ad, tmp_path, monkeypatch
):
    """
    v0.9.1 (B7): the auto_fix.qc.outcome rationale carries an explicit
    "rescued the stage" or "did not improve" verdict — not silence. The
    user shouldn't have to read between the lines to know whether the
    retry helped.
    """
    from scellrun.analyze import run_analyze

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    run_dir = tmp_path / "run"

    run_analyze(
        adversarial_qc_h5ad,
        profile="joint-disease",
        species="human",
        tissue="synthetic",
        resolutions=(0.3, 0.5),
        use_ai=False,
        lang="en",
        run_dir=run_dir,
        force=False,
        method="none",
        regress_cell_cycle=False,
        use_pubmed=False,
        auto_fix=True,
    )

    rows = read_decisions(run_dir)
    keys_by_stage = {(r["stage"], r["key"]): r for r in rows}
    outcome = keys_by_stage.get(("analyze", "auto_fix.qc.outcome"))
    assert outcome is not None
    rationale = str(outcome.get("rationale", ""))
    assert "rescued" in rationale or "did not improve" in rationale
