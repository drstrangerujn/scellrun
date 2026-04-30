"""
scRNA QC stage: compute per-cell metrics, flag (don't drop) outliers,
and emit:
    - report.html              human-readable summary with rationale text
    - per_cell_metrics.csv     per-cell numerics for downstream join
    - qc.h5ad                  AnnData with qc_pass column attached, ready for v0.2 integrate

Mapped from R AIO stage 1-2 (Read10X + qc subset). Filtering policy
diverges intentionally from AIO: AIO `subset()`s cells out at QC time,
scellrun only flags them (`obs.scellrun_qc_pass` boolean) so the user
picks the cuts after looking at the report.

obs/var column conventions
--------------------------
scellrun-specific columns are prefixed `scellrun_` to avoid colliding
with user metadata. scanpy-convention columns produced by
`sc.pp.calculate_qc_metrics` (n_genes_by_counts, total_counts,
pct_counts_*) keep their canonical names.

  obs["scellrun_qc_pass"]                bool, AND of all sub-tests
  obs["scellrun_qc_pass_min_genes"]      bool
  obs["scellrun_qc_pass_max_genes"]      bool
  obs["scellrun_qc_pass_min_counts"]     bool
  obs["scellrun_qc_pass_pct_mt"]         bool
  obs["scellrun_qc_pass_pct_ribo"]       bool
  obs["scellrun_qc_pass_pct_hb"]         bool
  obs["scellrun_doublet_flag"]           bool, optional (scrublet)

  var["mt"], var["ribo"], var["hb"]      scanpy convention; left unprefixed
                                          because calculate_qc_metrics derives
                                          column names from them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from scellrun.decisions import Decision, record_many
from scellrun.defaults import SCRNA_QC, SNRNA_QC, ScrnaQCThresholds
from scellrun.self_check import SelfCheckFinding, qc_self_check, record_findings


@dataclass
class QCResult:
    n_cells_in: int
    n_cells_pass: int
    n_genes_in: int
    n_genes_after_filter: int
    pct_mt_median: float
    pct_mt_p95: float
    pct_hb_median: float
    n_doublets_flagged: int
    metrics: pd.DataFrame
    thresholds: ScrnaQCThresholds
    raw_counts_check: str  # 'looks_like_raw' | 'looks_like_normalized' | 'unknown'
    flag_breakdown: dict[str, int]  # which threshold rejected how many cells
    top_flagged: pd.DataFrame  # up to 20 worst-offender cells with reasons
    sensitivity: dict[str, list[dict]]  # per-knob: [{threshold, n_pass, pct_pass}, ...]
    findings: list[SelfCheckFinding] = field(default_factory=list)
    """v0.8 self-check findings (low pass-rate etc.); empty when all looks fine."""
    doublet_threshold_collapsed: bool = False
    """True when scrublet ran but auto-threshold yielded zero predicted doublets
    despite a non-trivial doublet_score spread (max > 0.3). The QC report
    surfaces this so users don't read "0 doublets" as "data is clean"."""


class InvalidInputError(ValueError):
    """Raised when the input AnnData fails an upfront sanity check."""


def _annotate_gene_groups(adata: ad.AnnData) -> None:
    """
    Mark mt / ribo / hb genes regardless of species.

    Uppercasing gene symbols first means the same pattern hits both
    human (MT-CO1, HBB) and mouse (mt-Co1, Hbb-bs).
    """
    var_names_upper = adata.var_names.str.upper()
    adata.var["mt"] = var_names_upper.str.startswith(("MT-", "MT."))
    adata.var["ribo"] = var_names_upper.str.startswith(("RPS", "RPL"))
    # HB[^P] excludes HBP1 etc., matches HBA/HBB/HBM/HBQ/HBZ
    adata.var["hb"] = var_names_upper.str.match(r"^HB[^P]")


def _check_raw_counts(adata: ad.AnnData) -> str:
    """
    Heuristic raw-counts check. Returns one of:
      'looks_like_raw'         — integer-valued, max > 1, looks like UMI counts
      'looks_like_normalized'  — float-valued with non-integers OR small max
      'unknown'                — empty matrix or otherwise indeterminate

    scrublet expects raw counts; calculate_qc_metrics still works on
    normalized data but the resulting pct_* values become meaningless.
    """
    import scipy.sparse as sp

    if adata.n_obs == 0 or adata.n_vars == 0:
        return "unknown"

    X = adata.X
    rows = min(200, X.shape[0])
    sample = X[:rows]
    if sp.issparse(sample):
        sample = sample.toarray()
    sample = np.asarray(sample)
    if sample.size == 0:
        return "unknown"

    fmax = float(sample.max())
    is_integer_valued = np.allclose(sample, sample.astype(int))
    if is_integer_valued and fmax > 1:
        return "looks_like_raw"
    if fmax < 30 and not is_integer_valued:
        return "looks_like_normalized"
    return "unknown"


def run_qc(
    adata: ad.AnnData,
    *,
    assay: str = "scrna",
    flag_doublets: bool | None = None,
    thresholds: ScrnaQCThresholds | None = None,
    run_dir: Path | None = None,
    profile: str = "default",
    user_thresholds_overrides: dict[str, object] | None = None,
    lang: str = "en",
) -> QCResult:
    """
    Compute QC metrics in place on `adata` and return a summary.

    NOTE: this function FLAGS cells (sets `scellrun_qc_pass` and friends in
    `adata.obs`) but does NOT drop them. Filtering is the user's decision;
    we make the decision legible by writing a CSV + report alongside the
    annotated h5ad.

    When ``run_dir`` is provided, every non-trivial choice (profile, assay,
    species, mt% applied, doublet flag, raw-counts heuristic outcome) is
    appended to ``<run_dir>/00_decisions.jsonl`` for the v0.7 decision log.
    """
    if adata.n_obs == 0:
        raise InvalidInputError("input AnnData has 0 cells")
    if adata.n_vars == 0:
        raise InvalidInputError("input AnnData has 0 genes")

    if thresholds is None:
        thresholds = SNRNA_QC if assay == "snrna" else SCRNA_QC
    if flag_doublets is None:
        flag_doublets = thresholds.flag_doublets
    flag_doublets_user_choice = flag_doublets

    raw_check = _check_raw_counts(adata)
    raw_check_skipped_doublets = False
    if raw_check == "looks_like_normalized" and flag_doublets:
        # scrublet on normalized data is meaningless; skip rather than mislead.
        flag_doublets = False
        raw_check_skipped_doublets = True

    n_cells_in = adata.n_obs
    n_genes_in = adata.n_vars

    # Gene-level filter (AIO min.cells=3): drop genes detected in too few cells
    # before computing per-cell metrics. Reduces noise in pct_counts_* tails.
    if thresholds.min_cells_per_gene > 0:
        sc.pp.filter_genes(adata, min_cells=thresholds.min_cells_per_gene)

    n_genes_after_filter = adata.n_vars

    _annotate_gene_groups(adata)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    obs = adata.obs
    obs["scellrun_qc_pass_min_genes"] = obs["n_genes_by_counts"] >= thresholds.min_genes
    obs["scellrun_qc_pass_max_genes"] = obs["n_genes_by_counts"] <= thresholds.max_genes
    obs["scellrun_qc_pass_min_counts"] = obs["total_counts"] >= thresholds.min_counts
    obs["scellrun_qc_pass_pct_mt"] = obs["pct_counts_mt"] <= thresholds.max_pct_mt
    obs["scellrun_qc_pass_pct_ribo"] = obs["pct_counts_ribo"] <= thresholds.max_pct_ribo
    obs["scellrun_qc_pass_pct_hb"] = obs["pct_counts_hb"] <= thresholds.max_pct_hb
    obs["scellrun_qc_pass"] = (
        obs["scellrun_qc_pass_min_genes"]
        & obs["scellrun_qc_pass_max_genes"]
        & obs["scellrun_qc_pass_min_counts"]
        & obs["scellrun_qc_pass_pct_mt"]
        & obs["scellrun_qc_pass_pct_ribo"]
        & obs["scellrun_qc_pass_pct_hb"]
    )

    n_doublets = 0
    doublet_threshold_collapsed = False
    if flag_doublets:
        # scanpy 1.11 deprecated sc.external.pp.scrublet → sc.pp.scrublet.
        # Prefer the new path; fall back to the external path for older scanpy.
        scrublet_fn = getattr(sc.pp, "scrublet", None) or sc.external.pp.scrublet
        try:
            scrublet_fn(adata)
            obs["scellrun_doublet_flag"] = obs["predicted_doublet"]
            n_doublets = int(obs["scellrun_doublet_flag"].sum())
            # Scrublet's auto-threshold occasionally fails to find a bimodal
            # cut and returns "no doublets" even when the score distribution
            # is clearly heterogeneous. Surface that case so the QC summary
            # doesn't silently say "0 doublets" when the real story is
            # "auto-threshold collapsed". See ISSUES.md #004 from the v0.7
            # OA dogfood (12k cells, max doublet_score 0.52, predicted=0).
            if (
                n_doublets == 0
                and "doublet_score" in obs.columns
                and float(obs["doublet_score"].max()) > 0.3
            ):
                doublet_threshold_collapsed = True
        except Exception:
            obs["scellrun_doublet_flag"] = False

    metrics = obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
            "scellrun_qc_pass",
        ]
    ].describe(percentiles=[0.05, 0.5, 0.95]).T

    # Flag breakdown: how many cells failed each threshold
    flag_breakdown = {
        "min_genes":  int((~obs["scellrun_qc_pass_min_genes"]).sum()),
        "max_genes":  int((~obs["scellrun_qc_pass_max_genes"]).sum()),
        "min_counts": int((~obs["scellrun_qc_pass_min_counts"]).sum()),
        "max_pct_mt": int((~obs["scellrun_qc_pass_pct_mt"]).sum()),
        "max_pct_ribo": int((~obs["scellrun_qc_pass_pct_ribo"]).sum()),
        "max_pct_hb": int((~obs["scellrun_qc_pass_pct_hb"]).sum()),
    }

    # Top flagged: 20 worst offenders with the reason annotated
    failing = obs[~obs["scellrun_qc_pass"]].copy()
    if len(failing) > 0:
        reasons = []
        for _, row in failing.iterrows():
            why = []
            if not row["scellrun_qc_pass_min_genes"]:
                why.append(f"n_genes {int(row['n_genes_by_counts'])} < {thresholds.min_genes}")
            if not row["scellrun_qc_pass_max_genes"]:
                why.append(f"n_genes {int(row['n_genes_by_counts'])} > {thresholds.max_genes}")
            if not row["scellrun_qc_pass_min_counts"]:
                why.append(f"counts {int(row['total_counts'])} < {thresholds.min_counts}")
            if not row["scellrun_qc_pass_pct_mt"]:
                why.append(f"mt {row['pct_counts_mt']:.1f}% > {thresholds.max_pct_mt}%")
            if not row["scellrun_qc_pass_pct_ribo"]:
                why.append(f"ribo {row['pct_counts_ribo']:.1f}% > {thresholds.max_pct_ribo}%")
            if not row["scellrun_qc_pass_pct_hb"]:
                why.append(f"hb {row['pct_counts_hb']:.1f}% > {thresholds.max_pct_hb}%")
            reasons.append("; ".join(why))
        failing["why_flagged"] = reasons
        # rank by "how badly they failed": sum of normalized z-scores on the
        # offending metric. Cheap proxy: high pct_mt + low n_genes are the
        # usual culprits, just sort by pct_mt then n_genes.
        top_flagged = (
            failing[["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb", "why_flagged"]]
            .sort_values(["pct_counts_mt", "n_genes_by_counts"], ascending=[False, True])
            .head(20)
        )
    else:
        top_flagged = pd.DataFrame(
            columns=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb", "why_flagged"]
        )

    # Sensitivity sweep — show how many cells survive at alternative thresholds.
    # The user uses this to pick a threshold; we don't rerun, just count from
    # the existing per-cell metrics.
    n_total = len(obs)
    sensitivity: dict[str, list[dict]] = {
        "max_pct_mt": [
            {"threshold": t, "n_pass": int((obs["pct_counts_mt"] <= t).sum()),
             "pct_pass": float(100.0 * (obs["pct_counts_mt"] <= t).sum() / max(n_total, 1))}
            for t in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
        ],
        "max_genes": [
            {"threshold": t, "n_pass": int((obs["n_genes_by_counts"] <= t).sum()),
             "pct_pass": float(100.0 * (obs["n_genes_by_counts"] <= t).sum() / max(n_total, 1))}
            for t in (3000, 4000, 5000, 6000, 8000)
        ],
        "min_genes": [
            {"threshold": t, "n_pass": int((obs["n_genes_by_counts"] >= t).sum()),
             "pct_pass": float(100.0 * (obs["n_genes_by_counts"] >= t).sum() / max(n_total, 1))}
            for t in (100, 200, 300, 500, 800)
        ],
    }

    if run_dir is not None:
        _record_qc_decisions(
            run_dir=run_dir,
            profile=profile,
            assay=assay,
            thresholds=thresholds,
            base_default=SNRNA_QC if assay == "snrna" else SCRNA_QC,
            user_overrides=user_thresholds_overrides or {},
            flag_doublets=flag_doublets,
            flag_doublets_user_choice=flag_doublets_user_choice,
            raw_check=raw_check,
            raw_check_skipped_doublets=raw_check_skipped_doublets,
            lang=lang,
        )

    n_cells_pass = int(obs["scellrun_qc_pass"].sum())

    # v0.8 self-check: pass-rate guardrail. Findings are recorded to the
    # decision log (source="auto") and exposed on the result so the
    # orchestrator can apply auto-fixes.
    findings = qc_self_check(
        n_cells_in=n_cells_in,
        n_cells_pass=n_cells_pass,
        sensitivity=sensitivity,
        thresholds=thresholds,
    )
    if run_dir is not None:
        record_findings(run_dir, findings)

    return QCResult(
        n_cells_in=n_cells_in,
        n_cells_pass=n_cells_pass,
        n_genes_in=n_genes_in,
        n_genes_after_filter=n_genes_after_filter,
        pct_mt_median=float(np.median(obs["pct_counts_mt"])),
        pct_mt_p95=float(np.percentile(obs["pct_counts_mt"], 95)),
        pct_hb_median=float(np.median(obs["pct_counts_hb"])),
        n_doublets_flagged=n_doublets,
        metrics=metrics,
        thresholds=thresholds,
        raw_counts_check=raw_check,
        flag_breakdown=flag_breakdown,
        top_flagged=top_flagged,
        sensitivity=sensitivity,
        findings=findings,
        doublet_threshold_collapsed=doublet_threshold_collapsed,
    )


def _record_qc_decisions(
    *,
    run_dir: Path,
    profile: str,
    assay: str,
    thresholds: ScrnaQCThresholds,
    base_default: ScrnaQCThresholds,
    user_overrides: dict[str, object],
    flag_doublets: bool,
    flag_doublets_user_choice: bool,
    raw_check: str,
    raw_check_skipped_doublets: bool,
    lang: str,
) -> None:
    """Append all QC-stage decisions to the run-dir's decision log."""
    decisions: list[Decision] = [
        Decision(
            stage="qc",
            key="profile",
            value=profile,
            default="default",
            source="user" if profile != "default" else "auto",
            rationale=(
                f"profile {profile!r} chosen; profiles ship tissue-tuned thresholds "
                "(mt%, hb%) and marker panels."
            ),
        ),
        Decision(
            stage="qc",
            key="assay",
            value=assay,
            default="scrna",
            source="user" if assay != "scrna" else "auto",
            rationale=(
                "snRNA tightens mt% to 5% (nuclei should be ~0% mt by definition)"
                if assay == "snrna"
                else "fresh-tissue scRNA defaults; mt% ceiling 20% from PI cohort"
            ),
        ),
        Decision(
            stage="qc",
            key="species",
            value=thresholds.species,
            default="human",
            source="user" if thresholds.species != "human" else "auto",
            rationale=(
                "drives mt/ribo/hb gene-pattern matching (uppercased; same regex "
                "covers MT-CO1 and mt-Co1)"
            ),
        ),
        Decision(
            stage="qc",
            key="max_pct_mt",
            value=thresholds.max_pct_mt,
            default=base_default.max_pct_mt,
            source="user" if "max_pct_mt" in user_overrides else "auto",
            rationale=(
                f"mt% ceiling {thresholds.max_pct_mt}% — joint tissue is stress-prone, "
                "the textbook 10% silently drops real chondrocytes (PI cohort 2024-2026, "
                "AIO PM=20)"
            ),
        ),
        Decision(
            stage="qc",
            key="min_genes",
            value=thresholds.min_genes,
            default=base_default.min_genes,
            source="user" if "min_genes" in user_overrides else "auto",
            rationale=(
                f"min n_genes {thresholds.min_genes} — Ilicic 2016 / "
                "Luecken & Theis 2019 droplet-vs-cell floor"
            ),
        ),
        Decision(
            stage="qc",
            key="max_genes",
            value=thresholds.max_genes,
            default=base_default.max_genes,
            source="user" if "max_genes" in user_overrides else "auto",
            rationale=(
                f"max n_genes {thresholds.max_genes} — multiplet upper cap on 10x v3 chemistry; "
                "raise via --max-genes for high-RNA cell types (megakaryocytes, hepatocytes)"
            ),
        ),
        Decision(
            stage="qc",
            key="flag_doublets",
            value=flag_doublets,
            default=base_default.flag_doublets,
            source=(
                "auto"
                if raw_check_skipped_doublets
                else ("user" if flag_doublets_user_choice != base_default.flag_doublets else "auto")
            ),
            rationale=(
                "scrublet skipped — input matrix looks log-normalized, doublet scoring is meaningless on it"
                if raw_check_skipped_doublets
                else (
                    "scrublet on; report score distribution, never auto-drop (AIO leaves DF result attached)"
                    if flag_doublets
                    else "scrublet off by user request"
                )
            ),
        ),
        Decision(
            stage="qc",
            key="raw_counts_check",
            value=raw_check,
            default=None,
            source="auto",
            rationale=(
                "input matrix sampled (200 rows): integer-valued + max>1 = looks_like_raw; "
                "non-integer with small max = looks_like_normalized; otherwise unknown"
            ),
        ),
        Decision(
            stage="qc",
            key="lang",
            value=lang,
            default="en",
            source="user" if lang != "en" else "auto",
            rationale="report language; both EN and ZH templates kept in sync",
        ),
    ]
    record_many(run_dir, decisions)


def write_artifacts(
    result: QCResult,
    adata: ad.AnnData,
    out_dir: Path,
    *,
    write_h5ad: bool = True,
    lang: str = "en",
) -> dict[str, Path]:
    """
    Write report.html, per_cell_metrics.csv, and (by default) qc.h5ad.

    Returns a dict of artifact name → path. The h5ad is the canonical
    handoff for v0.2 `scellrun scrna integrate`.
    """
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)

    artifacts: dict[str, Path] = {}

    csv_path = out_dir / "per_cell_metrics.csv"
    adata.obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
            "scellrun_qc_pass",
        ]
    ].to_csv(csv_path)
    artifacts["per_cell_metrics"] = csv_path

    if write_h5ad:
        h5ad_path = out_dir / "qc.h5ad"
        adata.write_h5ad(h5ad_path)
        artifacts["qc_h5ad"] = h5ad_path

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_qc_zh.html.j2" if lang == "zh" else "scrna_qc.html.j2"
    template = env.get_template(template_name)
    html = template.render(
        result=result,
        thresholds=result.thresholds,
        metrics_html=result.metrics.to_html(float_format=lambda x: f"{x:.2f}"),
        n_cells_dropped=result.n_cells_in - result.n_cells_pass,
        raw_counts_check=result.raw_counts_check,
        flag_breakdown=result.flag_breakdown,
        top_flagged_html=result.top_flagged.to_html(float_format=lambda x: f"{x:.2f}"),
        n_top_flagged=len(result.top_flagged),
        sensitivity=result.sensitivity,
    )

    html_path = out_dir / "report.html"
    html_path.write_text(html, encoding="utf-8")
    artifacts["report"] = html_path
    return artifacts
