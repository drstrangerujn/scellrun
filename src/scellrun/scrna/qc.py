"""
scRNA QC stage: compute per-cell metrics, flag (don't drop) outliers,
emit a small Pandas summary and an HTML deliverable.

Mapped from R AIO stage 1-2 (Read10X + qc subset). Filtering policy
diverges intentionally from AIO: AIO `subset()`s cells out at QC time,
scellrun only flags them via `obs.qc_pass` so the user picks the cuts
after looking at the report.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from scellrun.defaults import SCRNA_QC, SNRNA_QC, ScrnaQCThresholds


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


def run_qc(
    adata: ad.AnnData,
    *,
    assay: str = "scrna",
    flag_doublets: bool | None = None,
    thresholds: ScrnaQCThresholds | None = None,
) -> QCResult:
    """
    Compute QC metrics in place on `adata` and return a summary.

    NOTE: this function FLAGS cells (sets columns in `adata.obs`) but does
    NOT drop them. Filtering is the user's decision; we make the decision
    legible.
    """
    if thresholds is None:
        thresholds = SNRNA_QC if assay == "snrna" else SCRNA_QC
    if flag_doublets is None:
        flag_doublets = thresholds.flag_doublets

    n_cells_in = adata.n_obs
    n_genes_in = adata.n_vars

    # Gene-level filter (AIO min.cells=3): drop genes detected in too few cells
    # before computing per-cell metrics. Keeps QC distributions interpretable.
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
    obs["qc_pass_min_genes"] = obs["n_genes_by_counts"] >= thresholds.min_genes
    obs["qc_pass_max_genes"] = obs["n_genes_by_counts"] <= thresholds.max_genes
    obs["qc_pass_min_counts"] = obs["total_counts"] >= thresholds.min_counts
    obs["qc_pass_pct_mt"] = obs["pct_counts_mt"] <= thresholds.max_pct_mt
    obs["qc_pass_pct_ribo"] = obs["pct_counts_ribo"] <= thresholds.max_pct_ribo
    obs["qc_pass_pct_hb"] = obs["pct_counts_hb"] <= thresholds.max_pct_hb
    obs["qc_pass"] = (
        obs["qc_pass_min_genes"]
        & obs["qc_pass_max_genes"]
        & obs["qc_pass_min_counts"]
        & obs["qc_pass_pct_mt"]
        & obs["qc_pass_pct_ribo"]
        & obs["qc_pass_pct_hb"]
    )

    n_doublets = 0
    if flag_doublets:
        try:
            sc.external.pp.scrublet(adata)
            obs["doublet_flag"] = obs["predicted_doublet"]
            n_doublets = int(obs["doublet_flag"].sum())
        except Exception:
            obs["doublet_flag"] = False

    metrics = obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
            "qc_pass",
        ]
    ].describe(percentiles=[0.05, 0.5, 0.95]).T

    return QCResult(
        n_cells_in=n_cells_in,
        n_cells_pass=int(obs["qc_pass"].sum()),
        n_genes_in=n_genes_in,
        n_genes_after_filter=n_genes_after_filter,
        pct_mt_median=float(np.median(obs["pct_counts_mt"])),
        pct_mt_p95=float(np.percentile(obs["pct_counts_mt"], 95)),
        pct_hb_median=float(np.median(obs["pct_counts_hb"])),
        n_doublets_flagged=n_doublets,
        metrics=metrics,
        thresholds=thresholds,
    )


def write_report(result: QCResult, adata: ad.AnnData, out_dir: Path) -> Path:
    """Render an HTML report and a CSV of per-cell metrics. Returns path to HTML."""
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "per_cell_metrics.csv"
    adata.obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
            "qc_pass",
        ]
    ].to_csv(csv_path)

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("scrna_qc.html.j2")
    html = template.render(
        result=result,
        thresholds=result.thresholds,
        metrics_html=result.metrics.to_html(float_format=lambda x: f"{x:.2f}"),
        n_cells_dropped=result.n_cells_in - result.n_cells_pass,
    )

    html_path = out_dir / "report.html"
    html_path.write_text(html, encoding="utf-8")
    return html_path
