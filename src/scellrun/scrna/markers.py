"""
scrna markers stage: per-cluster differential expression for each requested
resolution, written as one CSV per resolution + an HTML report listing top
markers per cluster.

Mapped from R AIO stage 5 (`fam()` → FindAllMarkers per resolution) +
Rmd § 10 first block. The test (Wilcoxon rank-sum) and the defaults
(`logfc_threshold=1.0`, `pct_min=0.25`, `only_positive=True`) match the
in-house Seurat conventions.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from scellrun.decisions import Decision, record_many

DEFAULT_LOGFC = 1.0
DEFAULT_PCT_MIN = 0.25
DEFAULT_TOP_N_PER_CLUSTER = 10


@dataclass
class MarkersResult:
    n_cells_used: int
    resolutions: tuple[float, ...]
    cluster_counts: dict[float, int]
    n_markers_per_resolution: dict[float, int]
    logfc_threshold: float
    pct_min: float
    only_positive: bool
    top_n: int


def _resolution_keys(adata: ad.AnnData) -> dict[float, str]:
    """Map resolution float → obs column name written by `scrna integrate`."""
    out: dict[float, str] = {}
    for col in adata.obs.columns:
        if not col.startswith("leiden_res_"):
            continue
        # leiden_res_0_5 → 0.5
        try:
            num = col.removeprefix("leiden_res_").replace("_", ".")
            out[float(num)] = col
        except ValueError:
            continue
    return out


def run_markers(
    adata: ad.AnnData,
    *,
    resolutions: tuple[float, ...] | None = None,
    logfc_threshold: float = DEFAULT_LOGFC,
    pct_min: float = DEFAULT_PCT_MIN,
    only_positive: bool = True,
    top_n_per_cluster: int = DEFAULT_TOP_N_PER_CLUSTER,
    run_dir: Path | None = None,
) -> tuple[MarkersResult, dict[float, pd.DataFrame]]:
    """
    Compute per-cluster markers at one or more clustering resolutions.

    `adata` should already have `leiden_res_*` columns from `scrna integrate`
    and a normalized + log-transformed expression matrix in `.X` (which is
    what `scrna integrate` leaves behind by default).

    Returns the MarkersResult summary plus a dict[resolution → markers DataFrame].
    Each DataFrame has columns: cluster, gene, log2fc, pval, pval_adj, pct.1, pct.2.
    """
    res_to_key = _resolution_keys(adata)
    if not res_to_key:
        raise ValueError(
            "No leiden_res_* columns found in adata.obs — "
            "did you forget to run `scellrun scrna integrate` first?"
        )

    resolutions_user_supplied = resolutions is not None
    if resolutions is None:
        resolutions = tuple(sorted(res_to_key))
    else:
        missing = [r for r in resolutions if r not in res_to_key]
        if missing:
            raise ValueError(
                f"Requested resolutions {missing} not present in adata. "
                f"Available: {sorted(res_to_key)}"
            )

    cluster_counts: dict[float, int] = {}
    n_markers_per_resolution: dict[float, int] = {}
    per_res_df: dict[float, pd.DataFrame] = {}

    for res in resolutions:
        key = res_to_key[res]
        groups = adata.obs[key].astype("category")
        cluster_counts[res] = int(groups.nunique())

        # Prefer .raw (log-normalized expression set by `scrna integrate`).
        # If no .raw is present, fall back to .X with a warning — log2fc on
        # scaled X is in scaled units, which is misleading.
        use_raw = adata.raw is not None
        if not use_raw:
            print(
                "[scellrun] WARNING: adata.raw is unset — running on .X. "
                "If this h5ad came from `scellrun scrna integrate` v0.3+, .raw "
                "should be populated. log2fc values may be in scaled units."
            )
        sc.tl.rank_genes_groups(
            adata,
            groupby=key,
            method="wilcoxon",
            use_raw=use_raw,
            pts=True,  # pct.1 / pct.2 columns
            key_added=f"rank_{key}",
        )

        rg = adata.uns[f"rank_{key}"]
        cluster_names = list(rg["names"].dtype.names)

        rows: list[dict] = []
        for c in cluster_names:
            for i in range(len(rg["names"][c])):
                gene = rg["names"][c][i]
                logfc = float(rg["logfoldchanges"][c][i])
                pval_adj = float(rg["pvals_adj"][c][i])
                pval = float(rg["pvals"][c][i])
                # pts is added when pts=True, dict-of-arrays
                pct_in = float(rg["pts"][c][i]) if "pts" in rg else float("nan")
                pct_out = float(rg["pts_rest"][c][i]) if "pts_rest" in rg else float("nan")
                rows.append({
                    "cluster": c,
                    "gene": gene,
                    "log2fc": logfc,
                    "pct_in_cluster": pct_in,
                    "pct_other": pct_out,
                    "pval": pval,
                    "pval_adj": pval_adj,
                })

        df = pd.DataFrame(rows)
        if only_positive:
            df = df[df["log2fc"] > 0]
        df = df[df["log2fc"].abs() >= logfc_threshold]
        df = df[df["pct_in_cluster"] >= pct_min]
        df = df.sort_values(["cluster", "log2fc"], ascending=[True, False]).reset_index(drop=True)

        per_res_df[res] = df
        n_markers_per_resolution[res] = len(df)

    result = MarkersResult(
        n_cells_used=adata.n_obs,
        resolutions=tuple(resolutions),
        cluster_counts=cluster_counts,
        n_markers_per_resolution=n_markers_per_resolution,
        logfc_threshold=logfc_threshold,
        pct_min=pct_min,
        only_positive=only_positive,
        top_n=top_n_per_cluster,
    )

    if run_dir is not None:
        decisions: list[Decision] = [
            Decision(
                stage="markers",
                key="resolutions",
                value=list(resolutions),
                default=None,
                source="user" if resolutions_user_supplied else "auto",
                rationale=(
                    "user-supplied resolutions list"
                    if resolutions_user_supplied
                    else "auto-detected from leiden_res_* columns in the integrated h5ad"
                ),
            ),
            Decision(
                stage="markers",
                key="logfc_threshold",
                value=logfc_threshold,
                default=DEFAULT_LOGFC,
                source="user" if logfc_threshold != DEFAULT_LOGFC else "auto",
                rationale=(
                    f"|log2fc| >= {logfc_threshold} — Seurat in-house default; "
                    "tightens to 1.0 to drop near-zero-effect markers"
                ),
            ),
            Decision(
                stage="markers",
                key="pct_min",
                value=pct_min,
                default=DEFAULT_PCT_MIN,
                source="user" if pct_min != DEFAULT_PCT_MIN else "auto",
                rationale=(
                    f"min fraction expressing >= {pct_min} — drops genes hit by dropout in the focal cluster"
                ),
            ),
            Decision(
                stage="markers",
                key="only_positive",
                value=only_positive,
                default=True,
                source="user" if not only_positive else "auto",
                rationale=(
                    "positive-only markers (cluster-up) — what FindAllMarkers ships by default"
                    if only_positive
                    else "both directions kept; user wants down-regulated markers too"
                ),
            ),
        ]
        record_many(run_dir, decisions)

    return result, per_res_df


def write_artifacts(
    result: MarkersResult,
    per_res_df: dict[float, pd.DataFrame],
    out_dir: Path,
    *,
    lang: str = "en",
) -> dict[str, Path]:
    """
    Write one CSV per resolution + HTML report listing top-N markers per cluster.
    """
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, Path] = {}

    # CSVs
    csv_paths: dict[float, Path] = {}
    for res, df in per_res_df.items():
        p = out_dir / f"markers_res_{res:g}.csv"
        df.to_csv(p, index=False)
        csv_paths[res] = p
        artifacts[f"markers_res_{res:g}"] = p

    # Top-N table per resolution for report
    top_per_res: dict[float, pd.DataFrame] = {}
    for res, df in per_res_df.items():
        top = (
            df.groupby("cluster", group_keys=False, observed=True)
            .head(result.top_n)
            .reset_index(drop=True)
        )
        top_per_res[res] = top

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_markers_zh.html.j2" if lang == "zh" else "scrna_markers.html.j2"
    template = env.get_template(template_name)

    html = template.render(
        result=result,
        per_res_top_html={
            res: top.to_html(float_format=lambda x: f"{x:.3g}", index=False)
            for res, top in top_per_res.items()
        },
        csv_filenames={res: p.name for res, p in csv_paths.items()},
    )
    html_path = out_dir / "report.html"
    html_path.write_text(html, encoding="utf-8")
    artifacts["report"] = html_path
    return artifacts
