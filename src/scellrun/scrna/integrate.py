"""
scrna integrate stage: read 01_qc/qc.h5ad, optionally regress out cell-cycle,
normalize + log1p + scale, run PCA, integrate with Harmony (default), then
sweep multi-resolution Leiden clustering and produce a UMAP per resolution.

Mapped from R AIO stage 4 + Rmd § 6-8. Choices:
- LogNormalize + ScaleData (no SCT) — SCT in Python ecosystem is unstable;
  PI's standard pipeline doesn't use it (sct=F in the example AIO call).
- Harmony default; rpca/cca/none also accepted.
- Leiden clustering at resolutions [0.1, 0.3, 0.5, 0.8, 1.0]; AIO's full
  13-resolution sweep available via --resolutions.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

DEFAULT_RESOLUTIONS: tuple[float, ...] = (0.1, 0.3, 0.5, 0.8, 1.0)
AIO_FULL_RESOLUTIONS: tuple[float, ...] = (
    0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
)

# Tirosh 2016 cell-cycle gene lists (human). For mouse, use lowercase first letter.
S_GENES_HUMAN: tuple[str, ...] = (
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
)
G2M_GENES_HUMAN: tuple[str, ...] = (
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
)


@dataclass
class IntegrateResult:
    n_cells_in: int
    n_cells_used: int
    method: str
    n_pcs: int
    resolutions: tuple[float, ...]
    cluster_counts: dict[float, int]
    sample_key: str | None
    cc_regressed: bool


def _to_mouse_case(genes: tuple[str, ...]) -> tuple[str, ...]:
    """Mouse cell-cycle genes are the same symbols, lowercase except first letter."""
    return tuple(g[0] + g[1:].lower() for g in genes)


class IntegrationError(RuntimeError):
    """Raised when an explicitly requested integration method can't be run."""


def run_integrate(
    adata: ad.AnnData,
    *,
    method: str = "harmony",
    sample_key: str | None = None,
    n_pcs: int = 30,
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS,
    regress_cell_cycle: bool = False,
    species: str = "human",
    drop_qc_fail: bool = True,
) -> IntegrateResult:
    """
    Integrate one or more samples and run multi-resolution clustering.

    AnnData expected to have already passed through `scellrun scrna qc`
    (carries `obs.scellrun_qc_pass`). If `drop_qc_fail` is True (default),
    cells failing QC are removed before integration.
    """
    if method in ("rpca", "cca"):
        raise NotImplementedError(
            f"--method {method!r} is not implemented in v0.3 (no stable scanpy port). "
            "Use --method harmony or --method none for now; rpca/cca tracked for v0.4+."
        )
    if method not in ("harmony", "none"):
        raise ValueError(f"unknown method {method!r}; expected harmony / none / rpca / cca.")

    n_cells_in = adata.n_obs

    if drop_qc_fail and "scellrun_qc_pass" in adata.obs.columns:
        adata = adata[adata.obs["scellrun_qc_pass"]].copy()

    # Sample-key auto-detection (Seurat convention: orig.ident first)
    candidates = [c for c in ("orig.ident", "sample", "batch", "donor") if c in adata.obs.columns]
    if sample_key is None and candidates:
        sample_key = candidates[0]
        if len(candidates) > 1:
            print(
                f"[scellrun] note: multiple potential sample keys detected: {candidates}. "
                f"Using {sample_key!r}; pass --sample-key to override."
            )

    n_cells_used = adata.n_obs

    # Cell-cycle scoring (optional regression target). Done BEFORE HVG so the
    # downstream regression target is well-defined.
    cc_regressed = False
    if regress_cell_cycle:
        s_genes = list(S_GENES_HUMAN if species == "human" else _to_mouse_case(S_GENES_HUMAN))
        g2m_genes = list(G2M_GENES_HUMAN if species == "human" else _to_mouse_case(G2M_GENES_HUMAN))
        # keep only genes actually present in adata
        s_genes_present = [g for g in s_genes if g in adata.var_names]
        g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]
        if not s_genes_present or not g2m_genes_present:
            print(
                f"[scellrun] WARNING: --regress-cell-cycle requested but only "
                f"{len(s_genes_present)}/{len(s_genes)} S genes and "
                f"{len(g2m_genes_present)}/{len(g2m_genes)} G2M genes "
                f"matched adata.var_names (species={species!r}). "
                "Cell-cycle regression SKIPPED — set --species correctly or check gene naming."
            )
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
            adata.obs["CC_difference"] = adata.obs["S_score"] - adata.obs["G2M_score"]
            cc_regressed = True

    # Normalize + log1p if not already done for CC scoring
    if not cc_regressed:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Stash the log-normalized state as .raw BEFORE selecting HVGs / regressing
    # / scaling. This is what v0.3 markers (and v0.4 annotate) will read for
    # rank_genes_groups, so log2fc numbers stay in normalized-expression units.
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

    if cc_regressed:
        sc.pp.regress_out(adata, ["CC_difference"])

    # Scale + PCA on HVGs only — matches Seurat's default and avoids paying the
    # HVG-selection cost for nothing.
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata,
        n_comps=min(n_pcs, adata.n_obs - 1, adata.n_vars - 1),
        use_highly_variable=True,
    )

    # Integration
    use_rep = "X_pca"
    if method == "harmony":
        if sample_key is None:
            raise IntegrationError(
                "--method harmony requires a sample/batch key, but none was provided "
                "and none of (orig.ident, sample, batch, donor) is in obs. "
                "Pass --sample-key explicitly or use --method none."
            )
        if adata.obs[sample_key].nunique() <= 1:
            raise IntegrationError(
                f"--method harmony requires multiple levels in obs[{sample_key!r}], "
                f"but only {adata.obs[sample_key].nunique()} found. "
                "Use --method none for a single-batch run."
            )
        try:
            sc.external.pp.harmony_integrate(adata, key=sample_key)
        except (ImportError, AttributeError) as e:
            raise IntegrationError(
                "--method harmony requested but harmonypy isn't installed in the runtime. "
                "Install it (`pip install harmonypy`) or use --method none."
            ) from e
        use_rep = "X_pca_harmony"

    # Neighbors + UMAP + Leiden sweep
    sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=n_pcs)
    sc.tl.umap(adata)

    cluster_counts: dict[float, int] = {}
    for res in resolutions:
        key = f"leiden_res_{res:g}".replace(".", "_")
        # flavor='igraph' + directed=False matches scanpy's stable default for
        # post-1.10 versions; suppresses the FutureWarning and is faster.
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        cluster_counts[res] = int(adata.obs[key].nunique())

    return IntegrateResult(
        n_cells_in=n_cells_in,
        n_cells_used=n_cells_used,
        method=method,
        n_pcs=n_pcs,
        resolutions=tuple(resolutions),
        cluster_counts=cluster_counts,
        sample_key=sample_key,
        cc_regressed=cc_regressed,
    ), adata


def write_artifacts(
    result: IntegrateResult,
    adata: ad.AnnData,
    out_dir: Path,
    *,
    write_h5ad: bool = True,
    lang: str = "en",
) -> dict[str, Path]:
    """Write integrated.h5ad + UMAP grid + a small HTML report."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, Path] = {}

    if write_h5ad:
        h5ad_path = out_dir / "integrated.h5ad"
        adata.write_h5ad(h5ad_path)
        artifacts["integrated_h5ad"] = h5ad_path

    # UMAP per resolution
    n = len(result.resolutions)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), squeeze=False)
    for i, res in enumerate(result.resolutions):
        key = f"leiden_res_{res:g}".replace(".", "_")
        sc.pl.umap(adata, color=key, ax=axes[0, i], show=False, frameon=False, legend_loc=None)
        axes[0, i].set_title(f"res={res:g} ({result.cluster_counts[res]} clusters)")
    fig.tight_layout()
    umap_path = out_dir / "umap_grid.png"
    fig.savefig(umap_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    artifacts["umap_grid"] = umap_path

    # Cluster x sample contingency for the middle resolution
    ratio_path: Path | None = None
    if result.sample_key is not None and len(result.resolutions) > 0:
        mid_res = result.resolutions[len(result.resolutions) // 2]
        mid_key = f"leiden_res_{mid_res:g}".replace(".", "_")
        ratio = pd.crosstab(adata.obs[mid_key], adata.obs[result.sample_key])
        ratio_path = out_dir / "cluster_by_sample.csv"
        ratio.to_csv(ratio_path)
        artifacts["cluster_by_sample"] = ratio_path

    # HTML report
    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_integrate_zh.html.j2" if lang == "zh" else "scrna_integrate.html.j2"
    template = env.get_template(template_name)
    html = template.render(
        result=result,
        cluster_counts=result.cluster_counts,
        umap_grid_filename=umap_path.name,
        ratio_filename=ratio_path.name if ratio_path else None,
    )
    html_path = out_dir / "report.html"
    html_path.write_text(html, encoding="utf-8")
    artifacts["report"] = html_path
    return artifacts
