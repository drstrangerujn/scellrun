"""
scellrun CLI entrypoint.

Top-level command groups:
    scellrun scrna ...     single-cell RNA-seq commands
    scellrun profiles ...  list / inspect profiles
"""
from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from scellrun import __version__

app = typer.Typer(
    no_args_is_help=True,
    help="scellrun — opinionated, report-first single-cell + multi-omics CLI.",
    add_completion=False,
)
scrna_app = typer.Typer(no_args_is_help=True, help="Single-cell RNA-seq commands.")
profiles_app = typer.Typer(no_args_is_help=True, help="Inspect available profiles.")
app.add_typer(scrna_app, name="scrna")
app.add_typer(profiles_app, name="profiles")

console = Console()


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"scellrun {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool | None = typer.Option(
        None, "--version", "-V", help="Show version and exit.", callback=_version_callback, is_eager=True
    ),
) -> None:
    pass


@scrna_app.command("convert")
def scrna_convert(
    input_path: Path = typer.Argument(..., exists=True, readable=True, help="Input file or directory (10x mtx dir / .h5 / .loom / .csv / .tsv / .h5ad)."),
    output: Path = typer.Option(..., "--out", "-o", help="Output .h5ad path."),
    fmt: str = typer.Option("auto", "--format", "-f", help="Input format. 'auto' tries to detect from extension/contents."),
) -> None:
    """
    Convert raw single-cell input (10x mtx dir, cellranger .h5, loom, csv, tsv)
    into the .h5ad format that scellrun expects.

    Most users have cellranger output sitting on disk and don't want to write
    Python to load it. This command is the first mile.
    """
    from scellrun.scrna.convert import UnsupportedInputError, convert

    valid = {"auto", "10x_mtx", "10x_h5", "loom", "csv", "tsv", "h5ad"}
    if fmt not in valid:
        console.print(f"[red]error:[/red] --format must be one of {sorted(valid)}")
        raise typer.Exit(2) from None

    console.print(f"[bold]scellrun scrna convert[/bold]  •  format=[cyan]{fmt}[/cyan]")
    console.print(f"reading [dim]{input_path}[/dim]")
    try:
        adata = convert(input_path, output, fmt=fmt)  # type: ignore[arg-type]
    except UnsupportedInputError as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from None
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    console.print(f"wrote: {output}")


@scrna_app.command("qc")
def scrna_qc(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input .h5ad file."),
    run_dir: Path | None = typer.Option(None, "--run-dir", help="Run root directory. Default: scellrun_out/run-YYYYMMDD-HHMMSS."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Override output dir for this stage. By default, writes to <run-dir>/01_qc/."),
    assay: str = typer.Option("scrna", "--assay", help="Assay flavor: 'scrna' or 'snrna'.", show_default=True),
    species: str = typer.Option("human", "--species", help="'human' or 'mouse'. Currently informational; v0.2+ uses it for cell-cycle gene lists."),
    profile: str = typer.Option("default", "--profile", "-p", help="Profile name (see `scellrun profiles list`)."),
    max_genes: int | None = typer.Option(None, "--max-genes", help="Override profile max_genes upper cap."),
    flag_doublets: bool = typer.Option(True, "--flag-doublets/--no-flag-doublets", help="Run scrublet for doublet flagging."),
    write_h5ad: bool = typer.Option(True, "--write-h5ad/--no-write-h5ad", help="Also write annotated qc.h5ad (canonical handoff to v0.2 integrate)."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing artifacts in the stage dir. Default: error if 01_qc/ already has output."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
) -> None:
    """
    Compute single-cell QC metrics, FLAG (don't drop) outlier cells, and emit:
      - report.html
      - per_cell_metrics.csv
      - qc.h5ad  (annotated AnnData, canonical handoff for v0.2 integrate)

    Output goes to <run-dir>/01_qc/ so the v0.4+ pipeline can pick it up.
    """
    import dataclasses

    import anndata as ad

    from scellrun.profiles import load as load_profile
    from scellrun.runlayout import (
        StageOutputExists,
        default_run_dir,
        stage_dir,
        write_run_meta,
    )
    from scellrun.scrna.qc import InvalidInputError, run_qc, write_artifacts

    if assay not in ("scrna", "snrna"):
        console.print(f"[red]error:[/red] --assay must be 'scrna' or 'snrna', got {assay!r}")
        raise typer.Exit(2) from None
    if species not in ("human", "mouse"):
        console.print(f"[red]error:[/red] --species must be 'human' or 'mouse', got {species!r}")
        raise typer.Exit(2) from None
    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None

    try:
        prof = load_profile(profile)
    except ModuleNotFoundError:
        console.print(f"[red]error:[/red] unknown profile {profile!r}. Try `scellrun profiles list`.")
        raise typer.Exit(2) from None

    base_thresholds = prof.snrna_qc if assay == "snrna" else prof.scrna_qc
    overrides: dict = {"species": species}
    if max_genes is not None:
        overrides["max_genes"] = max_genes
    thresholds = dataclasses.replace(base_thresholds, **overrides)

    if run_dir is None:
        run_dir = default_run_dir()
    if out is not None:
        out_dir = out
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        try:
            out_dir = stage_dir(run_dir, "qc", force=force)
        except StageOutputExists as e:
            console.print(f"[red]error:[/red] {e}")
            raise typer.Exit(2) from None

    console.print(
        f"[bold]scellrun scrna qc[/bold]  •  profile=[cyan]{profile}[/cyan]  "
        f"assay=[cyan]{assay}[/cyan]  species=[cyan]{species}[/cyan]"
    )
    console.print(f"run-dir: [dim]{run_dir}[/dim]")
    console.print(f"out-dir: [dim]{out_dir}[/dim]")
    console.print("[bold]effective thresholds:[/bold]")
    for k, v in dataclasses.asdict(thresholds).items():
        console.print(f"  {k:24s} {v}")
    console.print(f"reading [dim]{h5ad}[/dim]")

    adata = ad.read_h5ad(h5ad)
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    user_overrides_for_log = {k: v for k, v in overrides.items() if k != "species"}
    # Detect which CLI flags the caller actually supplied so the decision
    # log can tag them as "user" vs "auto" without false positives. typer
    # doesn't expose this directly; the cheap-and-correct approach is to
    # scan sys.argv for the matching long flags.
    import sys as _sys
    _argv = " ".join(_sys.argv)
    qc_flag_supplied = "--profile" in _argv or "-p" in _sys.argv
    qc_assay_supplied = "--assay" in _argv
    qc_species_supplied = "--species" in _argv
    qc_lang_supplied = "--lang" in _argv
    qc_doublets_supplied = (
        "--flag-doublets" in _argv or "--no-flag-doublets" in _argv
    )
    try:
        result = run_qc(
            adata,
            assay=assay,
            flag_doublets=flag_doublets,
            thresholds=thresholds,
            run_dir=run_dir,
            profile=profile,
            user_thresholds_overrides=user_overrides_for_log,
            profile_applied_thresholds=base_thresholds,
            profile_user_supplied=qc_flag_supplied,
            assay_user_supplied=qc_assay_supplied,
            species_user_supplied=qc_species_supplied,
            lang_user_supplied=qc_lang_supplied,
            flag_doublets_user_supplied=qc_doublets_supplied,
            lang=lang,
        )
    except InvalidInputError as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from None

    if result.raw_counts_check == "looks_like_normalized":
        console.print(
            "[yellow]warning:[/yellow] X looks log-normalized, not raw counts. "
            "pct_counts_* metrics may be misleading; doublet flagging skipped."
        )
    elif result.raw_counts_check == "unknown":
        console.print(
            "[yellow]note:[/yellow] could not confirm whether X is raw counts; "
            "proceeding anyway."
        )

    pct = 100 * result.n_cells_pass / max(result.n_cells_in, 1)
    console.print(
        f"QC pass: [bold green]{result.n_cells_pass:,}[/bold green] / {result.n_cells_in:,} "
        f"({pct:.1f}%) — pct_mt median {result.pct_mt_median:.2f}%, p95 {result.pct_mt_p95:.2f}%"
    )
    if flag_doublets and result.raw_counts_check != "looks_like_normalized":
        console.print(f"doublets flagged: {result.n_doublets_flagged:,}")

    artifacts = write_artifacts(result, adata, out_dir, write_h5ad=write_h5ad, lang=lang)
    write_run_meta(
        run_dir,
        command="scrna qc",
        params={
            "h5ad": h5ad,
            "out_dir": out_dir,
            "assay": assay,
            "species": species,
            "profile": profile,
            "thresholds": thresholds,
            "flag_doublets": flag_doublets,
            "write_h5ad": write_h5ad,
            "force": force,
            "raw_counts_check": result.raw_counts_check,
            "lang": lang,
        },
    )
    console.print(f"[bold]report:[/bold] [link=file://{artifacts['report'].resolve()}]{artifacts['report']}[/link]")
    console.print(f"per-cell CSV: {artifacts['per_cell_metrics']}")
    if "qc_h5ad" in artifacts:
        console.print(f"annotated h5ad: {artifacts['qc_h5ad']}")


@scrna_app.command("integrate")
def scrna_integrate(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input .h5ad (typically <run-dir>/01_qc/qc.h5ad)."),
    run_dir: Path | None = typer.Option(None, "--run-dir", help="Run root directory. Default: same as the qc.h5ad's parent run-dir."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Override output dir for this stage. By default, writes to <run-dir>/02_integrate/."),
    method: str = typer.Option("harmony", "--method", help="Integration method: harmony / rpca / cca / none."),
    sample_key: str | None = typer.Option(None, "--sample-key", help="obs column naming the sample/batch (default: auto-detect orig.ident / sample / batch / donor)."),
    n_pcs: int = typer.Option(30, "--n-pcs", help="Number of PCA components."),
    resolutions: str = typer.Option(
        "0.1,0.3,0.5,0.8,1.0",
        "--resolutions",
        help="Comma-separated Leiden resolutions, or 'aio' for the AIO 13-step sweep.",
    ),
    regress_cell_cycle: bool = typer.Option(False, "--regress-cell-cycle/--no-regress-cell-cycle", help="Score and regress out S/G2M cell-cycle effect."),
    species: str = typer.Option("human", "--species", help="'human' or 'mouse' (drives cell-cycle gene-list selection)."),
    drop_qc_fail: bool = typer.Option(True, "--drop-qc-fail/--keep-qc-fail", help="Drop cells with scellrun_qc_pass=False before integration."),
    write_h5ad: bool = typer.Option(True, "--write-h5ad/--no-write-h5ad", help="Write integrated.h5ad with PCA/UMAP/clusters attached."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing 02_integrate/ artifacts."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
    use_ai: bool = typer.Option(False, "--ai/--no-ai", help="Get an LLM resolution recommendation in the report. Requires ANTHROPIC_API_KEY."),
    ai_model: str = typer.Option("claude-haiku-4-5-20251001", "--ai-model", help="Anthropic model for the resolution recommender."),
    tissue: str | None = typer.Option(None, "--tissue", help="Tissue context (e.g. 'osteoarthritis cartilage'). Used by the AI recommender."),
) -> None:
    """
    Cross-sample integration with multi-resolution Leiden clustering.

    Reads a QC'd h5ad, normalizes + scales, runs PCA, integrates (Harmony
    by default) if a sample key is detected, then sweeps clustering at
    several resolutions. Outputs integrated.h5ad + a UMAP grid + report.

    Maps to AIO stage 4 + Rmd § 6-8.
    """
    import anndata as ad

    from scellrun.runlayout import (
        StageOutputExists,
        default_run_dir,
        stage_dir,
        write_run_meta,
    )
    from scellrun.scrna.integrate import (
        AIO_FULL_RESOLUTIONS,
        DEFAULT_RESOLUTIONS,
        IntegrationError,
        run_integrate,
        write_artifacts,
    )

    if method not in ("harmony", "rpca", "cca", "none"):
        console.print(f"[red]error:[/red] --method must be one of harmony/rpca/cca/none, got {method!r}")
        raise typer.Exit(2) from None
    if species not in ("human", "mouse"):
        console.print(f"[red]error:[/red] --species must be 'human' or 'mouse', got {species!r}")
        raise typer.Exit(2) from None
    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None

    if resolutions == "aio":
        res_tuple = AIO_FULL_RESOLUTIONS
    elif resolutions == "default":
        res_tuple = DEFAULT_RESOLUTIONS
    else:
        try:
            res_tuple = tuple(float(x) for x in resolutions.split(",") if x.strip())
        except ValueError:
            console.print(f"[red]error:[/red] could not parse --resolutions {resolutions!r}")
            raise typer.Exit(2) from None
    if not res_tuple:
        console.print("[red]error:[/red] --resolutions must list at least one value")
        raise typer.Exit(2) from None

    # Default run-dir: assume h5ad lives under <run-dir>/01_qc/, walk up two levels
    if run_dir is None:
        if h5ad.parent.name == "01_qc":
            run_dir = h5ad.parent.parent
        else:
            run_dir = default_run_dir()

    if out is not None:
        out_dir = out
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        try:
            out_dir = stage_dir(run_dir, "integrate", force=force)
        except StageOutputExists as e:
            console.print(f"[red]error:[/red] {e}")
            raise typer.Exit(2) from None

    console.print(
        f"[bold]scellrun scrna integrate[/bold]  •  method=[cyan]{method}[/cyan]  "
        f"species=[cyan]{species}[/cyan]"
    )
    console.print(f"run-dir: [dim]{run_dir}[/dim]")
    console.print(f"out-dir: [dim]{out_dir}[/dim]")
    console.print(f"reading [dim]{h5ad}[/dim]")

    adata = ad.read_h5ad(h5ad)
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if resolutions == "aio":
        resolutions_source = "aio"
    elif resolutions == "default" or resolutions == "0.1,0.3,0.5,0.8,1.0":
        resolutions_source = "auto"
    else:
        resolutions_source = "user"

    import sys as _sys
    _argv = " ".join(_sys.argv)
    integ_method_supplied = "--method" in _argv
    integ_n_pcs_supplied = "--n-pcs" in _argv
    integ_cc_supplied = (
        "--regress-cell-cycle" in _argv or "--no-regress-cell-cycle" in _argv
    )

    try:
        result, integrated = run_integrate(
            adata,
            method=method,
            method_user_supplied=integ_method_supplied,
            sample_key=sample_key,
            n_pcs=n_pcs,
            n_pcs_user_supplied=integ_n_pcs_supplied,
            resolutions=res_tuple,
            regress_cell_cycle=regress_cell_cycle,
            regress_cell_cycle_user_supplied=integ_cc_supplied,
            species=species,
            drop_qc_fail=drop_qc_fail,
            use_ai=use_ai,
            ai_model=ai_model,
            tissue=tissue,
            run_dir=run_dir,
            sample_key_user_supplied=sample_key is not None,
            resolutions_source=resolutions_source,
        )
    except (IntegrationError, NotImplementedError) as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from None

    console.print(f"used [bold green]{result.n_cells_used:,}[/bold green] / {result.n_cells_in:,} cells")
    if result.sample_key:
        console.print(f"sample key: {result.sample_key} ({integrated.obs[result.sample_key].nunique()} levels)")
    else:
        console.print("sample key: [dim](none — single-batch run)[/dim]")
    console.print(f"clusters per resolution: {dict(result.cluster_counts)}")

    artifacts = write_artifacts(result, integrated, out_dir, write_h5ad=write_h5ad, lang=lang)
    write_run_meta(
        run_dir,
        command="scrna integrate",
        params={
            "h5ad": h5ad,
            "out_dir": out_dir,
            "method": method,
            "sample_key": result.sample_key,
            "n_pcs": n_pcs,
            "resolutions": list(res_tuple),
            "regress_cell_cycle": regress_cell_cycle,
            "species": species,
            "drop_qc_fail": drop_qc_fail,
            "force": force,
            "lang": lang,
        },
    )

    console.print(f"[bold]report:[/bold] [link=file://{artifacts['report'].resolve()}]{artifacts['report']}[/link]")
    console.print(f"UMAP grid: {artifacts['umap_grid']}")
    if "integrated_h5ad" in artifacts:
        console.print(f"integrated h5ad: {artifacts['integrated_h5ad']}")


@scrna_app.command("markers")
def scrna_markers(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input integrated .h5ad (typically <run-dir>/02_integrate/integrated.h5ad)."),
    run_dir: Path | None = typer.Option(None, "--run-dir", help="Run root directory. Default: parent of <h5ad>'s 02_integrate/ dir, else fresh."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Override output dir. By default writes to <run-dir>/03_markers/."),
    resolutions: str | None = typer.Option(None, "--resolutions", help="Comma-separated resolutions (e.g. '0.3,0.5'). Default: all leiden_res_* columns present in the h5ad."),
    logfc_threshold: float = typer.Option(1.0, "--logfc-threshold", help="Minimum |log2 fold-change| to keep a marker."),
    pct_min: float = typer.Option(0.25, "--pct-min", help="Minimum fraction of cells in a cluster expressing the gene."),
    only_positive: bool = typer.Option(True, "--only-positive/--all-direction", help="Keep only positive (cluster-up) markers."),
    top_n: int = typer.Option(10, "--top-n", help="How many top markers per cluster to show in the HTML report (CSV always full)."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing 03_markers/ artifacts."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
) -> None:
    """
    Per-cluster differential markers across all (or specified) clustering
    resolutions in an integrated h5ad. Wilcoxon rank-sum, logfc>=1, pct>=0.25,
    positive-only by default — same defaults the in-house pipeline uses.

    Maps to AIO stage 5 (FindAllMarkers) + Rmd § 10 first block.
    """
    import anndata as ad

    from scellrun.runlayout import (
        StageOutputExists,
        default_run_dir,
        stage_dir,
        write_run_meta,
    )
    from scellrun.scrna.markers import run_markers, write_artifacts

    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None

    if resolutions is not None:
        try:
            res_tuple: tuple[float, ...] | None = tuple(
                float(x) for x in resolutions.split(",") if x.strip()
            )
        except ValueError:
            console.print(f"[red]error:[/red] could not parse --resolutions {resolutions!r}")
            raise typer.Exit(2) from None
        if not res_tuple:
            console.print("[red]error:[/red] --resolutions must list at least one value")
            raise typer.Exit(2) from None
    else:
        res_tuple = None

    if run_dir is None:
        if h5ad.parent.name == "02_integrate":
            run_dir = h5ad.parent.parent
        else:
            run_dir = default_run_dir()

    if out is not None:
        out_dir = out
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        try:
            out_dir = stage_dir(run_dir, "markers", force=force)
        except StageOutputExists as e:
            console.print(f"[red]error:[/red] {e}")
            raise typer.Exit(2) from None

    console.print(
        f"[bold]scellrun scrna markers[/bold]  •  "
        f"logfc≥{logfc_threshold}  pct≥{pct_min}  "
        f"only_positive={only_positive}"
    )
    console.print(f"run-dir: [dim]{run_dir}[/dim]")
    console.print(f"out-dir: [dim]{out_dir}[/dim]")
    console.print(f"reading [dim]{h5ad}[/dim]")

    adata = ad.read_h5ad(h5ad)
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    import sys as _sys
    _argv = " ".join(_sys.argv)
    markers_logfc_supplied = "--logfc-threshold" in _argv
    markers_pct_supplied = "--pct-min" in _argv
    markers_only_pos_supplied = (
        "--only-positive" in _argv or "--all-direction" in _argv
    )

    try:
        result, per_res_df = run_markers(
            adata,
            resolutions=res_tuple,
            logfc_threshold=logfc_threshold,
            pct_min=pct_min,
            only_positive=only_positive,
            top_n_per_cluster=top_n,
            run_dir=run_dir,
            logfc_threshold_user_supplied=markers_logfc_supplied,
            pct_min_user_supplied=markers_pct_supplied,
            only_positive_user_supplied=markers_only_pos_supplied,
        )
    except ValueError as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from None

    for res, n in result.n_markers_per_resolution.items():
        console.print(f"  res={res:g}: {n:,} markers across {result.cluster_counts[res]} clusters")

    artifacts = write_artifacts(result, per_res_df, out_dir, lang=lang)
    write_run_meta(
        run_dir,
        command="scrna markers",
        params={
            "h5ad": h5ad,
            "out_dir": out_dir,
            "resolutions": list(result.resolutions),
            "logfc_threshold": logfc_threshold,
            "pct_min": pct_min,
            "only_positive": only_positive,
            "top_n": top_n,
            "force": force,
            "lang": lang,
        },
    )

    console.print(f"[bold]report:[/bold] [link=file://{artifacts['report'].resolve()}]{artifacts['report']}[/link]")
    for key, p in artifacts.items():
        if key.startswith("markers_res_"):
            console.print(f"  {p}")


@scrna_app.command("annotate")
def scrna_annotate(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input integrated .h5ad (typically <run-dir>/02_integrate/integrated.h5ad)."),
    run_dir: Path | None = typer.Option(None, "--run-dir", help="Run root directory. Default: parent of h5ad's 02_integrate/."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Override output dir. By default writes to <run-dir>/04_annotate/."),
    profile: str = typer.Option("default", "--profile", "-p", help="Profile name (must define a panel: chondrocyte_markers or celltype_broad)."),
    panel: str | None = typer.Option(None, "--panel", help="Specific panel name in the profile (e.g. chondrocyte_markers). Default: auto-pick."),
    resolution: float = typer.Option(0.5, "--resolution", help="Which leiden_res_<X> column to annotate."),
    use_ai: bool = typer.Option(False, "--ai/--no-ai", help="Get an LLM second opinion alongside the panel match. Requires ANTHROPIC_API_KEY."),
    ai_model: str = typer.Option("claude-haiku-4-5-20251001", "--ai-model", help="Anthropic model for the AI second opinion."),
    use_pubmed: bool = typer.Option(False, "--pubmed/--no-pubmed", help="Fetch top PubMed papers per top marker, embed as evidence in the report."),
    tissue: str | None = typer.Option(None, "--tissue", help="Tissue context (e.g. 'osteoarthritis cartilage'). Used to scope PubMed and AI prompt."),
    top_n: int = typer.Option(30, "--top-n", help="How many top markers per cluster to consider for matching."),
    write_h5ad: bool = typer.Option(True, "--write-h5ad/--no-write-h5ad", help="Write annotated.h5ad with celltype obs columns."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing 04_annotate/ artifacts."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
) -> None:
    """
    Two-tier cluster annotation: deterministic panel match + optional LLM
    second opinion + optional PubMed evidence column.

    Maps to Rmd § 9-10. The panel-match step uses overlap-fraction with the
    profile's marker dict. The AI step (--ai) is opt-in and never overrides
    the deterministic call — both are reported side-by-side so the user
    chooses.
    """
    import anndata as ad

    from scellrun.profiles import load as load_profile
    from scellrun.runlayout import (
        StageOutputExists,
        default_run_dir,
        stage_dir,
        write_run_meta,
    )
    from scellrun.scrna.annotate import run_annotate, write_artifacts

    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None

    try:
        prof = load_profile(profile)
    except ModuleNotFoundError:
        console.print(f"[red]error:[/red] unknown profile {profile!r}. Try `scellrun profiles list`.")
        raise typer.Exit(2) from None

    if run_dir is None:
        if h5ad.parent.name == "02_integrate":
            run_dir = h5ad.parent.parent
        else:
            run_dir = default_run_dir()

    if out is not None:
        out_dir = out
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        try:
            out_dir = stage_dir(run_dir, "annotate", force=force)
        except StageOutputExists as e:
            console.print(f"[red]error:[/red] {e}")
            raise typer.Exit(2) from None

    console.print(
        f"[bold]scellrun scrna annotate[/bold]  •  profile=[cyan]{profile}[/cyan]  "
        f"resolution=[cyan]{resolution}[/cyan]  ai=[cyan]{use_ai}[/cyan]  "
        f"pubmed=[cyan]{use_pubmed}[/cyan]"
    )
    console.print(f"run-dir: [dim]{run_dir}[/dim]")
    console.print(f"out-dir: [dim]{out_dir}[/dim]")
    console.print(f"reading [dim]{h5ad}[/dim]")

    adata = ad.read_h5ad(h5ad)
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    import sys as _sys
    _argv = " ".join(_sys.argv)
    annot_profile_supplied = "--profile" in _argv or "-p" in _sys.argv
    annot_ai_supplied = "--ai" in _argv or "--no-ai" in _argv
    annot_pubmed_supplied = "--pubmed" in _argv or "--no-pubmed" in _argv

    try:
        result = run_annotate(
            adata,
            prof,
            resolution=resolution,
            panel_name=panel,
            use_ai=use_ai,
            ai_model=ai_model,
            use_pubmed=use_pubmed,
            tissue=tissue,
            top_n_markers=top_n,
            run_dir=run_dir,
            profile_user_supplied=annot_profile_supplied,
            resolution_user_supplied=True,  # CLI always exposes --resolution
            use_ai_user_supplied=annot_ai_supplied,
            use_pubmed_user_supplied=annot_pubmed_supplied,
            tissue_user_supplied=tissue is not None,
        )
    except ValueError as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from None

    console.print(f"annotated [bold green]{len(result.annotations)}[/bold green] clusters using panel [cyan]{result.panel_name}[/cyan]")
    for a in result.annotations[:10]:
        ai_str = f" | AI: {a.ai_label}" if a.ai_label else ""
        console.print(f"  cluster {a.cluster:>3}: {a.panel_label} (score {a.panel_score:.2f}){ai_str}")
    if len(result.annotations) > 10:
        console.print(f"  ... ({len(result.annotations) - 10} more)")

    artifacts = write_artifacts(result, adata, out_dir, write_h5ad=write_h5ad, lang=lang)
    write_run_meta(
        run_dir,
        command="scrna annotate",
        params={
            "h5ad": h5ad,
            "out_dir": out_dir,
            "profile": profile,
            "panel_name": result.panel_name,
            "resolution": resolution,
            "use_ai": use_ai,
            "ai_model": ai_model if use_ai else None,
            "use_pubmed": use_pubmed,
            "tissue": tissue,
            "top_n": top_n,
            "force": force,
            "lang": lang,
        },
    )

    console.print(f"[bold]report:[/bold] [link=file://{artifacts['report'].resolve()}]{artifacts['report']}[/link]")
    console.print(f"annotations CSV: {artifacts['annotations']}")
    if "annotated_h5ad" in artifacts:
        console.print(f"annotated h5ad: {artifacts['annotated_h5ad']}")


@app.command("analyze")
def analyze_cmd(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input .h5ad file."),
    profile: str = typer.Option("default", "--profile", "-p", help="Profile name (see `scellrun profiles list`)."),
    species: str = typer.Option("human", "--species", help="'human' or 'mouse'."),
    tissue: str | None = typer.Option(None, "--tissue", help="Tissue context (e.g. 'OA cartilage'). Drives the AI recommender prompt and PubMed scoping."),
    resolutions: str = typer.Option(
        "default",
        "--resolutions",
        help="Comma-separated Leiden resolutions, 'default' (0.1,0.3,0.5,0.8,1.0), or 'aio' for the 13-step sweep.",
    ),
    use_ai: bool | None = typer.Option(None, "--ai/--no-ai", help="Use Anthropic LLM for resolution recommendation + annotation second opinion. Default: True if ANTHROPIC_API_KEY is set, else False."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
    run_dir: Path | None = typer.Option(None, "--run-dir", help="Run root directory. Default: scellrun_out/run-YYYYMMDD-HHMMSS."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing stage subdirs."),
    max_genes: int | None = typer.Option(None, "--max-genes", help="Override profile max_genes upper cap (QC stage)."),
    method: str = typer.Option("harmony", "--method", help="Integration method: harmony / none. (rpca/cca planned, currently raises NotImplementedError.)"),
    regress_cell_cycle: bool = typer.Option(False, "--regress-cell-cycle/--no-regress-cell-cycle", help="Score and regress out S/G2M cell-cycle effect during integrate."),
    use_pubmed: bool = typer.Option(False, "--pubmed/--no-pubmed", help="Annotate stage: fetch top PubMed papers per top marker."),
    auto_fix: bool = typer.Option(
        False,
        "--auto-fix/--no-auto-fix",
        help=(
            "Re-run a stage once (capped) if its self-check finds a fixable problem "
            "(e.g. low QC pass-rate -> relax mt%, <=2 clusters -> wider sweep). Costs: "
            "the failed first-pass artifacts are preserved at NN_stage.failed-1/, and "
            "the retry runs the FULL stage again (integrate/annotate may be slow). "
            "The retry may still fail to rescue the stage. Default: off."
        ),
    ),
) -> None:
    """
    One-shot pipeline: qc → integrate → markers → annotate → report.

    A single command takes a raw .h5ad and produces a populated run-dir
    with one final `index.html` link. New users do not have to learn the
    stage layout.

    Param surface: this command exposes the most common knobs as flat
    options (`--profile`, `--max-genes`, `--method`, `--regress-cell-cycle`,
    `--resolutions`, `--ai`, `--pubmed`, `--lang`, `--run-dir`, `--force`).
    Stage-specific deep customization (e.g. setting `--qc.flag-doublets=false`
    or `--integrate.n-pcs 50`) is intentionally NOT plumbed through here:
    if you need that level of control you should be running the per-stage
    commands directly. The same defaults that ship with the per-stage
    commands apply here.
    """
    import os
    import sys

    from scellrun.analyze import StageFailure, first_run_hint, run_analyze
    from scellrun.scrna.integrate import AIO_FULL_RESOLUTIONS, DEFAULT_RESOLUTIONS

    method_user_supplied = any(a == "--method" or a.startswith("--method=") for a in sys.argv)

    hint = first_run_hint()
    if hint:
        console.print(f"[dim]{hint}[/dim]")

    if species not in ("human", "mouse"):
        console.print(f"[red]error:[/red] --species must be 'human' or 'mouse', got {species!r}")
        raise typer.Exit(2) from None
    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None
    if method not in ("harmony", "rpca", "cca", "none"):
        console.print(f"[red]error:[/red] --method must be one of harmony/rpca/cca/none, got {method!r}")
        raise typer.Exit(2) from None

    if resolutions == "aio":
        res_value: tuple[float, ...] | str = AIO_FULL_RESOLUTIONS
    elif resolutions == "default":
        res_value = DEFAULT_RESOLUTIONS
    else:
        try:
            res_value = tuple(float(x) for x in resolutions.split(",") if x.strip())
        except ValueError:
            console.print(f"[red]error:[/red] could not parse --resolutions {resolutions!r}")
            raise typer.Exit(2) from None
        if not res_value:
            console.print("[red]error:[/red] --resolutions must list at least one value")
            raise typer.Exit(2) from None

    # Default --ai: True iff an API key is set in the environment.
    if use_ai is None:
        use_ai = bool(os.environ.get("ANTHROPIC_API_KEY"))

    console.print(
        f"[bold]scellrun analyze[/bold]  •  profile=[cyan]{profile}[/cyan]  "
        f"species=[cyan]{species}[/cyan]  tissue=[cyan]{tissue or '(none)'}[/cyan]  "
        f"ai=[cyan]{use_ai}[/cyan]"
    )

    try:
        result = run_analyze(
            h5ad,
            profile=profile,
            species=species,
            tissue=tissue,
            resolutions=res_value,
            use_ai=use_ai,
            lang=lang,
            run_dir=run_dir,
            force=force,
            max_genes=max_genes,
            method=method,
            method_user_supplied=method_user_supplied,
            regress_cell_cycle=regress_cell_cycle,
            use_pubmed=use_pubmed,
            auto_fix=auto_fix,
            on_progress=lambda s: console.print(s),
        )
    except StageFailure as e:
        console.print(f"[red]error:[/red] pipeline aborted at stage {e.stage!r}: {e.original}")
        raise typer.Exit(1) from None

    if result.report_index is not None:
        console.print(
            f"[bold]index:[/bold] [link=file://{result.report_index.resolve()}]{result.report_index}[/link]"
        )


@app.command("report")
def report_cmd(
    run_dir: Path = typer.Argument(..., exists=True, file_okay=False, readable=True, help="Run directory (e.g. scellrun_out/run-...)."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Override output dir. Default: <run-dir>/05_report/."),
    force: bool = typer.Option(False, "--force/--no-force", help="Overwrite existing 05_report/ artifacts."),
    lang: str = typer.Option("en", "--lang", help="Report language: 'en' or 'zh'."),
) -> None:
    """
    Stitch QC + integrate + markers + annotate into a single HTML report
    (index.html) with the three-tier provenance trail. Open in a browser
    to read; print-to-PDF for a publication-style deliverable.
    """
    from scellrun.report import build_report
    from scellrun.runlayout import StageOutputExists, stage_dir, write_run_meta

    if lang not in ("en", "zh"):
        console.print(f"[red]error:[/red] --lang must be 'en' or 'zh', got {lang!r}")
        raise typer.Exit(2) from None

    if out is not None:
        out_dir = out
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        try:
            out_dir = stage_dir(run_dir, "report", force=force)
        except StageOutputExists as e:
            console.print(f"[red]error:[/red] {e}")
            raise typer.Exit(2) from None

    console.print(f"[bold]scellrun report[/bold]  •  run-dir=[dim]{run_dir}[/dim]  out=[dim]{out_dir}[/dim]")
    artifacts = build_report(run_dir, out_dir, lang=lang)
    write_run_meta(
        run_dir,
        command="report",
        params={"out_dir": out_dir, "lang": lang, "force": force},
    )
    console.print(f"[bold]index:[/bold] [link=file://{artifacts['index'].resolve()}]{artifacts['index']}[/link]")


@profiles_app.command("list")
def profiles_list() -> None:
    """List available profiles."""
    from scellrun.profiles import list_profiles

    for name in list_profiles():
        console.print(f"  • {name}")


@profiles_app.command("show")
def profiles_show(
    name: str = typer.Argument("default", help="Profile name."),
) -> None:
    """Show the thresholds and panels (if any) in a profile."""
    from scellrun.profiles import load as load_profile

    try:
        mod = load_profile(name)
    except ModuleNotFoundError:
        console.print(f"[red]error:[/red] unknown profile {name!r}.")
        raise typer.Exit(2) from None

    console.print(f"[bold]profile:[/bold] {name}")
    for attr in ("scrna_qc", "snrna_qc"):
        if hasattr(mod, attr):
            console.print(f"  [cyan]{attr}[/cyan]: {getattr(mod, attr)}")

    for panel_attr in ("chondrocyte_markers", "celltype_broad"):
        if hasattr(mod, panel_attr):
            panel = getattr(mod, panel_attr)
            console.print(f"\n  [cyan]{panel_attr}[/cyan] ({len(panel)} groups):")
            for label, genes in panel.items():
                console.print(f"    {label:10s} {', '.join(genes)}")


if __name__ == "__main__":
    app(prog_name="scellrun")
