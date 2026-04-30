"""
scellrun CLI entrypoint.

Top-level command groups:
    scellrun scrna ...     single-cell RNA-seq commands
    scellrun profiles ...  list / inspect profiles
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

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
    version: Optional[bool] = typer.Option(
        None, "--version", "-V", help="Show version and exit.", callback=_version_callback, is_eager=True
    ),
) -> None:
    pass


@scrna_app.command("qc")
def scrna_qc(
    h5ad: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True, help="Input .h5ad file."),
    run_dir: Optional[Path] = typer.Option(None, "--run-dir", help="Run root directory. Default: scellrun_out/run-YYYYMMDD-HHMMSS."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Override output dir for this stage. By default, writes to <run-dir>/01_qc/."),
    assay: str = typer.Option("scrna", "--assay", help="Assay flavor: 'scrna' or 'snrna'.", show_default=True),
    species: str = typer.Option("human", "--species", help="'human' or 'mouse' (drives downstream gene-list selection)."),
    profile: str = typer.Option("default", "--profile", "-p", help="Profile name (see `scellrun profiles list`)."),
    max_genes: Optional[int] = typer.Option(None, "--max-genes", help="Override profile max_genes upper cap."),
    flag_doublets: bool = typer.Option(True, "--flag-doublets/--no-flag-doublets", help="Run scrublet for doublet flagging."),
) -> None:
    """
    Compute single-cell QC metrics, FLAG (don't drop) outlier cells, and emit
    an HTML report + per-cell CSV.

    Output goes to <run-dir>/01_qc/ so the v0.4+ pipeline can pick it up.
    """
    import dataclasses

    import anndata as ad

    from scellrun.profiles import load as load_profile
    from scellrun.runlayout import default_run_dir, stage_dir, write_run_meta
    from scellrun.scrna.qc import run_qc, write_report

    if assay not in ("scrna", "snrna"):
        console.print(f"[red]error:[/red] --assay must be 'scrna' or 'snrna', got {assay!r}")
        raise typer.Exit(2)
    if species not in ("human", "mouse"):
        console.print(f"[red]error:[/red] --species must be 'human' or 'mouse', got {species!r}")
        raise typer.Exit(2)

    try:
        prof = load_profile(profile)
    except ModuleNotFoundError:
        console.print(f"[red]error:[/red] unknown profile {profile!r}. Try `scellrun profiles list`.")
        raise typer.Exit(2)

    base_thresholds = prof.snrna_qc if assay == "snrna" else prof.scrna_qc
    overrides: dict = {"species": species}
    if max_genes is not None:
        overrides["max_genes"] = max_genes
    thresholds = dataclasses.replace(base_thresholds, **overrides)

    if run_dir is None:
        run_dir = default_run_dir()
    out_dir = out if out is not None else stage_dir(run_dir, "qc")
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print(
        f"[bold]scellrun scrna qc[/bold]  •  profile=[cyan]{profile}[/cyan]  "
        f"assay=[cyan]{assay}[/cyan]  species=[cyan]{species}[/cyan]"
    )
    console.print(f"run-dir: [dim]{run_dir}[/dim]")
    console.print(f"reading [dim]{h5ad}[/dim]")

    adata = ad.read_h5ad(h5ad)
    console.print(f"loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    result = run_qc(adata, assay=assay, flag_doublets=flag_doublets, thresholds=thresholds)

    pct = 100 * result.n_cells_pass / max(result.n_cells_in, 1)
    console.print(
        f"QC pass: [bold green]{result.n_cells_pass:,}[/bold green] / {result.n_cells_in:,} "
        f"({pct:.1f}%) — pct_mt median {result.pct_mt_median:.2f}%, p95 {result.pct_mt_p95:.2f}%"
    )
    if flag_doublets:
        console.print(f"doublets flagged: {result.n_doublets_flagged:,}")

    html_path = write_report(result, adata, out_dir)
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
        },
    )
    console.print(f"[bold]report:[/bold] [link=file://{html_path.resolve()}]{html_path}[/link]")
    console.print(f"per-cell CSV: {out_dir / 'per_cell_metrics.csv'}")


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
        raise typer.Exit(2)

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
