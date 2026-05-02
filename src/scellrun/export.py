"""
v1.3.0 Part B: `scellrun export <run-dir> --format pdf`

Renders ``<run-dir>/05_report/index.html`` to a print-ready PDF via
WeasyPrint. Embedded PNGs resolve via relative paths from the source
HTML, so we render in place rather than copying the report tree.

Optional dep: WeasyPrint is shipped under the ``[export]`` extra.
``pip install scellrun[export]``. If WeasyPrint isn't importable,
``run_export`` raises ``ExportDepMissing`` with the install hint.

The ``--format`` flag is plural-ready so future formats (markdown,
epub) can slot in without a CLI surface change.
"""
from __future__ import annotations

from pathlib import Path

from scellrun.runlayout import STAGE_DIRS

PDF_MAGIC = b"%PDF-"


class ExportDepMissing(RuntimeError):
    """WeasyPrint (or another optional export dep) isn't installed."""


class ExportError(RuntimeError):
    """Something broke during the conversion (HTML missing, render failure)."""


def _default_out(run_dir: Path, fmt: str) -> Path:
    sub = run_dir / STAGE_DIRS["report"]
    sub.mkdir(parents=True, exist_ok=True)
    if fmt == "pdf":
        return sub / "index.pdf"
    raise ExportError(f"unsupported --format {fmt!r}")


def run_export(
    run_dir: Path,
    *,
    fmt: str = "pdf",
    out: Path | None = None,
    landscape: bool = False,
) -> Path:
    """
    Render the run's main report into the requested format. Returns the
    path to the written file.

    Currently only ``pdf`` is implemented; the dispatch is structured so
    additional formats can land without breaking callers.
    """
    run_dir = Path(run_dir)
    if not run_dir.exists():
        raise ExportError(f"run-dir does not exist: {run_dir}")
    src_html = run_dir / STAGE_DIRS["report"] / "index.html"
    if not src_html.exists():
        raise ExportError(
            f"no report at {src_html}; run `scellrun report {run_dir}` first"
        )

    if out is None:
        out = _default_out(run_dir, fmt)
    else:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "pdf":
        return _render_pdf(src_html, out, landscape=landscape)
    raise ExportError(f"unsupported --format {fmt!r}")


def _render_pdf(src_html: Path, out_path: Path, *, landscape: bool) -> Path:
    """Run WeasyPrint over the report HTML."""
    try:
        from weasyprint import CSS, HTML  # type: ignore[import-not-found]
    except ImportError as e:
        raise ExportDepMissing(
            "scellrun export needs the [export] extra: pip install scellrun[export]"
        ) from e

    base_url = str(src_html.parent.resolve()) + "/"
    css_extra = []
    if landscape:
        # `@page { size: landscape; }` flips the page orientation; default
        # is portrait so the no-css path stays simple.
        css_extra.append(CSS(string="@page { size: A4 landscape; }"))

    html = HTML(filename=str(src_html), base_url=base_url)
    html.write_pdf(target=str(out_path), stylesheets=css_extra or None)
    return out_path


__all__ = [
    "ExportDepMissing",
    "ExportError",
    "PDF_MAGIC",
    "run_export",
]
