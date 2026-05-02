"""
Tests for the v1.3.0 export module.

If WeasyPrint is importable (the [export] extra is installed), run a
real conversion and assert the output starts with %PDF-. If not,
skip with a clear message — the CLI surface still gets exercised
by the missing-dep error test.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from scellrun.export import PDF_MAGIC, ExportError, run_export

WEASYPRINT_AVAILABLE = importlib.util.find_spec("weasyprint") is not None


def _seed_run_dir(tmp_path: Path) -> Path:
    run = tmp_path / "run-export-test"
    (run / "05_report").mkdir(parents=True)
    # Minimal valid HTML; WeasyPrint accepts even a near-empty body.
    (run / "05_report" / "index.html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>scellrun export test</title></head>"
        "<body><h1>scellrun export</h1>"
        "<p>Smoke-test report for the v1.3.0 export module.</p>"
        "</body></html>",
        encoding="utf-8",
    )
    return run


def test_export_errors_on_missing_run_dir(tmp_path: Path) -> None:
    with pytest.raises(ExportError):
        run_export(tmp_path / "no-such-run", fmt="pdf")


def test_export_errors_on_missing_report(tmp_path: Path) -> None:
    run = tmp_path / "run-no-report"
    run.mkdir()
    with pytest.raises(ExportError):
        run_export(run, fmt="pdf")


def test_export_errors_on_unsupported_format(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    with pytest.raises(ExportError):
        run_export(run, fmt="markdown")


@pytest.mark.skipif(
    not WEASYPRINT_AVAILABLE,
    reason="weasyprint not installed (skip; install with `pip install scellrun[export]`)",
)
def test_export_pdf_writes_pdf_magic_header(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    out = run_export(run, fmt="pdf")
    assert out.exists()
    assert out.stat().st_size > 0
    head = out.read_bytes()[:5]
    assert head == PDF_MAGIC, f"expected %PDF- header, got {head!r}"


@pytest.mark.skipif(
    not WEASYPRINT_AVAILABLE,
    reason="weasyprint not installed",
)
def test_export_pdf_custom_out_path(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    custom = tmp_path / "custom" / "report.pdf"
    out = run_export(run, fmt="pdf", out=custom)
    assert out == custom
    assert custom.exists()
    assert custom.read_bytes()[:5] == PDF_MAGIC


@pytest.mark.skipif(
    not WEASYPRINT_AVAILABLE,
    reason="weasyprint not installed",
)
def test_export_pdf_landscape(tmp_path: Path) -> None:
    run = _seed_run_dir(tmp_path)
    out = run_export(run, fmt="pdf", landscape=True)
    assert out.exists()
    assert out.read_bytes()[:5] == PDF_MAGIC


def test_export_dep_missing_message_carries_install_hint() -> None:
    """
    The error message must point the user at the correct install incantation.
    """
    from scellrun.export import ExportDepMissing

    e = ExportDepMissing("scellrun export needs the [export] extra: pip install scellrun[export]")
    msg = str(e)
    assert "[export]" in msg
    assert "pip install scellrun[export]" in msg
