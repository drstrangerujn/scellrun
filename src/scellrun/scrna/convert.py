"""
scrna convert: turn raw single-cell outputs into a `.h5ad` file scellrun can read.

This is the "first mile" most users skip past — research scRNA-seq outputs
arrive as 10x cellranger directories, .h5 files, .loom files, or expression
matrices. scellrun assumes .h5ad downstream, so we wrap the standard scanpy
readers under one entrypoint with auto-detection.
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal

import anndata as ad

InputFormat = Literal["auto", "10x_mtx", "10x_h5", "loom", "csv", "tsv", "h5ad"]


class UnsupportedInputError(ValueError):
    """Raised when a file/dir doesn't match any reader scellrun knows about."""


def detect_format(path: Path) -> InputFormat:
    """
    Inspect the path and guess what reader to use.

    A 10x mtx output looks like a directory with `barcodes.tsv[.gz]`,
    `features.tsv[.gz]` or `genes.tsv[.gz]`, and `matrix.mtx[.gz]`.
    """
    if path.is_dir():
        names = {f.name for f in path.iterdir()}
        has_barcodes = any(n.startswith("barcodes.tsv") for n in names)
        has_matrix = any(n.startswith("matrix.mtx") for n in names)
        if has_barcodes and has_matrix:
            return "10x_mtx"
        raise UnsupportedInputError(
            f"{path} is a directory but doesn't contain 10x mtx files "
            "(expected barcodes.tsv*, matrix.mtx*, features.tsv*)."
        )
    suffix = path.suffix.lower()
    if suffix == ".h5":
        return "10x_h5"
    if suffix == ".loom":
        return "loom"
    if suffix == ".h5ad":
        return "h5ad"
    if suffix == ".csv":
        return "csv"
    if suffix in (".tsv", ".txt"):
        return "tsv"
    raise UnsupportedInputError(
        f"can't auto-detect format for {path}. "
        "Pass --format explicitly (10x_mtx / 10x_h5 / loom / csv / tsv / h5ad)."
    )


def read_any(path: Path, fmt: InputFormat = "auto") -> ad.AnnData:
    """
    Read a single-cell input under any supported format and return AnnData.
    """
    import scanpy as sc

    if fmt == "auto":
        fmt = detect_format(path)

    if fmt == "10x_mtx":
        return sc.read_10x_mtx(path)
    if fmt == "10x_h5":
        return sc.read_10x_h5(path)
    if fmt == "loom":
        return sc.read_loom(path)
    if fmt == "h5ad":
        return ad.read_h5ad(path)
    if fmt == "csv":
        # scanpy reads with samples as rows by default; most expression matrices
        # are genes × cells, so transpose
        a = sc.read_csv(path)
        return a.T if a.n_obs > a.n_vars else a
    if fmt == "tsv":
        a = sc.read_text(path, delimiter="\t")
        return a.T if a.n_obs > a.n_vars else a
    raise UnsupportedInputError(f"unknown format {fmt!r}")


def convert(input_path: Path, output_path: Path, fmt: InputFormat = "auto") -> ad.AnnData:
    """
    Read `input_path` (any supported format), make var_names unique, and write
    to `output_path` as h5ad. Returns the loaded AnnData.
    """
    adata = read_any(input_path, fmt=fmt)
    adata.var_names_make_unique()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    return adata
