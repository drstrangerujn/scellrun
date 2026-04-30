"""Tests for scellrun.scrna.convert auto-detection + read pipeline."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

from scellrun.scrna.convert import (
    UnsupportedInputError,
    convert,
    detect_format,
    read_any,
)


def _make_10x_mtx(dir_: Path, n_cells: int = 50, n_genes: int = 30) -> None:
    """Drop a minimal valid 10x mtx output into dir_ (gzipped, like cellranger emits)."""
    import gzip

    rng = np.random.default_rng(0)
    mat = sp.csr_matrix(rng.poisson(lam=1.0, size=(n_genes, n_cells)).astype(np.int32))
    dir_.mkdir(parents=True, exist_ok=True)

    mtx_path = dir_ / "matrix.mtx"
    scipy.io.mmwrite(str(mtx_path), mat)
    with open(mtx_path, "rb") as f_in, gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
        f_out.write(f_in.read())
    mtx_path.unlink()

    barcodes = "\n".join(f"BC{i}-1" for i in range(n_cells)) + "\n"
    with gzip.open(dir_ / "barcodes.tsv.gz", "wt") as f:
        f.write(barcodes)

    features = "\n".join(f"ENSG{i}\tGENE{i}\tGene Expression" for i in range(n_genes)) + "\n"
    with gzip.open(dir_ / "features.tsv.gz", "wt") as f:
        f.write(features)


def test_detect_10x_mtx(tmp_path):
    d = tmp_path / "raw"
    _make_10x_mtx(d)
    assert detect_format(d) == "10x_mtx"


def test_detect_h5ad(tmp_path):
    p = tmp_path / "x.h5ad"
    p.write_bytes(b"fake")
    assert detect_format(p) == "h5ad"


def test_detect_loom(tmp_path):
    p = tmp_path / "x.loom"
    p.write_bytes(b"fake")
    assert detect_format(p) == "loom"


def test_detect_unknown_extension_raises(tmp_path):
    p = tmp_path / "x.weird"
    p.write_bytes(b"fake")
    with pytest.raises(UnsupportedInputError):
        detect_format(p)


def test_detect_dir_without_mtx_raises(tmp_path):
    d = tmp_path / "empty_dir"
    d.mkdir()
    with pytest.raises(UnsupportedInputError):
        detect_format(d)


def test_read_any_10x_mtx(tmp_path):
    d = tmp_path / "raw"
    _make_10x_mtx(d, n_cells=20, n_genes=15)
    a = read_any(d)
    assert a.n_obs == 20
    assert a.n_vars == 15


def test_convert_writes_h5ad_roundtrip(tmp_path):
    d = tmp_path / "raw"
    _make_10x_mtx(d, n_cells=30, n_genes=20)
    out = tmp_path / "converted.h5ad"
    a = convert(d, out)
    assert out.exists()
    reloaded = ad.read_h5ad(out)
    assert reloaded.n_obs == 30
    assert reloaded.n_vars == 20
    assert a.n_obs == reloaded.n_obs
