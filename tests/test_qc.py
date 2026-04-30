import anndata as ad
import numpy as np
import pytest

from scellrun.scrna.qc import InvalidInputError, run_qc, write_artifacts


def test_qc_runs_on_synthetic(synthetic_h5ad):
    a = ad.read_h5ad(synthetic_h5ad)
    result = run_qc(a, flag_doublets=False)
    assert result.n_cells_in == 600
    assert result.n_cells_pass <= result.n_cells_in
    assert "scellrun_qc_pass" in a.obs.columns
    assert a.obs["scellrun_qc_pass"].dtype == bool


def test_qc_namespaced_columns(synthetic_h5ad):
    """Our boolean columns are prefixed scellrun_ so they don't shadow user metadata."""
    a = ad.read_h5ad(synthetic_h5ad)
    run_qc(a, flag_doublets=False)
    expected = {
        "scellrun_qc_pass",
        "scellrun_qc_pass_min_genes",
        "scellrun_qc_pass_max_genes",
        "scellrun_qc_pass_min_counts",
        "scellrun_qc_pass_pct_mt",
        "scellrun_qc_pass_pct_ribo",
        "scellrun_qc_pass_pct_hb",
    }
    assert expected.issubset(set(a.obs.columns))


def test_qc_empty_anndata_raises():
    a = ad.AnnData(X=np.zeros((0, 10)))
    with pytest.raises(InvalidInputError):
        run_qc(a, flag_doublets=False)


def test_qc_zero_genes_raises():
    a = ad.AnnData(X=np.zeros((10, 0)))
    with pytest.raises(InvalidInputError):
        run_qc(a, flag_doublets=False)


def test_write_artifacts_produces_three_files(synthetic_h5ad, tmp_path):
    a = ad.read_h5ad(synthetic_h5ad)
    result = run_qc(a, flag_doublets=False)
    out = tmp_path / "out"
    artifacts = write_artifacts(result, a, out)
    assert artifacts["report"].exists()
    assert artifacts["per_cell_metrics"].exists()
    assert artifacts["qc_h5ad"].exists()


def test_write_artifacts_skip_h5ad(synthetic_h5ad, tmp_path):
    a = ad.read_h5ad(synthetic_h5ad)
    result = run_qc(a, flag_doublets=False)
    out = tmp_path / "out"
    artifacts = write_artifacts(result, a, out, write_h5ad=False)
    assert "qc_h5ad" not in artifacts
    assert artifacts["report"].exists()


def test_write_artifacts_zh_lang(synthetic_h5ad, tmp_path):
    """Chinese report renders when lang='zh'."""
    a = ad.read_h5ad(synthetic_h5ad)
    result = run_qc(a, flag_doublets=False)
    out = tmp_path / "out_zh"
    artifacts = write_artifacts(result, a, out, lang="zh")
    html = artifacts["report"].read_text(encoding="utf-8")
    assert "单细胞 QC 报告" in html
    assert "总览" in html


def test_qc_h5ad_roundtrip(synthetic_h5ad, tmp_path):
    """qc.h5ad written by scellrun must reload with the qc_pass column intact —
    this is the contract for v0.2 integrate."""
    a = ad.read_h5ad(synthetic_h5ad)
    result = run_qc(a, flag_doublets=False)
    out = tmp_path / "out"
    write_artifacts(result, a, out)
    reloaded = ad.read_h5ad(out / "qc.h5ad")
    assert "scellrun_qc_pass" in reloaded.obs.columns
    assert int(reloaded.obs["scellrun_qc_pass"].sum()) == result.n_cells_pass
