import pytest

from scellrun.runlayout import (
    STAGE_DIRS,
    StageOutputExists,
    default_run_dir,
    stage_dir,
    write_run_meta,
)


def test_stage_dirs_have_expected_keys():
    assert set(STAGE_DIRS) == {"qc", "integrate", "markers", "annotate", "report"}
    assert STAGE_DIRS["qc"] == "01_qc"


def test_default_run_dir_format():
    p = default_run_dir()
    assert p.parent.name == "scellrun_out"
    assert p.name.startswith("run-")


def test_stage_dir_creates_subdir(tmp_path):
    sd = stage_dir(tmp_path, "qc")
    assert sd == tmp_path / "01_qc"
    assert sd.is_dir()


def test_stage_dir_unknown_stage_raises(tmp_path):
    with pytest.raises(ValueError):
        stage_dir(tmp_path, "not-a-stage")


def test_stage_dir_blocks_overwrite(tmp_path):
    sd = stage_dir(tmp_path, "qc")
    (sd / "report.html").write_text("dummy")
    with pytest.raises(StageOutputExists):
        stage_dir(tmp_path, "qc")  # second call without force


def test_stage_dir_force_overrides(tmp_path):
    sd = stage_dir(tmp_path, "qc")
    (sd / "report.html").write_text("dummy")
    sd2 = stage_dir(tmp_path, "qc", force=True)
    assert sd == sd2


def test_write_run_meta_appends(tmp_path):
    write_run_meta(tmp_path, command="scrna qc", params={"foo": 1})
    write_run_meta(tmp_path, command="scrna qc", params={"foo": 2})
    import json

    meta = json.loads((tmp_path / "00_run.json").read_text())
    assert len(meta["stages"]) == 2
    assert meta["stages"][0]["params"]["foo"] == 1
    assert meta["stages"][1]["params"]["foo"] == 2
