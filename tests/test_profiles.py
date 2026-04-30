from scellrun.profiles import list_profiles, load


def test_list_profiles_includes_defaults():
    profiles = list_profiles()
    assert "default" in profiles
    assert "joint-disease" in profiles


def test_load_default_profile_has_thresholds():
    mod = load("default")
    assert hasattr(mod, "scrna_qc")
    assert hasattr(mod, "snrna_qc")
    assert mod.scrna_qc.max_pct_mt == 20.0
    assert mod.snrna_qc.max_pct_mt == 5.0


def test_load_joint_disease_has_panel():
    mod = load("joint-disease")
    assert hasattr(mod, "chondrocyte_markers")
    assert "ProC" in mod.chondrocyte_markers
    assert "InfC" in mod.chondrocyte_markers
    # Fan 2024 panel: 11 subtypes
    assert len(mod.chondrocyte_markers) == 11


def test_load_unknown_profile_raises():
    import pytest

    with pytest.raises(ModuleNotFoundError):
        load("not-a-real-profile")
