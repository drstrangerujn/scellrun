import pytest

from scellrun.defaults import ScrnaQCThresholds
from scellrun.profiles import list_profiles, load


def test_list_profiles_includes_defaults():
    profiles = list_profiles()
    assert "default" in profiles
    assert "joint-disease" in profiles
    # v1.2.0: cross-disease profiles ship alongside joint-disease.
    assert "tumor" in profiles
    assert "brain" in profiles
    assert "kidney" in profiles


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
    with pytest.raises(ModuleNotFoundError):
        load("not-a-real-profile")


# v1.2.0: cross-disease profiles ship alongside joint-disease.
# Sanity-check each new profile's threshold dataclass + broad panel
# shape; end-to-end dogfooding is a separate cold-validation pass.
@pytest.mark.parametrize("name", ["tumor", "brain", "kidney"])
def test_v120_cross_disease_profile_shape(name: str) -> None:
    mod = load(name)

    # species constant
    assert getattr(mod, "species", None) == "human"

    # scrna_qc is a ScrnaQCThresholds (subclasses ok — snRNA variant inherits).
    assert hasattr(mod, "scrna_qc")
    assert isinstance(mod.scrna_qc, ScrnaQCThresholds)

    # celltype_broad: dict[str, list[str]] with at least 5 cell types,
    # all gene names uppercase strings.
    assert hasattr(mod, "celltype_broad")
    panel = mod.celltype_broad
    assert isinstance(panel, dict)
    assert len(panel) >= 5
    for label, genes in panel.items():
        assert isinstance(label, str) and label, f"empty/non-string label in {name}"
        assert isinstance(genes, list) and genes, f"empty gene list for {label} in {name}"
        for g in genes:
            assert isinstance(g, str), f"non-string gene {g!r} in {name}/{label}"
            assert g == g.upper(), f"gene {g!r} not uppercase in {name}/{label}"


def test_brain_profile_tightens_mt_to_10() -> None:
    mod = load("brain")
    assert mod.scrna_qc.max_pct_mt == 10.0


def test_kidney_profile_uses_mt_15() -> None:
    mod = load("kidney")
    assert mod.scrna_qc.max_pct_mt == 15.0


def test_tumor_profile_keeps_default_mt() -> None:
    mod = load("tumor")
    assert mod.scrna_qc.max_pct_mt == 20.0
