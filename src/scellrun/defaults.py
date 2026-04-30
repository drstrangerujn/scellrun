"""
Default thresholds with one-line rationale per knob.

Every value here MUST be explainable in one sentence to a reviewer.
If you can't, take it out. These are the opinionated defaults that
make scellrun scellrun.

Defaults trace back to the in-house R AIO pipeline (Liu lab, 2024-2025):
    /home/rpackage/R/aio/1.0.3.2/scRNA_scripts/aio.R
…translated into Python idioms (no eval(parse), no SCT fallback chain;
scanpy/scrublet equivalents instead of DoubletFinder/decontX where
those have stable Python ports).
"""
from dataclasses import dataclass
from typing import Literal

Species = Literal["human", "mouse"]


@dataclass(frozen=True)
class ScrnaQCThresholds:
    """Per-cell filters for fresh-tissue dissociated scRNA-seq (10x v3 chemistry)."""

    species: Species = "human"
    """Drives mt/ribo/hb gene-pattern selection. AIO uses ^MT-/^HB for human, ^mt-/^Hb for mouse."""

    min_genes: int = 200
    """Cells with <200 detected genes are usually droplets/debris (Ilicic 2016, Luecken & Theis 2019; AIO NRMI=200)."""

    max_genes: int = 4000
    """Cells with >4000 genes are likely multiplets in 10x v3 chemistry (AIO NRMA=4000). Adjust upward for high-RNA cell types."""

    min_counts: int = 500
    """Floor for total UMI; below this, expression is too sparse for downstream tests."""

    max_pct_mt: float = 20.0
    """Synovium and cartilage are stress-prone tissues; the textbook 10% silently kills real chondrocytes (PI cohort 2024-2026, AIO PM=20). Flag, don't drop, anything above."""

    max_pct_ribo: float = 50.0
    """High ribo fraction is informative, not always artifactual; flag, don't drop."""

    max_pct_hb: float = 5.0
    """Hemoglobin gene fraction >5% indicates RBC contamination from incomplete lysis (AIO `ht` knob). Off by default in AIO; on by default here, since most modern protocols can hit it."""

    min_cells_per_gene: int = 3
    """Drop genes detected in <3 cells before computing QC (AIO min.cells=3 in CreateSeuratObject). Reduces noise + speeds downstream."""

    flag_doublets: bool = True
    """Run scrublet by default; report score distribution, do NOT auto-filter (AIO leaves DF result attached as DF_hi.lo so user can choose)."""


@dataclass(frozen=True)
class ScrnaQCThresholdsSN(ScrnaQCThresholds):
    """Single-nucleus variant: stricter mt% — nuclei should have ~0% mt by definition."""

    max_pct_mt: float = 5.0
    """>5% mt in snRNA means cytoplasmic carry-over — investigate, don't blindly include."""


SCRNA_QC = ScrnaQCThresholds()
SNRNA_QC = ScrnaQCThresholdsSN()
