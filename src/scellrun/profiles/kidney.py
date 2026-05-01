"""
kidney profile — biopsy / cortex / medulla scRNA / snRNA.

Tightens mt% to 15 — kidney tubular cells are mt-rich (proximal tubule
in particular sits between brain and joint on baseline mt fraction),
so the default 20% is too loose and the textbook 10% would lose real
tubular cells. KPMP cohort working practice anchors the 15% choice.

Panel covers the major nephron compartments (proximal tubule, distal
tubule, thick ascending limb, collecting duct principal vs
intercalated, podocytes, parietal cells), vascular (endothelial,
mesangial), and broad immune.

References:
- Stewart BJ, Ferdinand JR, Young MD, et al. Spatiotemporal immune
  zonation of the human kidney. Science. 2019;365(6460):1461-1466.
  PMID: 31604275.
- Kidney Precision Medicine Project (KPMP) reference atlas:
  https://www.kpmp.org/  (Lake et al. 2023, Nature 619:585-594,
  PMID: 37468583, "An atlas of healthy and injured cell states and
  niches in the human kidney").
"""
from scellrun.defaults import ScrnaQCThresholds, ScrnaQCThresholdsSN

species: str = "human"

scrna_qc = ScrnaQCThresholds(
    species="human",
    max_pct_mt=15.0,
)
"""Proximal tubule is mt-rich; 10% loses real tubular cells, 20% lets through stressed/dying ones — split the difference at 15%."""

snrna_qc = ScrnaQCThresholdsSN(
    species="human",
)
"""snRNA inherits the default 5% mt cap; nuclei should have ~0% mt by definition."""


# Major nephron + vascular + immune compartments per Stewart 2019 / KPMP.
# SLC12A1 marks both DCT and TAL; we keep them as separate entries with
# UMOD anchoring TAL (the canonical TAL marker per KPMP) and SLC12A1
# alone marking DCT — so on a real dataset the cluster scoring TAL
# higher will pick up UMOD and the DCT-only cluster will pick up
# SLC12A1.
celltype_broad: dict[str, list[str]] = {
    "Proximal tubule":                   ["LRP2", "SLC34A1", "SLC5A12"],
    "Distal tubule":                     ["SLC12A1", "EGF"],
    "Thick ascending limb":              ["UMOD", "SLC12A1", "CASR"],
    "Collecting duct principal":         ["AQP2", "AQP3", "FXYD4"],
    "Collecting duct intercalated":      ["SLC4A1", "SLC26A4", "ATP6V0D2"],
    "Podocytes":                         ["NPHS1", "NPHS2", "PODXL"],
    "Parietal cells":                    ["CD24", "PAX8", "CLDN1"],
    "Endothelial cells":                 ["PECAM1", "EHD3", "CDH5"],
    "Mesangial cells":                   ["PDGFRB", "ITGA8", "GATA3"],
    "Macrophages":                       ["CD68", "PTPRC", "C1QA"],
    "T cells":                           ["CD3D", "PTPRC", "CD8A"],
}
