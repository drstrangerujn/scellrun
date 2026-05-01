"""
tumor profile — pan-cancer / tumor microenvironment.

Single coarse celltype_broad panel covering the major TME compartments
(malignant proxy via proliferation, TAMs, CAFs, T/B/NK/plasma cells,
endothelial, neutrophils). No fine-subtype panel: tumor heterogeneity is
disease-specific (ccRCC vs HNSCC vs melanoma cells share little), and
the right move when the user wants malignant-cell substructure is to
score copy-number alterations (inferCNV) or run a disease-specific
panel — both out of scope for this profile.

QC thresholds left at the default profile's values; tumor dissociations
share the stress-prone characteristics that motivate the 20% mt cap.

Reference panel: TISCH 2.0 (Han et al. 2023, NAR) + the canonical major
TME marker review:
- Sun D, Wang J, Han Y, et al. TISCH: a comprehensive web resource
  enabling interactive single-cell transcriptome visualization of
  tumor microenvironment. Nucleic Acids Res. 2021;49(D1):D1420-D1430.
  PMID: 33179754.
"""
from scellrun.defaults import ScrnaQCThresholds

species: str = "human"

scrna_qc = ScrnaQCThresholds(
    species="human",
)
"""Tumor dissociations are stress-prone; default thresholds (mt 20%, hb 5%) are appropriate."""


# Pan-cancer TME panel. Order is roughly: malignant proxy → stroma →
# myeloid → lymphoid → endothelial → granulocyte. Per TISCH 2.0 (Sun
# 2021) and field-standard major-class markers.
celltype_broad: dict[str, list[str]] = {
    "Malignant (proliferating)":         ["MKI67", "TOP2A", "PCNA"],
    "Tumor-associated macrophages":      ["CD68", "CD163", "MRC1"],
    "Cancer-associated fibroblasts":     ["ACTA2", "FAP", "COL1A1"],
    "T cells":                           ["CD3D", "CD8A", "CD4"],
    "Regulatory T cells":                ["FOXP3", "IL2RA", "CTLA4"],
    "B cells":                           ["MS4A1", "CD79A", "CD19"],
    "Plasma cells":                      ["MZB1", "JCHAIN", "IGHG1"],
    "NK cells":                          ["NKG7", "KLRD1", "GNLY"],
    "Endothelial cells":                 ["PECAM1", "VWF", "CDH5"],
    "Neutrophils":                       ["CSF3R", "FCGR3B", "S100A8"],
}
