"""
joint-disease profile — articular cartilage, synovium, subchondral bone.

Tightens a few QC thresholds for stress-prone joint tissues, and ships the
chondrocyte subtype marker panel that downstream `scellrun scrna annotate`
will use (v0.4+).

Reference: Fan et al. 2024, eleven chondrocyte subtypes from human OA cartilage
(see `reference_chondrocyte_subtypes` in PI memory; gene panel mirrored in the
in-house Rmd at ~/Gtest/R_data/南方医科大软骨下骨/scRNA_study.Rmd § 10).
"""
from scellrun.defaults import ScrnaQCThresholds, ScrnaQCThresholdsSN

scrna_qc = ScrnaQCThresholds(
    species="human",
    max_pct_mt=20.0,
    max_pct_hb=2.0,
)
"""Cartilage dissociations rarely have RBCs (avascular tissue). Tighten hb to 2% — anything above suggests synovium contamination."""

snrna_qc = ScrnaQCThresholdsSN(
    species="human",
    max_pct_hb=2.0,
)


# Fan 2024 chondrocyte subtype marker panel (human articular cartilage).
# Order is biologically meaningful: progenitor → effector → terminal stages.
chondrocyte_markers: dict[str, list[str]] = {
    "ProC":    ["C11orf96", "BMP2", "HMGA1"],         # proliferating chondrocytes
    "EC":      ["CHRDL2", "FRZB", "CYTL1"],           # effector chondrocytes
    "RegC":    ["CHI3L1", "CHI3L2"],                  # regulatory chondrocytes
    "RepC":    ["CILP2", "CILP", "OGN"],              # reparative chondrocytes
    "HomC":    ["HSPA1B", "HSPA1A", "HSPA6", "DDIT3", "JUN"],  # homeostatic / stress-response
    "preHTC":  ["PRG4", "ABI3BP", "CRTAC1"],          # pre-hypertrophic
    "HTC":    ["SPP1", "IBSP", "COL10A1"],            # hypertrophic
    "preFC":   ["COL27A1", "PLCG2", "WWP2"],          # pre-fibro
    "FC":      ["MMP2", "COL1A1", "COL1A2"],          # fibro
    "preInfC": ["IFI16", "IFI27"],                    # pre-inflammatory
    "InfC":    ["CXCL8", "CD74", "GPR183"],           # inflammatory
}


# Coarse parent labels for two-tier annotation (Rmd § 10 first-pass).
celltype_broad: dict[str, list[str]] = {
    "Chondrocytes":              ["COL2A1", "ACAN", "SOX9"],
    "Fibroblasts":               ["COL1A1", "COL1A2", "DCN"],
    "Endothelial cells":         ["PECAM1", "VWF", "CDH5"],
    "Pericytes":                 ["RGS5", "MCAM", "ACTA2"],
    "Smooth muscle cells":       ["MYH11", "ACTA2", "TAGLN"],
    "Mesenchymal stromal cells": ["NT5E", "ENG", "THY1"],
    "Macrophages":               ["CD68", "CD163", "LYZ"],
    "Monocytes":                 ["CD14", "FCN1", "S100A8"],
    "T cells":                   ["CD3D", "CD3E", "CD8A"],
    "B cells":                   ["CD19", "MS4A1", "CD79A"],
    "Plasma cells":              ["JCHAIN", "MZB1", "XBP1"],
    "NK cells":                  ["NKG7", "GNLY", "KLRD1"],
    "Dendritic cells":           ["CLEC9A", "FCER1A", "CST3"],
    "Osteoblasts":               ["RUNX2", "BGLAP", "SP7"],
    "Osteoclasts":               ["CTSK", "ACP5", "MMP9"],
}
