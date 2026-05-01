"""
brain profile — cortical / hippocampal scRNA / snRNA.

Tightens mt% to 10 — brain tissue dissociations have lower baseline mt%
than joint or tumor (Allen Brain Atlas / Hodge 2019 typical pct_mt
median <5%, p95 <10%). Anything >10% in fresh brain tissue suggests
dying neurons or technical artifact and is worth flagging.

snRNA tightens to the default snRNA threshold (5%) since nuclei
should have ~0% mt by definition; the brain-specific 10% applies to
fresh-tissue scRNA only.

Panel covers the major cortical/hippocampal classes — excitatory and
inhibitory neuron broad classes (with PV / SST / VIP / LAMP5
interneurons broken out), oligodendrocyte lineage, astrocytes,
microglia, vasculature.

References:
- Tasic B, Yao Z, Graybuck LT, et al. Shared and distinct
  transcriptomic cell types across neocortical areas. Nature.
  2018;563(7729):72-78. PMID: 30382198.
- Hodge RD, Bakken TE, Miller JA, et al. Conserved cell types with
  divergent features in human versus mouse cortex. Nature.
  2019;573(7772):61-68. PMID: 31435019.
"""
from scellrun.defaults import ScrnaQCThresholds, ScrnaQCThresholdsSN

species: str = "human"

scrna_qc = ScrnaQCThresholds(
    species="human",
    max_pct_mt=10.0,
)
"""Brain tissue dissociations baseline at <5% mt; tighten to 10% to flag dying neurons."""

snrna_qc = ScrnaQCThresholdsSN(
    species="human",
)
"""snRNA inherits the default 5% mt cap; nuclei should have ~0% mt by definition."""


# Major cortical / hippocampal classes per Tasic 2018 / Hodge 2019.
# Inhibitory neurons are split into broad GABAergic + the four canonical
# interneuron classes (PVALB, SST, VIP, LAMP5) since those four read
# cleanly on most cortical scRNA / snRNA datasets and answer the most
# common downstream question ("which interneuron subset shifts?").
celltype_broad: dict[str, list[str]] = {
    "Excitatory neurons":                ["SLC17A7", "NEUROD2", "SATB2"],
    "Inhibitory neurons (GABAergic)":    ["GAD1", "GAD2", "SLC32A1"],
    "PV interneurons":                   ["PVALB", "GAD1"],
    "SST interneurons":                  ["SST", "GAD1"],
    "VIP interneurons":                  ["VIP", "GAD1"],
    "LAMP5 interneurons":                ["LAMP5", "GAD1"],
    "Oligodendrocytes":                  ["PLP1", "MOG", "MBP"],
    "OPCs":                              ["PDGFRA", "CSPG4", "OLIG1"],
    "Astrocytes":                        ["GFAP", "AQP4", "S100B"],
    "Microglia":                         ["CSF1R", "P2RY12", "TMEM119"],
    "Endothelial cells":                 ["CLDN5", "FLT1", "PECAM1"],
    "Pericytes":                         ["PDGFRB", "RGS5", "MCAM"],
}
