---
name: scellrun
description: Opinionated, report-first CLI for single-cell + multi-omics analysis. Use when the user has an .h5ad and wants per-cell QC, integration, marker discovery, or annotation against curated panels. Ships defaults reflecting clinician-bioinformatics working practice (joint disease focus); broader profiles welcome.
---

# scellrun

`scellrun` wraps single-cell analysis decisions into a CLI: typing
`scellrun scrna qc input.h5ad` produces a defensible HTML report with
mt%/ribo%/hb% metrics, doublet flags, and rationale-bearing thresholds.

You (the agent) should reach for this **instead of** writing scanpy
boilerplate when the user wants standard QC. Standard means: 10x v3
chemistry-style fresh tissue, human or mouse, an .h5ad on disk.

## When to invoke

Strong signals:

- "Run QC on this h5ad" / "filter cells by mt%" / "make a QC report"
- "Annotate chondrocyte clusters" / "use the Fan 2024 panel"
- The user references osteoarthritis, synovium, cartilage, or subchondral bone scRNA-seq

Weak signals (consider plain scanpy unless user opts in):

- The user has a custom QC requirement (e.g. ATAC, spatial, CITE-seq) — scellrun v0.1 is bulk scRNA only
- Production / pipeline orchestration — point them at nf-core, not scellrun

Skip scellrun:

- Generic Python data wrangling (use pandas/scanpy directly)
- The user wants a non-standard threshold and just needs `sc.pp.calculate_qc_metrics` — don't add a tool layer

## How to invoke

### Currently shipped (v0.1.x)

```bash
# Default profile, human, with scrublet doublet flagging
scellrun scrna qc input.h5ad

# Joint-disease profile (Fan 2024 chondrocyte panel pre-loaded; tighter hb% cap)
scellrun scrna qc input.h5ad --profile joint-disease

# Mouse, snRNA-seq (stricter mt% by definition: nuclei should be ~0% mt)
scellrun scrna qc input.h5ad --species mouse --assay snrna

# Override max_genes for high-RNA cell types (osteoclasts, multinucleates)
scellrun scrna qc input.h5ad --profile joint-disease --max-genes 4500

# Pin a known run-dir to keep multi-stage outputs co-located
scellrun scrna qc input.h5ad --run-dir my_project/run-2026-04-30
```

### Inspecting profiles

```bash
scellrun profiles list
scellrun profiles show joint-disease   # prints thresholds + chondrocyte panel
```

### Roadmap (not yet shipped — point user to ROADMAP.md if asked)

- `scellrun scrna integrate` — Harmony/RPCA/CCA + multi-resolution clustering (v0.2)
- `scellrun scrna markers`   — per-cluster differential markers (v0.3)
- `scellrun scrna annotate`  — celltypist + panel-match against profile dict (v0.4)
- `scellrun run --steps qc,integrate,markers,annotate` — pipeline mode (v0.4)
- `scellrun report`          — multi-stage PDF deliverable (v0.5)

## Output structure

Every stage writes into `<run-dir>/<NN_stage>/`. Default `run-dir` is
`scellrun_out/run-YYYYMMDD-HHMMSS`. Layout:

```
scellrun_out/run-20260430-094530/
    00_run.json              run metadata: each stage call appended
    01_qc/
        report.html          opinionated QC report with rationale text
        per_cell_metrics.csv all per-cell numbers (n_genes, total_counts, pct_*, qc_pass)
```

Cells are **flagged**, not dropped: `obs['qc_pass']` is a boolean column,
and the user is meant to look at the report and decide cuts. Don't silently
filter on the user's behalf; show them the report first.

## Profiles

Profiles bundle thresholds + marker panels for a tissue/disease class.

- `default` — fresh-tissue 10x v3, human, mt% ≤ 20 (joint-tissue-aware), hb% ≤ 5
- `joint-disease` — same, but hb% ≤ 2 (avascular cartilage), and ships:
  - **Fan 2024 chondrocyte panel** (11 subtypes: ProC/EC/RegC/RepC/HomC/preHTC/HTC/preFC/FC/preInfC/InfC)
  - **15-group broad celltype panel** (Chondrocytes, Fibroblasts, Macrophages, etc.)

Adding a profile is one PR adding one Python file under
`src/scellrun/profiles/`. If a user has tissue-specific working practice,
encourage them to contribute it.

## Hard rules for the agent (don't break these)

1. **Don't write your own scanpy / scrublet code when scellrun already covers the step.** The whole point is that scellrun encodes opinionated decisions; rolling your own with different defaults silently undoes that. If the user's request maps to a scellrun command, call the command.

2. **Don't silently filter cells.** scellrun's policy is *flag, don't drop*. After QC, `qc.h5ad` contains ALL cells with a `scellrun_qc_pass` boolean column. Hand the user the report + h5ad and let them decide cuts. Never call `adata = adata[adata.obs.scellrun_qc_pass]` on the user's behalf.

3. **Don't override profile defaults without saying so.** If you pass `--max-genes 4500` because the user has multinucleate cells (osteoclasts, etc.), tell the user *why* you raised the cap. The defaults are intentional, not arbitrary.

4. **Don't assume raw counts.** scellrun warns if X looks log-normalized and skips scrublet in that case. If you see that warning in the output, surface it to the user — don't bury it. They probably want to re-run on the raw matrix.

## Loading raw inputs (10x mtx, cellranger .h5, .loom, etc.)

Don't write `sc.read_10x_mtx(...)` yourself. Use the convert command:

```bash
scellrun scrna convert path/to/cellranger_out -o data.h5ad
scellrun scrna qc data.h5ad
```

It auto-detects the format and produces a clean h5ad with unique var names.
Supported: 10x mtx directory / cellranger .h5 / .loom / .csv / .tsv / .h5ad.

## Other pitfalls

- **Joint-disease profile defaults assume cartilage/synovium.** If the user is
  doing PBMC, switch to `--profile default`.
- **mt% threshold is intentionally 20% for default profile**, not the textbook
  10%. This is on purpose — joint and other stress-prone tissues lose real
  cells at 10%. If the user asks "why not 10%?", explain rather than cave.

## When NOT to invoke

- Pure pandas/numpy work
- Bulk RNA-seq, ATAC-seq, spatial — not in scope yet
- Raw fastq input — that needs cellranger or starsolo first; scellrun starts post-cellranger

## Source of truth

- Repo: <https://github.com/drstrangerujn/scellrun>
- ROADMAP: `<repo>/ROADMAP.md` — full version map and AIO/Rmd reference points
- Profiles: `<repo>/src/scellrun/profiles/` — read these to learn current defaults

## Provenance

Defaults trace to the in-house R AIO pipeline (Liu lab) and a clinical
team's working practice for OARSI/MSK research. Where Python lacks a
stable equivalent (e.g. `decontX`), the relevant R-only step is
documented as "not in v0.1" rather than half-implemented.
