# scellrun roadmap

The shape of scellrun is set by two reference workflows:

- **AIO** — `aio.R` from the in-house Liu-lab R package
  (`/home/rpackage/R/aio/1.0.3.2/scRNA_scripts/aio.R`, 337 lines, 5 stages).
  A function library that runs from `Read10X` to `FindAllMarkers` end-to-end.
- **Rmd** — `scRNA_study.Rmd` (subchondral bone study, 2,795 lines, 10 sections).
  The narrative end-to-end notebook including SingleR auto-annotation and the
  Fan 2024 chondrocyte subtype panel.

Together they define the spine. scellrun is the Python re-implementation of the
*decisions* — not the *code*. R idioms (`eval(parse(...))`, hard-coded `setwd`,
SCT vs LogNormalize forking) don't carry over; the working *practice* does.

---

## Version map

| version | command                  | scope                                                                         | reference                                            |
| ------- | ------------------------ | ----------------------------------------------------------------------------- | ---------------------------------------------------- |
| **v0.1** ✓ | `scellrun scrna qc`     | Per-cell metrics (mt%, ribo%, hb%, n_genes, doublets), flag-don't-drop policy. HTML report.    | AIO stage 1-2 (`Read10X`, `cs()` subset block)       |
| v0.2    | `scellrun scrna integrate` | Cross-sample merge, optional cell-cycle scoring, normalization, PCA, Harmony/RPCA, multi-resolution clustering sweep. Integration QC report. | AIO stage 3-4 + Rmd § 6-8                           |
| **v0.3** ✓ | `scellrun scrna markers`  | Per-cluster differential markers (Wilcoxon, logfc≥1, pct≥0.25, only-positive) across all resolutions. One CSV per resolution + HTML report with top-N per cluster. | AIO stage 5 (`fam()`) + Rmd § 10 (`FindAllMarkers`) |
| **v0.4** ✓ | `scellrun scrna annotate` | Two-tier annotation: deterministic panel match (overlap-fraction with profile's marker dict) + optional `--ai` LLM second opinion via Anthropic API + optional `--pubmed` evidence column. Reports both calls side-by-side; user makes the final pick. | Rmd § 9-10 |
| **v0.5** ✓ | `scellrun report` + integrate quality scoring | Multi-stage HTML report aggregator (`05_report/index.html`) linking QC + integrate + markers + annotate with the three-tier provenance trail. Integrate also gains per-resolution quality metrics + cluster-size chart + optional `--ai` resolution recommender. Print-to-PDF for publication. | Rmd full structure |
| **v0.6** ✓ | `scellrun analyze` | One-shot pipeline: qc → integrate → markers → annotate → report in a single command. New users do not have to know the stage layout. Stage-specific deep customization still goes through the per-stage commands. | PLAN.md v0.6 |
| **v0.7** ✓ | decision log + AI summary | Every non-trivial choice (profile, mt%, sample-key auto-detect, resolution pick, panel pick, AI calls) appended to `00_decisions.jsonl`; the report's first page renders them grouped by stage with rationale + AI badges. The user reads one screen and answers "why mt% 20?" / "why res 0.5?". | PLAN.md v0.7 |
| **v0.8** ✓ | error self-heal + actionable suggestions | Each stage runs a `self_check()`: QC suggests the cheapest threshold relaxation when pass-rate drops below 30%, integrate suggests `--resolutions aio` when every resolution yields ≤2 clusters or `--regress-cell-cycle` when the largest cluster dominates >50%, annotate suggests a different `--panel` (or `celltype_broad` for immune-dominated data) when every cluster's panel margin is <0.05. New `--auto-fix` flag at `analyze` level applies the suggestion and re-runs that stage once. Findings recorded as `source="auto"` rows in the decision log (trigger + suggestion). | PLAN.md v0.8 |
| v0.9+   | Streamlit web UI + distribution polish | v0.9 Streamlit web UI; v1.0 conda-forge + Docker + frozen API. Beyond v1.0: bulk RNA-seq, metabolomics composite scoring, proteomics integration. | PLAN.md v0.9 / v1.0 |

---

## v0.1 ↔ AIO mapping (shipped)

| AIO knob (R)            | scellrun field (Python)            | Notes                                                     |
| ----------------------- | --------------------------------- | --------------------------------------------------------- |
| `NRMI = 200`            | `min_genes = 200`                 | Same.                                                     |
| `NRMA = 4000`           | `max_genes = 4000`                | Multiplet upper cap.                                      |
| `PM = 20`               | `max_pct_mt = 20.0`               | Joint-tissue-aware. snRNA tightens to 5%.                 |
| `min.cells = 3`         | `min_cells_per_gene = 3`          | Gene-level filter applied before per-cell metrics.        |
| `sp = "H"/"M"`          | `species: "human" \| "mouse"`     | Drives downstream gene-list selection (CC genes, panels). |
| `ht` (optional)         | `max_pct_hb = 5.0` (always on)    | Default tightened in `joint-disease` profile to 2%.       |
| `dx` (decontX)          | (not in v0.1)                     | No stable Python port; revisit when celda matures.        |
| `db = T` (DoubletFinder)| `flag_doublets = True` (scrublet) | Same intent, equivalent tool in Python.                   |
| `subset(...)` cells out | `obs.qc_pass` flag column         | Policy divergence: scellrun flags, never silently drops.  |

---

## v0.2 ↔ AIO mapping (shipped)

| AIO/Rmd knob (R)              | scellrun field (Python)             | Notes                                                |
| ----------------------------- | ----------------------------------- | ---------------------------------------------------- |
| `cc=T` (CellCycleScoring)     | `--regress-cell-cycle`              | Tirosh genes, regress out S - G2M difference.        |
| `redu = "harmony"`            | `--method harmony`                  | Default. Falls back to PCA if no sample key.         |
| `redu = "rpca"/"cca"`         | `--method rpca/cca`                 | Currently no-op + note in report (no stable scanpy port). |
| `r = c(0.01, …, 2.0)`         | `--resolutions aio`                 | 13-step sweep available; default is shorter [0.1, 0.3, 0.5, 0.8, 1.0]. |
| `di = 30` (PCs)               | `--n-pcs 30`                        | Same.                                                |
| `sct = F` (no SCTransform)    | (default)                           | LogNormalize + scale; SCT not supported.             |
| `subset(DF_hi.lo == "Singlet")` | `--drop-qc-fail`                  | Drop cells with `scellrun_qc_pass=False`. Same intent. |

---

## v0.2 plan — `scellrun scrna integrate`

Targets Rmd § 6-8 + AIO stage 4. Subcommand contract:

```
scellrun scrna integrate <h5ad>... [--method harmony|rpca|none] \
    [--regress-cc/--no-regress-cc] \
    [--resolutions 0.1,0.3,0.5,0.8,1.0] \
    [--profile default|joint-disease]
```

Decisions to bake in:
- LogNormalize + scale (no SCT in v0.2; SCT in scanpy ecosystem is awkward,
  add later as opt-in if PI sees missed signal).
- Cell-cycle scoring uses Tirosh genes (scanpy ships the list); subtract
  `S_score - G2M_score` like AIO.
- Harmony default integration (most common in PI workflow).
- Multi-resolution sweep at `[0.1, 0.3, 0.5, 0.8, 1.0]` by default;
  AIO's full 13-resolution sweep available via `--resolutions all`.
- HTML report: per-resolution UMAP grid + cluster-vs-sample balloon
  table (the Rmd panel that always lands in figures).

---

## v0.3-v0.5 plan — annotation + report

- v0.3 produces a per-resolution markers.csv (positive markers,
  logfc>1, pct.exp>0.25) — same defaults Rmd uses.
- v0.4 SingleR via `celltypist` (Python port of the same idea, no R
  dependency); panel-match step computes mean-z and dropout-aware
  score per panel group, picks best label per cluster, surfaces ties
  in the report. Embedded PubMed evidence column: top 5 papers per
  marker for the user's tissue keyword.
- v0.5 ties everything: ingests outputs from QC + integrate + markers
  + annotate, renders a single PDF report with the three-tier
  provenance trail (`feedback_report_rules`).

---

## Pipeline mode (v0.4+)

The R AIO call signature is the canonical one-liner:

```r
AIO("path/to/seurat",
    db = T, NRMI = 200, NRMA = 4500, PM = 20,
    di = 19, procession = c(1,2,3,4,5),
    redu = "cca", sct = F)
```

…runs five stages in order, where `procession` selects which to run.
scellrun's Python equivalent (lands in v0.4):

```bash
scellrun run path/to/data \
    --profile joint-disease \
    --steps qc,integrate,markers,annotate \
    --species human \
    --max-genes 4500 \
    --pca-dims 19 \
    --integrate-method cca \
    --no-flag-doublets
```

### Artifact handoff protocol (live since v0.1.2)

Every stage writes into `<run_dir>/<NN_stage>/`. `--run-dir` is shared
across stages; default is `scellrun_out/run-YYYYMMDD-HHMMSS`. Layout:

```
scellrun_out/run-20260430-094530/
    00_run.json          run metadata: every stage call appends here
    01_qc/
        report.html
        per_cell_metrics.csv
    02_integrate/        (v0.2)
    03_markers/          (v0.3)
    04_annotate/         (v0.4)
    05_report/           (v0.5)
```

This is in place *now* so `scellrun run` is later just an orchestration
shell over the existing single-stage commands. No retrofit needed.

### `--steps` (= AIO `procession`)

- `--steps qc` — single stage, equivalent to the current `scellrun scrna qc`
- `--steps integrate,markers` — assumes upstream `01_qc/` exists in the run-dir
- `--steps qc,integrate,markers,annotate` — full chain

### Parameter passing

- High-frequency knobs surface at top level (`--max-genes`, `--integrate-method`,
  `--resolutions`)
- Stage-specific overrides via dotted form (`--qc.flag-doublets=false`)
- Heavy customization → YAML config: `scellrun run --config pipeline.yaml`
  (added only when there's real demand; ROADMAP placeholder)

---

## Non-goals

- Not a workflow manager. nf-core / Snakemake exist; scellrun is a CLI you
  call from one of those if you want orchestration.
- Not a replacement for scanpy/scrublet/celltypist/decoupler. Calls into
  them with sensible parameters, returns artifacts. Thin shim with thick
  defaults.
- Not an LLM agent. Optional `scellrun ai *` namespace later, but the
  spine is deterministic.
