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
| **v0.9** ✓ | SKILL.md v2: agent skill realignment + agent dialogue example | `skills/scellrun/SKILL.md` rewritten from a v0.1 QC-only doc into agent-facing operational guidance for the full v0.8 surface: decision tree for picking commands, the `analyze` one-shot, reading `00_decisions.jsonl`, surfacing self-check trigger/suggest pairs, profile selection by tissue, env hygiene, an end-to-end agent dialogue example, and explicit anti-patterns. No code change; the CLI is unchanged. | this repo |
| **v0.9.1** ✓ | patch bundle: P1 fixes from v0.7 codex / v0.8 codex / hospital dogfooding | Group A — Decision schema honesty: profile-applied default replaces library default in the report's "default" column; `source` derived from explicit `is_user_override` flags, not value-vs-default diffs; `--force` reruns truncate prior decisions of the same stage; JSONL appends use `fcntl.flock`; new `schema_version` + `attempt_id` + `fix_payload` fields. Group B — Self-check tightening: QC trigger raised to <60% pass-rate; integrate "too few clusters" silenced for n_cells<500; auto-fix retry preserves failed first-pass at `NN_stage.failed-N/`; `--auto-fix` help text clarified for clinician use; `SelfCheckFinding.fix` now persisted as structured `fix_payload` on the suggest row; orchestrator logs every fired finding (`<stage>.skipped_findings`), not just the one it acted on. Group C — Dogfooding maintainer-side: single-sample input auto-downgrades `--method harmony` to `none` instead of dying; `_pick_best_resolution` now prefers fewest singletons / smallest largest_pct on the all-fragmented branch; integrated/annotated h5ads cast to float32 (annotated swaps in `.raw` log-normalized); index.html gets an "At a glance" block with QC pass rate + cluster count + top labels + self-check codes; `analyze` auto-resumes incomplete prior runs without `--force`; first-run install heaviness hint via `~/.cache/scellrun/installed.touch`; auto-pick chondrocyte-vs-broad panel by data top-marker overlap. | PLAN.md v0.9 + ISSUES.md |
| **v1.0** ✓ | agent demo + distribution polish | The v1.0 surface = the v0.9.1 codepath, dogfooded a second time end-to-end on real OA cartilage data by an LLM agent (this repo's `docs/agent-demo.md` is the verbatim transcript). Adds: a Dockerfile that builds a single-stage conda+pip image, a `docs/quickstart.md` clinician walkthrough, a `docs/contributing.md` covering profile / panel / stage additions, and a top-level `pyproject` development status bump from Alpha → Production/Stable. The CLI surface is frozen for the v1.x series; v0.x churn ends here. | PLAN.md v1.0 |
| **v1.0.1** ✓ | SKILL.md sync to v0.9.1 behaviors + verification-status section | `skills/scellrun/SKILL.md` rewritten to match the actual v0.9.1 + v1.0 surface: QC self-check trigger ceiling reflected as 60% (was stale at 30%), decision-log JSONL sample updated with `schema_version` / `attempt_id` / `fix_payload`, new `.failed-N/` retry-artifact section, harmony→none single-sample auto-degrade documented at the analyze entrypoint, joint-disease panel auto-pick documented in profile guidance, and a top-level "Verification status" pin to v1.0.0 anchoring the staleness window to `docs/agent-demo.md`. `tests/test_skill_md.py` extended with assertions for each of these landmarks. No code change. | this repo |
| **v1.0.2** ✓ | README rewrite around standardization / audit story | `README.md` rewritten from the differentiation angle: leads with the "two analysts → two answers" framing, quotes the 20% reanalysis-divergence finding from PMC9122178, contrasts vanilla LLM agent + scanpy with scellrun on decision log + tested defaults + community profiles + self-check, and reframes profiles as community-encoded working practice. Quick-start moves below the differentiation pitch. New `tests/test_readme.py` enforces grep anchors (`00_decisions.jsonl`, `Fan 2024`, the 20% figure, install fenced code block, "Why this exists" / "Who this is for" H2s). No code change. | `docs/differentiation_research.md` |
| **v1.1.0** | cold-validation gap fixes (panel auto-pick tie-break, panel-tissue-mismatch self-check wired, orchestrator-passed flag source semantics) | Three behavioral fixes for the gaps the v1.0 cold-agent validation surfaced on BML_1 (`docs/v1demo/cold_validation.md`): (1) `_autopick_panel_for_data` now requires chondrocyte hits >= 1.5x broad hits to keep `chondrocyte_markers`; tie or smaller margin swaps to `celltype_broad`; rationale string carries the hit counts. (2) `IMMUNE_MARKER_HINTS` broadened with HLA class II / plasma IGs / mast tryptases / B CD37 / dendritic CD1C; per-cluster threshold lowered from >=2 to >=1 hit; cluster-fraction trigger lowered from >50% to >=40%. (3) New `panel_name_user_supplied` parameter on `run_annotate`: orchestrator passes `False` even when injecting an auto-picked or self-check-fix panel name; the per-stage CLI passes `True` only when `--panel` actually appears on argv. Decision log `annotate.panel` row now correctly reads `source="auto"` for orchestrator-passed values. SKILL.md adds a chondrocyte panel label glossary. | `docs/v1demo/cold_validation.md` |
| **v1.1.1** | SKILL.md patch — version-sync frontmatter, When-NOT-to-use section, Environment hygiene expansion | Doc-only patch addressing the PI's 5-point evaluation of the v1.1.0 SKILL.md. (1) Frontmatter gains `min_scellrun_version`, `tested_against_version`, `schema_version` keys; Verification status section bumps to v1.1.0 and tells the agent to compare installed `scellrun --version` against the pin. (2) New "When NOT to use scellrun" H2 covering Seurat workflows, non-h5ad inputs, custom analysis, orchestrator double-processing, rare cell types, and refusal of non-defensible answers. (3) Environment hygiene expanded with subsections for ssh patterns (login-shell conda activation), multi-sample anndata.concat merge recipe, `.failed-N/` retry-artifact handoff hygiene, and the canonical `05_report/index.html` artifact-handoff path including scp-back-to-laptop. No source code change. | this repo |
| post-v1.0 | conda-forge + registry image + multi-omics | conda-forge feedstock (`scellrun-feedstock`) so install isn't pip-only on locked-down clusters; tagged Docker push to ghcr.io (registry hand-off); bulk RNA-seq subcommand (`scellrun bulk *`); metabolomics composite scoring (composite-vs-N-FDRs, per `feedback_composite_score_vs_multiple_testing`); proteomics integration. Each lands in its own minor version. | — |

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
