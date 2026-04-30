# scellrun v0.7 dogfood — issues from real OA data

Test environment: hospital server, conda env `scellrun-dogfood` (python 3.11.15),
`pip install scellrun` from PyPI (installed v0.7.0). Sample: BML_1
(`~/Gtest/R_data/南方医科大软骨下骨/seurat/BML_1`), 12,451 cells × 33,538 genes,
10x v3 mtx output. One-shot `scellrun analyze --tissue "osteoarthritis cartilage"
--profile joint-disease --no-ai`. Successful run dir on hospital:
`/tmp/scellrun_dogfood/scellrun_out/run-20260430-141432`.

## Table of contents

- [#001] analyze --method default = harmony, single-sample dies — Severity P1
- [#002] markers CSV / report top rows are zero-vs-near-zero noise (log2fc=21, pval~1) — Severity P0
- [#003] joint-disease panel auto-pick assigns "ProC" to every score-0 cluster — Severity P0
- [#004] Scrublet auto-threshold collapse → "0 doublets" silently — Severity P1
- [#005] _pick_best_resolution falls through to fragmented resolution when every res has >1 singleton — Severity P1
- [#006] integrated.h5ad / annotated.h5ad ~2 GB each (dense float64 scaled X) — Severity P2
- [#007] index.html lacks any analysis summary (no QC numbers, no annotation table) — Severity P2
- [#008] Empty `<h2>Decision summary</h2>` headers and JSON-Path leak in report — Severity P2
- [#009] Pipeline noise: scanpy `use_highly_variable` deprecated, scrublet API moved, pandas FutureWarning — Severity P2
- [#010] Re-running `analyze` after a mid-pipeline failure dies on StageOutputExists — Severity P2
- [#011] Install: 1m16s, no progress hint about heaviness; no smoke message for `scellrun --version` — Severity P2

---

## [#001] analyze --method default = harmony, single-sample dies

**Stage:** integrate
**Severity:** P1
**Reproducer:**
```
scellrun analyze /tmp/scellrun_dogfood/dogfood.h5ad \
  --tissue "osteoarthritis cartilage" --profile joint-disease --no-ai
```
**What happened:**
```
[1/5] qc: 10,835 / 12,451 cells passed (87.0%)
error: pipeline aborted at stage 'integrate': --method harmony requires a
sample/batch key, but none was provided and none of (orig.ident, sample,
batch, donor) is in obs. Pass --sample-key explicitly or use --method none.
```
**What should happen:** if the user converted a single 10x dir, `analyze` should
either auto-fall-back to `--method none` (and log it as a `source="auto"`
decision) or do the no-sample-key check up front before running QC, so the user
isn't already 1m49s in when it dies.
**Suspected cause:** `analyze.run_analyze` defaults `method="harmony"` and only
fails inside `run_integrate`'s `if sample_key is None: raise IntegrationError`
guard.
**Fix:** needs design — either auto-degrade to "none" when no sample_key
candidate is found (note in decision log), or surface the issue at the start
of `analyze` before QC runs. The error message also points at `--sample-key`,
which `scellrun analyze` does not expose (only `scellrun scrna integrate` does).
At minimum, soften the error to mention `--method none`.

---

## [#002] markers CSV / report top rows are zero-vs-near-zero noise

**Stage:** markers
**Severity:** P0 (silent wrongness — what a clinician reads first is the wrong gene list)
**Reproducer:** end-to-end analyze run; open `03_markers/report.html`.
**What happened:** every cluster's top 5–10 marker rows look like:
```
cluster gene log2fc  pct_in_cluster pct_other pval   pval_adj
0       C11orf91 21.3   0.861       0.054     0.914  1.0
0       GABBR2   21.3   0.796       0.652     0.914  1.0
0       SEC14L5  21.3   0.351       0.232     0.914  1.0
0       CLDN34   21.2   0.357       0.260     0.914  1.0
0       AL033504.1 20.8 0.401       0.401     0.914  1.0
0       PTPRB    7.9    ...         ...       ~1e-187 ~1e-185   <-- real marker
0       RAMP3    7.7    ...         ...       0.0    0.0        <-- real marker
0       CD34     7.4    ...         ...       ~1e-270 ~1e-267   <-- real marker
```
The genuine cluster-0 markers (PTPRB / RAMP3 / CD34 — endothelial signature)
are buried below 9–15 rows of nonsense.

**What should happen:** top-of-report markers should be statistically significant
genes, not log2-of-near-zero ratios with NS p-values.

**Suspected cause:** `src/scellrun/scrna/markers.py:148` sorts the per-cluster
DataFrame by `["cluster", "log2fc"]` desc. Scanpy's `rank_genes_groups` with
log-normalized `.raw` produces enormous log2fc values for genes where
`mean_out` is ~0 (log2 of mean+1 of mean+1 of small values). These pass the
default `logfc_threshold=1.0` and `pct_min=0.25` filters, then sort to the
top by log2fc, hiding the real markers.

**Fix:** easy: sort by `pval_adj` ascending (with `log2fc` desc as tiebreak),
which matches the de-facto Seurat "top markers" definition and pushes
statistically-significant genes up. Alternatively, sort by scanpy's `scores`
(Wilcoxon z-statistic) which is what scanpy itself ranks by — and which the
annotate stage already implicitly uses (annotate's panel match works correctly,
because it pulls `rg["names"]` in scanpy's native order, not log2fc order).

A complementary fix would be to add a `pval_adj < 0.05` filter alongside
`logfc>=1` and `pct>=0.25`, but that's a default-threshold change so I'm
deferring it (per "do not fix yourself: anything that changes default
thresholds").

---

## [#003] joint-disease panel auto-pick assigns "ProC" to every score-0 cluster

**Stage:** annotate
**Severity:** P0 (silent wrongness — labels are confidently wrong for >50% of clusters)
**Reproducer:** same end-to-end run, open `04_annotate/annotations.csv`.
**What happened:** clusters whose top markers are clearly non-chondrocyte
(NK, T, B, plasma, mast, macrophage, pericyte, endothelial) are silently
labeled "ProC" with score 0.0 / margin 0.0:
```
cluster top_markers                                     panel_label score
11      RGS5, MCAM, ACTA2, TAGLN, MYL9, NOTCH3 ...     ProC        0.33  ← pericyte
13      MT1X, MT2A, RPS27, IL32, ZFP36L2 ...           ProC        0.00  ← stress / T-cell
19      NKG7, GNLY, GZMB, PRF1, GZMH ...               ProC        0.00  ← NK cell
21      KLRB1, IL7R, LTB, IL32, CD69 ...               ProC        0.00  ← T cell
23      MZB1, JCHAIN, IGHG1, CD79A, XBP1 ...           ProC        0.00  ← plasma cell
24      CPA3, TPSAB1, TPSB2, MS4A2, CTSG ...           ProC        0.00  ← mast cell
```
The deterministic-panel match scored every cluster against
`chondrocyte_markers` (Fan 2024 11 chondrocyte subtypes only). For non-
chondrocyte clusters, every panel group scores 0.0, and
`_best_panel_match` returns `sorted_labels[0]` — which because of dict
iteration order and `sorted(..., reverse=True)` on a tie of zeros is
deterministically "ProC" (first key in the panel dict).

**What should happen:** when every panel group scores 0 for a cluster, the
label should be "Unassigned" (or similar), not the alphabetically/insertion-
order first key. AND — for a subchondral-bone or synovium dataset, the
auto-picked panel `chondrocyte_markers` is the wrong panel: `celltype_broad`
is what should run.

**Suspected cause:** two separate things compound:
1. `_best_panel_match` in `src/scellrun/scrna/annotate.py:134-141` doesn't
   check whether `best_score == 0.0`. It just picks the first by `sorted(...,
   reverse=True)`. With all zeros, dict insertion order wins → `ProC`.
2. `_select_panel` always prefers `chondrocyte_markers` over `celltype_broad`
   when both are present, which is wrong for subchondral-bone or synovium
   tissue. This is correctly flagged by the v0.8 `annotate_self_check` (immune
   marker hint counter), but v0.8 isn't released — v0.7 silently picks the
   chondrocyte panel.

**Fix:** easy for #1: `_best_panel_match` should return `("Unassigned", 0.0,
0.0)` when `best_score == 0.0`. That's a cosmetic-but-honest change, no
default thresholds touched. For #2, defer to maintainer since it changes the
panel-pick policy (and v0.8's self-check already addresses it via
`annotate_panel_tissue_mismatch`).

---

## [#004] Scrublet auto-threshold collapse → "0 doublets" silently

**Stage:** qc
**Severity:** P1 (silent wrongness — user thinks doublets are clean)
**Reproducer:** end-to-end run, look at `01_qc/report.html` Summary panel.
**What happened:**
```
Cells in              12,451
Cells passing QC      10,835 (87.0%)
Doublets flagged      0
```
With a 12k-cell 10x v3 sample, expected doublet rate is 6–10%, so we'd expect
~750–1,200 flagged. Instead 0. Scrublet runs successfully; its `doublet_score`
distribution is sane (median 0.04, 95th pct 0.11, max 0.52); the issue is
scrublet's `predicted_doublet` thresholding. Scrublet's auto-threshold uses a
bimodal-mixture detector that quietly degrades to "threshold > max" when
the simulated-doublet score distribution doesn't form a clean second peak,
which happens routinely on heterogeneous well-separated cell mixes.

**What should happen:** when `predicted_doublet.sum() == 0`, the QC report
should call out that scrublet's auto-threshold collapsed (either no threshold
was found, or all cells fell below it) and suggest a manual cutoff at the
N-th percentile of `doublet_score`. Or just print the score distribution next
to the count.

**Suspected cause:** `src/scellrun/scrna/qc.py:196-202`:
```python
if flag_doublets:
    try:
        sc.external.pp.scrublet(adata)
        obs["scellrun_doublet_flag"] = obs["predicted_doublet"]
        n_doublets = int(obs["scellrun_doublet_flag"].sum())
    except Exception:
        obs["scellrun_doublet_flag"] = False
```
- Exception handler silently swallows everything, including the
  `auto-threshold-not-found` UserWarning scrublet emits.
- No surfaced indication when `n_doublets == 0` whether scrublet ran clean
  (real zero) or auto-threshold collapsed.

**Fix:** easy: when `n_doublets == 0` AND the doublet_score distribution has
non-trivial spread (e.g. max > 0.3), add an `auto_threshold_failed` flag to
the QCResult and surface it in the report's Summary as "Doublets flagged: 0
(scrublet auto-threshold did not converge — set --doublet-threshold manually;
see score distribution histogram)". Bare-except → log it, don't swallow.

---

## [#005] _pick_best_resolution falls through to fragmented resolution

**Stage:** analyze (resolution picker)
**Severity:** P1
**Reproducer:** end-to-end run; open `00_decisions.jsonl`, last `analyze`-stage
row.
**What happened:** decision log says
```json
{"stage": "analyze",
 "key": "chosen_resolution_for_annotate",
 "value": 1.0,
 "rationale": "...picked res=1: n_clusters=25, largest=9.9%, smallest=0.2%, singletons=7"}
```
i.e. the orchestrator picked res=1.0 with **7 fragmented clusters** for
annotation. The integrate quality table shows every resolution has at least 2
singletons (0.1→2, 0.3→2, 0.5→4, 0.8→6, 1.0→7), so the heuristic's
`n_singletons <= 1` filter excludes ALL of them and falls through to "any
with >=2 clusters", picks max n_clusters → res=1.0.
**What should happen:** prefer the resolution with fewest singletons, not most
clusters. With this data res=0.3 (13 clusters, 2 singletons, largest 31.5%)
or res=0.5 (19 clusters, 4 singletons) is the right pick. res=1.0 is the
worst available choice for downstream annotation.
**Suspected cause:** `src/scellrun/analyze.py:80-115` — the fallback branch
ignores singletons entirely and just picks max n_clusters, so the user gets
the most-fragmented resolution when none of them clear the strict bar.
**Fix:** easy: rewrite the fallback to "pick the resolution with the fewest
singletons; on ties, pick the one with the smallest largest_pct" (i.e., least
fragmented + best-balanced). Not a default-threshold change — it's a tweak
to the picker logic.

---

## [#006] integrated.h5ad / annotated.h5ad ~2 GB each

**Stage:** integrate, annotate
**Severity:** P2
**Reproducer:** `du -h <run_dir>/02_integrate/integrated.h5ad`
**What happened:** integrated.h5ad = 1.9 GB, annotated.h5ad = 1.9 GB. Total
run dir ≈ 4 GB on disk for ONE 12k-cell sample.
```
adata.X dtype: float64, shape: (10835, 20713), size: 1795 MB
adata.raw.X dtype: float32 (sparse, smaller)
```
**What should happen:** the integrated h5ad doesn't need the full
20k-gene scaled float64 dense matrix attached. PCA + neighbours + UMAP +
leiden are all stored in `obsm` / `obsp` / `uns` and used downstream;
`adata.X` could either be discarded (annotate / markers re-read `.raw`) or
kept sparse log-normalized (matching `.raw`).
**Suspected cause:** `src/scellrun/scrna/integrate.py:320` `sc.pp.scale(adata,
max_value=10)` densifies to float64 across all 20k genes (only 2k HVGs are
needed for PCA) and `write_artifacts` writes the whole thing.
**Fix:** needs design — either restrict scaling to HVGs (`subset=True`),
or after PCA replace `adata.X` with `adata.raw.X` before writing (recovers
the log-normalized sparse representation), or just write float32 instead of
float64 for the scaled matrix. v0.7 ships as-is is a footgun for any
multi-sample run on a constrained disk.

---

## [#007] index.html lacks any analysis summary

**Stage:** report
**Severity:** P2
**Reproducer:** open `05_report/index.html`.
**What happened:** the top-level index page contains:
- a Decision summary table (good)
- a Run manifest table with timestamps and CLI params (good)
- a Stages list with links to per-stage `report.html` (ok)
- a generic Provenance trail blurb

It does NOT contain:
- QC pass-rate numbers
- cluster counts per resolution
- the chosen-resolution + n_clusters + n_markers headline
- the annotation table (cluster → label) — for which we already wrote
  `annotations.csv` and the per-stage HTML

A clinician opening the index expects to see "X cells passed QC, Y clusters,
Z labeled cell types" without clicking 4 levels deep.

**What should happen:** the top-level index should embed at-a-glance numbers
from each stage (3-4 lines per stage) and the annotation cluster-label table.

**Suspected cause:** `src/scellrun/report.py` builds the index from the
`00_run.json` manifest only; it doesn't re-read each stage's artifacts.

**Fix:** needs design — adding stage-summary blocks on the index requires
the report builder to read each stage's CSV / JSON artifacts and pick the
right numbers. Defer to maintainer.

---

## [#008] Empty `<h2>Decision summary</h2>` and JSON-Path leak in report

**Stage:** report
**Severity:** P2
**Reproducer:** open `05_report/index.html`, look at top:
**What happened:**
```
<p style="color:#666;">
  Run directory: <code>/tmp/scellrun_dogfood/scellrun_out/run-20260430-141432</code>
</p>
<h2>Decision summary</h2>
<p class="summary-counts">
  25 decisions logged across 5 stage(s).
  Rows shaded amber were user overrides; <span class="badge badge-ai">ai</span>
  tags an AI-driven choice.
</p>
```
The "Run directory" prints the **absolute disk path on the hospital server**,
which:
- leaks the user's filesystem layout into a deliverable
- becomes meaningless once the report is moved (printed-to-PDF, emailed, etc.)
- hardcodes `/tmp/...` even though most users will run from `~/scellrun_out/`

Also: there's no `<h2>` for the per-stage decision tables (just
`<div class="stage-head">qc</div>`), so the document outline reads
[Decision summary, Stages, Run manifest, Provenance] with the per-stage
decisions invisible to a screen reader / TOC tool.

**What should happen:** show the run-dir as a relative path ("this report is
in `05_report/`; sister stages are in `01_qc/` ... `04_annotate/`") or hide
the absolute path entirely.

**Suspected cause:** `src/scellrun/report.py` template uses `run_dir` as-is.

**Fix:** easy — render only `run_dir.name` in the index (e.g. "run-20260430-141432")
plus a relative `..` reference to the parent. No design change.

---

## [#009] Pipeline noise: deprecation warnings dump to stdout/stderr

**Stage:** all
**Severity:** P2 (friction — user thinks the run is broken)
**Reproducer:** any `scellrun analyze` run on scanpy 1.11.5.
**What happened:** the run log is dominated by:
```
FutureWarning: Function scrublet is deprecated; Import from sc.pp instead
FutureWarning: Argument `use_highly_variable` is deprecated, consider using mask_var
UserWarning: zero-centering a sparse array/matrix densifies it.
PerformanceWarning: DataFrame is highly fragmented [×~50 occurrences from rank_genes_groups]
FutureWarning: Series.__getitem__ treating keys as positions is deprecated
  pct_in = float(rg["pts"][c][i]) ...
```
The `pts`/`pts_rest` warning fires from `markers.py:131-132` — every cluster ×
every gene = thousands of warnings. The PerformanceWarning fires from inside
scanpy itself but happens because of how rank_genes_groups stores results;
not user-actionable but spams the log.

**What should happen:** clean stdout/stderr — a tired postdoc on a deadline
shouldn't have to scroll past 2,000 lines of warnings to find the "report:
scellrun_out/..." line.

**Suspected cause:** various scanpy 1.11+ API churn that the v0.7 code
hasn't tracked yet; plus the `rg["pts"][c][i]` integer-index access in
markers.py.

**Fix:** easy:
- swap `sc.external.pp.scrublet(adata)` → `sc.pp.scrublet(adata)` in qc.py
- swap `use_highly_variable=True` → `mask_var="highly_variable"` in
  integrate.py PCA call
- in markers.py:131-132, switch `rg["pts"][c][i]` → use the dict-of-DataFrames
  version: scanpy returns `rg["pts"]` as a DataFrame indexed by gene, so
  `rg["pts"].loc[gene, c]` (or compute pct from groupby once and lookup) is
  the modern path.

I'll ship those three swaps as the safe fixes from this dogfood run.

---

## [#010] Re-running `analyze` after a mid-pipeline failure

**Stage:** runlayout / analyze
**Severity:** P2
**Reproducer:** run `scellrun analyze ...` (default --method harmony, single
sample → fails at integrate). Re-run the same command without `--force`.
**What happened:** `01_qc/` already exists from the failed run, so the second
attempt would die at QC with `StageOutputExists`. We worked around with
`rm -rf scellrun_out`.
**What should happen:** if the previous run aborted mid-pipeline, `analyze`
should either:
- detect the partial run-dir and resume from the last completed stage, or
- write to a NEW run-dir with a fresh timestamp (current behavior on
  default `run_dir` is "use timestamp" → already creates a new dir; but
  if the user reuses an explicit `--run-dir`, no resume).
**Suspected cause:** the orchestrator passes `force` literally to each
`stage_dir(run_dir, ..., force=force)` call. Without `--force`, it re-collides.
Default behavior of `default_run_dir()` IS new timestamp, so this only bites
explicit `--run-dir`.
**Fix:** needs design — auto-resume vs. always-rerun is a UX call.

---

## [#011] Install: 1m16s, no progress hint about heaviness

**Stage:** install
**Severity:** P2
**Reproducer:** `pip install scellrun` in a fresh py3.11 env.
**What happened:** install pulled scanpy + scrublet + scikit-image + zarr +
numba + llvmlite + harmonypy + a long tail of deps. Total: 1m16s on a fast
campus link. No README warning that scrublet → scikit-image → tifffile
chain is heavy. A non-engineer PI would think it's stuck after 30s.
**What should happen:** README mention the heaviness; offer a `pip install
scellrun --no-deps` + conda-recipe alternative; or document the conda-first
install path the PI memo'd in `feedback_install` for slow links.
**Fix:** easy: add a one-line note in the README install section:
"Note: scellrun depends on the scanpy stack + scrublet + harmonypy. First
install can take 1–3 minutes on a typical Linux box; please be patient."

---

## Summary — verdict on v0.7

The pipeline runs end-to-end on real OA data and produces deliverables. But
two P0s mean the deliverables lie:

1. The marker CSVs and HTML report show nonsense at the top of every cluster
   (#002), so a postdoc reading the report walks away with the wrong genes.
2. The annotation step silently labels half the clusters "ProC" when they're
   plainly NK / T / B / plasma / mast / pericyte / endothelial (#003).

Either of these alone is enough to make v0.7 unsafe for a clinician
deliverable without a maintainer review. The deterministic provenance trail
(decision log, per-cell metrics CSV) is honest, but the user-facing surfaces
that summarize them aren't.
