# Quickstart — scellrun in five minutes

This is for the postdoc or clinician who has cellranger output (or an h5ad)
and wants a single-command analysis with a defensible report. If you'd
rather see how an LLM agent drives scellrun for you, read
[`agent-demo.md`](agent-demo.md) first — most users never type any of these
commands themselves.

## What you need

- A Linux or macOS box with conda (or `uv`).
- 10x cellranger output **or** a `.h5ad`.
- ~3-5 minutes for the first install (scanpy + scrublet + harmonypy stack
  is heavy on the first pull).

## 1. Install

```bash
conda create -n scellrun python=3.11 -y
conda activate scellrun
pip install scellrun
scellrun --version   # expect 1.0.0+
```

If `pip install` looks stuck after 30s, it isn't — scellrun's deps include
scikit-image and tifffile, which are large. First pull is 1-3 min on a
typical campus link.

If your environment uses a regional pip mirror (Tsinghua tuna, Aliyun,
…), the mirror may lag PyPI by ~1 hour after a fresh release. If
`scellrun --version` shows an older version than what's on
[PyPI](https://pypi.org/project/scellrun/), force the canonical index
once:

```bash
pip install --upgrade --index-url https://pypi.org/simple/ scellrun
```

## 2. Convert (skip if you already have an h5ad)

```bash
scellrun scrna convert path/to/cellranger_outs -o data.h5ad
```

`convert` accepts the `outs/filtered_feature_bc_matrix/` directory layout,
the parent `outs/` directory, or a Seurat-style mtx triplet
(`matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz`).

## 3. Analyze

```bash
scellrun analyze data.h5ad \
    --profile joint-disease \
    --tissue "OA cartilage" \
    --no-ai \
    --auto-fix \
    --lang en
```

Flags:

- `--profile` picks the threshold + marker-panel bundle. v1.0 ships
  `default` and `joint-disease`; see `scellrun profiles list`.
- `--tissue` is a free-text string. It's forwarded to the AI summariser
  (when `--ai` is on) and rendered into the report header. Pick the
  cleanest two-to-four-word description of your sample.
- `--no-ai` runs deterministic-only. Set `ANTHROPIC_API_KEY` and drop
  this flag for an LLM second-pass annotation + summary.
- `--auto-fix` lets the orchestrator apply the cheapest stage self-check
  suggestion once (e.g. relax mt% if pass-rate falls below 60%, swap
  panels if every cluster matches at margin <0.05). The retry is logged
  as a `source="auto"` row in `00_decisions.jsonl`.
- `--lang en` or `--lang zh` for English / Chinese reports.

## 4. Read the report

```bash
ls scellrun_out/run-*/
# 00_run.json        run manifest
# 00_decisions.jsonl decision log (one JSON per row)
# 01_qc/             per-cell metrics + qc.html
# 02_integrate/      umaps + cluster table + integrated.h5ad
# 03_markers/        per-resolution markers.csv + report.html
# 04_annotate/       annotations.csv + report.html
# 05_report/         index.html  ←  this is the one to open
```

Open `05_report/index.html` in a browser. The first block is **At a
glance** — QC pass rate, integration method, panel, chosen resolution,
top cell-type labels, any self-check codes that fired. Below that:
**Decision summary** — every non-trivial choice scellrun made, grouped
by stage, with rationale strings and `auto` / `user` / `ai` source
badges. Below that: links to the per-stage reports.

For a printable PDF, use the browser's print-to-PDF on `05_report/index.html`.

## 5. What "defensible" looks like in the report

Every numeric threshold in the QC stage prints with its rationale text
(e.g. "max_pct_mt=20.0 — joint-tissue tolerance per Liu-lab AIO"). The
integrate stage prints the per-resolution cluster-quality table
(n_clusters, largest cluster %, smallest %, n_singletons). The annotate
stage prints, per cluster, the panel match label, score, margin, top
markers, and — if `--ai` was on — the LLM second opinion side-by-side
with the deterministic call. **scellrun never silently picks one over
the other**; the user sees both and decides.

## A note on what "agent-driven" looks like

If you're an LLM agent driving scellrun for a human user — read
[`agent-demo.md`](agent-demo.md) and `skills/scellrun/SKILL.md`. The
SKILL.md is the operational guide; the agent-demo is a worked example
on real OA cartilage data.

## When to call for help

- The QC stage flags a pass-rate <30%: that's a sample-quality issue,
  not a threshold issue. `--auto-fix` will not save it.
- Annotate labels every cluster `Unassigned`: the panel doesn't match
  the tissue; pass `--panel celltype_broad` (or build a profile module
  for your tissue, see [`contributing.md`](contributing.md)).
- The integrated.h5ad is enormous (>1 GB per 10k cells): expected on the
  first run; cast to float32 happens at write time as of v0.9.1, so
  re-run if you're on an older version.

If something doesn't fit any of the above, the run-dir is the bug
report — `00_run.json` + `00_decisions.jsonl` + `v1demo.log` is
sufficient context for a maintainer.
