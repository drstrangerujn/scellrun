# scellrun

[![CI](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml/badge.svg)](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Your LLM agent uses this.** scellrun is a CLI an agent (Claude Code,
> Hermes, Codex, …) drives end-to-end on a researcher's behalf — convert
> cellranger output, QC, integrate, cluster, annotate, render an HTML
> report — and quotes the decision log back when the researcher asks
> "why mt% 20?" or "why res 0.5?". The deliverable is a defensible
> analysis the user can read in five minutes, not a pipeline they have
> to learn.

The agent-driven story is the canonical one. See
[`docs/agent-demo.md`](docs/agent-demo.md) for a verbatim transcript on
real OA cartilage scRNA data. The agent's operational guide is
[`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md).

You can also drive scellrun yourself if you prefer:

```bash
# 1. one-time: a clean Python environment so scellrun's deps don't collide with anything else on your machine
conda create -n scellrun python=3.11 -y
conda activate scellrun

# 2. install
pip install scellrun

# 3. run the full pipeline in one shot — qc → integrate → markers → annotate → report
scellrun analyze data.h5ad --tissue "OA cartilage"
# → opens with a single index.html link you can drop into a browser
```

Don't have a `.h5ad`? Cellranger output works directly:

```bash
scellrun scrna convert path/to/cellranger_outs -o data.h5ad
scellrun analyze data.h5ad --tissue "OA cartilage"
```

Want a Chinese report? Add `--lang zh`. Want a clinician walkthrough?
See [`docs/quickstart.md`](docs/quickstart.md). Want to add a profile or
a stage? See [`docs/contributing.md`](docs/contributing.md).

**Status:** v1.0.0. The CLI surface is frozen for the v1.x series; new
stages and profiles land additively.

## Who this is for

- A clinician or rotating student looking at scRNA-seq data for the first time
- A postdoc on a deadline who doesn't want to write QC boilerplate
- A bioinformatics core that needs every project to look the same in a report
- An LLM coding agent that should reach for a real tool instead of re-deriving thresholds

If you already have your own pipeline, you don't need scellrun. If you don't,
or you're tired of re-litigating the same decisions every project, this is for you.

## Why this exists

Tools like `scanpy`, `pyDESeq2`, `decoupler`, `SCENIC` and `nf-core` are excellent
and intentionally unopinionated. Every new project re-litigates the same decisions:

- Which mt% threshold? Doublet filter before or after batch correction?
- The composite score correlates with disease — Spearman ρ on each of 18 clinical features, or one well-chosen dichotomy?
- What does the QC report actually need to look like for a reviewer to stop asking?

scellrun answers these once, in code, by encoding the working practice of a
clinician + bioinformatics team. Defaults aren't "neutral" — they're a real
position, made for real reasons, with the rationale rendered into every report.

You can override anything. But if you don't override, you get a defensible
analysis on day one.

## Where this fits

| Tool                          | What it is                                                         |
| ----------------------------- | ------------------------------------------------------------------ |
| Anthropic scientific skills   | Prompt scaffolds that teach an LLM how to call `scanpy` etc.        |
| Bioconda recipes              | Packaging                                                          |
| `nf-core` pipelines           | General-purpose, infrastructure-heavy, vendor-neutral workflows     |
| **`scellrun`**                | An *opinionated, report-first* CLI optimized for low barrier to entry |

## What v1.0 ships

- `scellrun analyze <h5ad>` — one-shot pipeline (qc → integrate →
  markers → annotate → report). Writes a deterministic decision log
  (`00_decisions.jsonl`) and a top-level `05_report/index.html` with an
  At-a-glance block + per-stage decision tables.
- `scellrun scrna {qc,integrate,markers,annotate}` — per-stage
  commands for deep customization.
- `scellrun scrna convert <cellranger_dir> -o data.h5ad` — convert
  10x cellranger / Seurat-mtx output.
- Self-check + `--auto-fix` — each stage flags actionable issues
  (low QC pass-rate, panel-tissue mismatch, fragmented clustering)
  and the orchestrator can apply the cheapest fix once.
- Profiles: `default` and `joint-disease` (Fan 2024 chondrocyte
  11-subtype + 15-group celltype_broad panel).
- Optional `--ai`: Anthropic-API LLM second opinion on annotation +
  resolution recommendation.
- Distribution: PyPI wheel, Dockerfile, `skills/scellrun/SKILL.md`
  for agent harnesses.

See [`ROADMAP.md`](ROADMAP.md) for the post-v1.0 plan (conda-forge
feedstock, registry-pushed Docker image, bulk RNA-seq, metabolomics
composite scoring, proteomics integration).

## Install

### For users — clean conda env (recommended)

```bash
conda create -n scellrun python=3.11 -y
conda activate scellrun
pip install scellrun
```

If you don't have conda, [miniconda](https://docs.anaconda.com/miniconda/)
takes about 2 minutes to install. Or use [`uv`](https://docs.astral.sh/uv/)
for a faster venv:

```bash
uv venv .scellrun-env && source .scellrun-env/bin/activate
uv pip install scellrun
```

Either way, scellrun ends up in its own environment — your other Python
projects (an old scanpy, Seurat-via-rpy2, etc.) won't be touched.

### For contributors — editable from source

```bash
git clone https://github.com/drstrangerujn/scellrun.git
cd scellrun
conda create -n scellrun-dev python=3.11 -y && conda activate scellrun-dev
pip install -e ".[dev]"
pytest -q          # full suite should be green
```

## Profiles

Different tissue domains have different defaults. v1.0 ships two profiles:

- `default` — fresh-tissue 10x v3 chemistry, joint-tissue-aware mt% ceiling
- `joint-disease` — tighter hb% for avascular cartilage; ships the **Fan 2024 chondrocyte 11-subtype panel** + 15-group broad celltype panel used by `scrna annotate`

```bash
scellrun profiles list
scellrun profiles show joint-disease   # prints thresholds + panels
```

Adding a profile = one Python file under `src/scellrun/profiles/`. If your
community has working practice, [contribute a profile](src/scellrun/profiles/).

## For LLM agents

`skills/scellrun/SKILL.md` is a portable instruction document that teaches
Claude Code, Hermes, Codex, or any markdown-aware agent harness when and how
to invoke scellrun. Install it by symlinking into your agent's skills
directory — see [`skills/README.md`](skills/README.md) for one-liner install
commands per agent.

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgements

Defaults reflect the working practice of a clinician + bioinformatics team
that has shipped these analyses for the OARSI / MSK community. The R "AIO"
pipeline that prefigured scellrun's design is documented in
[`ROADMAP.md`](ROADMAP.md).

Built with assistance from Claude (Anthropic).
