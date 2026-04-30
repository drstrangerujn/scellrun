# scellrun

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **scellrun lowers the bar to running a defensible single-cell analysis.**
> One command, a publication-quality report. You don't need to learn scanpy first.

```bash
pip install scellrun
scellrun scrna qc data.h5ad
# → opens a browser-readable QC report in seconds
```

**Status:** v0.1.0, early alpha. APIs will change without warning. Don't pin against a tag yet.

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

## v0.1 scope

- `scellrun scrna qc <h5ad>` — single-cell QC (mt%, ribo%, hb%, n_genes, doublets)
  with defensible thresholds. Outputs an HTML report with rationale text, a
  per-cell CSV, and an annotated `qc.h5ad` ready for the v0.2 integration step.
  Cells are *flagged* (`obs.scellrun_qc_pass`), never silently dropped.

That's it. See [`ROADMAP.md`](ROADMAP.md) for v0.2-v0.5 (integrate, markers,
annotate, multi-stage report) and the planned `scellrun run --steps ...`
pipeline mode.

## Install (dev)

```bash
git clone https://github.com/drstrangerujn/scellrun.git
cd scellrun
pip install -e ".[dev]"
scellrun --help
scellrun scrna qc /path/to/data.h5ad
```

## Profiles

Different tissue domains have different defaults. v0.1 ships two profiles:

- `default` — fresh-tissue 10x v3 chemistry, joint-tissue-aware mt% ceiling
- `joint-disease` — tighter hb% for avascular cartilage; ships the **Fan 2024 chondrocyte 11-subtype panel** + 15-group broad celltype panel for v0.4 annotation

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
