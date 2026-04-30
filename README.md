# scellrun

> An opinionated, report-first CLI for single-cell and multi-omics analysis. Sane defaults baked in, deliverables out the other end, no twelve-knob refactor needed for every project.

**Status:** v0.1, early alpha. APIs will change without warning. Don't pin against a tag yet.

## Why this exists

The existing tooling — `scanpy`, `pyDESeq2`, `decoupler`, `SCENIC`, `nf-core` — is excellent and intentionally unopinionated. Every new project re-litigates the same decisions:

- Which mt% threshold? Doublet filter before or after batch correction?
- score_genes is significant — should I trust it on a sparse target gene?
- The composite score correlates with disease — do I report Spearman ρ on each of 18 clinical features, or one well-chosen dichotomy?
- What does the QC report actually need to look like for a reviewer to stop asking?

`scellrun` answers these once, in code, by encoding the working practice of a clinician + bioinformatics team. The defaults are not "neutral" — they reflect a real position, made for real reasons, with the rationale rendered into every report.

You can override anything. But if you don't override, you get a defensible analysis on day one.

## Where this fits

| Tool                          | What it is                                                         |
| ----------------------------- | ------------------------------------------------------------------ |
| Anthropic scientific skills   | Prompt scaffolds that teach an LLM how to call `scanpy` etc.        |
| Bioconda recipes              | Packaging                                                          |
| `nf-core` pipelines           | General-purpose, infrastructure-heavy, vendor-neutral workflows     |
| **`scellrun`**                 | A *narrow, opinionated, report-first* CLI for human research data   |

LLM coding agents can call `scellrun` rather than re-deriving QC thresholds. Researchers get a reproducible HTML/PDF report rather than a folder of disconnected plots.

## Profiles (planned)

Different domains have different defaults. v0.1 ships one profile (`default`), good for fresh-tissue 10x v3 chemistry. The architecture is built so additional profiles can be contributed in a single PR — e.g. `--profile=tumor`, `--profile=joint-disease`, `--profile=brain-snrna`. Each profile is a small Python file with documented thresholds; no plugin system, no DSL, just a dataclass.

If your community has working practice, contribute a profile.

## v0.1 scope

- `scellrun scrna qc <h5ad>` — single-cell QC (mt%, ribo%, n_genes, doublets) with annotated, defensible thresholds. Outputs an HTML report with rationale text and a per-cell CSV. Cells are *flagged*, never silently dropped.

That's it. The roadmap is below.

## Install (dev)

```bash
git clone https://github.com/drstrangerujn/scellrun.git
cd scellrun
pip install -e ".[dev]"
scellrun --help
scellrun scrna qc /path/to/data.h5ad
```

## Roadmap

- v0.2 — `scellrun scrna integrate` (Harmony with sane batch QC, integration-quality report).
- v0.3 — `scellrun scrna annotate` (panel-based annotation with marker-verification step that catches the "neutrophil cluster lacking MPO" failure mode).
- v0.4 — `scellrun stats composite` (composite-score scaffolding for clinical–molecular dichotomies, replacing N×FDR with one defensible test).
- v0.5 — `scellrun report` (publication-quality PDF deliverable with the three-tier provenance trail: data / inference / literature).
- v0.6+ — additional profiles, additional assays (bulk RNA-seq, metabolomics, proteomics integration).

## Contributing

The fastest way to be useful: open an issue with a default you'd change and a one-line citation. The defaults document is the spine of the project.

## For LLM agents

`skills/scellrun/SKILL.md` is a portable instruction document that teaches
an agent (Claude Code, Hermes, Codex, etc.) when and how to invoke scellrun.
Install it by symlinking into your agent's skills directory — see
[`skills/README.md`](skills/README.md) for one-liner install commands per agent.

If you'd rather extend scellrun than teach an agent to use it, contribute a
**profile** under `src/scellrun/profiles/` instead.

## License

MIT.

## Acknowledgements

The defaults baked into this tool are the residue of a lot of failed analyses and one or two saved nights. If they save you a similar evening, that's the goal.

Built with assistance from Claude (Anthropic).
