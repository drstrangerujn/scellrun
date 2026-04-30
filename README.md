# scellrun

[![CI](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml/badge.svg)](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/scellrun.svg)](https://pypi.org/project/scellrun/)

scellrun stops two analysts — or two LLM agents — from getting two different answers on the same single-cell data.

## Why this exists

Single-cell analysis has a documented reproducibility problem. From the field's own retrospective ([Perspectives on rigor and reproducibility in single cell genomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC9122178/)): *"in my group's experience, it is not unusual for reanalysis to find 20% fewer or more clusters in datasets"* — same raw data, different analyst, different answer. The same review notes that of ~50 high-impact single-cell papers surveyed, *"just a handful"* reported any external validation. Most of the choices that drive that 20% divergence — mt% ceiling, HVG count, integration method, clustering resolution, panel pick — are made ad-hoc in a notebook and never written down.

A modern LLM agent handed an `.h5ad` and `scanpy` will write working code and produce a report. That solves *"can the work happen"*. It does not solve *"will two agents on the same data produce the same answer"*, *"can a reviewer reconstruct why mt% was 20 and not 10 six months later"*, or *"is the panel choice consistent with this lab's working practice on this tissue"*. Vanilla agents improvise thresholds, do not record the rationale in any machine-readable form, and have no way to encode the consensus a clinical-bioinformatics team has built over years of dogfooding.

scellrun fills that gap. Every threshold has a tested default with a one-sentence rationale; every choice the pipeline makes — auto, user-override, or LLM-recommended — is appended to a `00_decisions.jsonl` file you can grep; tissue-specific working practice ships as `profiles/` (cartilage today, contribute yours); each stage runs a self-check against PI-defined trigger thresholds and surfaces an actionable suggestion before the user sees the downstream finding. Different layer from a workflow manager: if you need orchestration across a cluster, use [nf-core](https://nf-co.re/scrnaseq); scellrun is what you call from inside one of those. Not a replacement for scanpy either — it calls scanpy under the hood, with opinionated parameters and a decision log on top.

## Who this is for

- **The LLM agent (Claude Code, Hermes, Codex) handling a clinician's request.** This is the primary user. The agent ssh's to the data, runs `scellrun analyze`, reads the artifacts, and translates. `skills/scellrun/SKILL.md` is the operational guide it reads.
- **The clinician-bioinformatics team that wants every project to look the same in a report.** Same QC layout, same decision table, same provenance trail across samples, students, and rotations.
- **The reviewer asking "why mt% 20?"** The answer is `00_decisions.jsonl` line 14, verbatim: *"mt% ceiling 20.0% — joint tissue is stress-prone, the textbook 10% silently drops real chondrocytes (PI cohort 2024-2026, AIO PM=20)"*.

## Quick start

```bash
conda create -n scellrun python=3.11 -y
conda activate scellrun
pip install scellrun

scellrun analyze data.h5ad --tissue "OA cartilage"
# → ./scellrun_out/run-<ts>/05_report/index.html
```

Don't have an `.h5ad`? Cellranger output works directly:

```bash
scellrun scrna convert path/to/cellranger_outs -o data.h5ad
scellrun analyze data.h5ad --tissue "OA cartilage"
```

Add `--lang zh` for a Chinese report. Add `--profile joint-disease` if the tissue is cartilage / synovium / subchondral bone (auto-loads the Fan 2024 chondrocyte panel and tightens hb% for avascular cartilage). Walkthrough in [`docs/quickstart.md`](docs/quickstart.md); contribution notes in [`docs/contributing.md`](docs/contributing.md).

## How an agent uses this

Drop [`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md) into your agent's skills directory and the agent will know which command maps to which user intent, how to read the decision log, when to surface a self-check trigger before answering, and which profile to pick by tissue keyword. [`docs/agent-demo.md`](docs/agent-demo.md) is a verbatim transcript of a Claude Code agent running scellrun end-to-end on real OA cartilage scRNA data — including the agent quoting the decision log when the user asks "why resolution 0.3?" and switching panels when the deterministic call is wrong.

## What's in the decision log

`00_decisions.jsonl` is the single source of truth for every non-trivial choice the pipeline made. One JSON object per line. Sample shape (real, from `docs/v1demo/decisions.jsonl`):

```jsonl
{"schema_version":1,"stage":"qc","key":"max_pct_mt","value":20.0,"default":20.0,"source":"auto","rationale":"mt% ceiling 20.0% — joint tissue is stress-prone, the textbook 10% silently drops real chondrocytes (PI cohort 2024-2026, AIO PM=20)","fix_payload":null,"attempt_id":"cae89793d0f9470a9c7f38894928f304","ts":"2026-04-30T15:20:32+00:00"}
{"schema_version":1,"stage":"analyze","key":"method_downgrade","value":"none","default":"harmony","source":"auto","rationale":"no sample/batch column in obs — single-sample input; auto-downgraded --method from harmony to none","fix_payload":null,"attempt_id":"cae89793d0f9470a9c7f38894928f304","ts":"2026-04-30T15:18:51+00:00"}
{"schema_version":1,"stage":"analyze","key":"chosen_resolution_for_annotate","value":0.3,"default":null,"source":"auto","rationale":"fewest singletons → most balanced (every resolution fragmented) — picked res=0.3: n_clusters=13, largest=31.5%, smallest=0.2%, singletons=2","fix_payload":null,"attempt_id":"cae89793d0f9470a9c7f38894928f304","ts":"2026-04-30T15:20:52+00:00"}
{"schema_version":1,"stage":"analyze","key":"annotate.auto_panel","value":"chondrocyte_markers","default":null,"source":"auto","rationale":"fine-subtype panel preferred (chondrocyte hits dominate or are tied)","fix_payload":null,"attempt_id":"cae89793d0f9470a9c7f38894928f304","ts":"2026-04-30T15:24:19+00:00"}
{"schema_version":1,"stage":"annotate","key":"panel","value":"chondrocyte_markers","default":null,"source":"user","rationale":"--panel 'chondrocyte_markers' forced","fix_payload":null,"attempt_id":"cae89793d0f9470a9c7f38894928f304","ts":"2026-04-30T15:24:49+00:00"}
```

Every choice the pipeline made, with a one-sentence rationale, in a file you can grep. `source` is one of `auto` (a built-in heuristic), `user` (a CLI override), or `ai` (an LLM call). `fix_payload` is non-null only on self-check `*.suggest` rows — it carries the structured fix the orchestrator can mechanically apply when `--auto-fix` is on. `attempt_id` groups rows by invocation. The full schema is in [`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md).

## Profiles

A profile is community-encoded working practice for a tissue domain — defaults plus marker panels — in one Python file. v1.0 ships two:

- **`default`** — fresh-tissue 10x v3 chemistry, joint-tissue-aware mt% ceiling at 20% (the textbook 10% silently drops real chondrocytes; the OARSI working group ceiling is 20%).
- **`joint-disease`** — same QC plus tighter hb% (cartilage is avascular), the Fan 2024 chondrocyte 11-subtype panel, and a 15-group `celltype_broad` panel. Auto-swaps from chondrocyte to broad when the data is immune-rich (subchondral bone, infiltrated synovium, joint fluid) so the report doesn't blindly mis-label pericytes / plasmacytoid DCs / osteoclasts as chondrocyte subtypes.

```bash
scellrun profiles list
scellrun profiles show joint-disease   # prints thresholds + panels
```

If your tissue or disease has working practice that diverges from the defaults, contribute a profile — one Python file under `src/scellrun/profiles/`.

## Roadmap

v0.1 → v1.0.1 has shipped: per-stage QC / integrate / markers / annotate, the `analyze` one-shot, the decision log, self-check + `--auto-fix`, the joint-disease profile, panel auto-pick, single-sample auto-downgrade, agent demo, Dockerfile, and the v1.0.1 SKILL.md sync. The CLI surface is frozen for the v1.x series; new stages and profiles land additively. Post-v1.0 directions tracked in [`ROADMAP.md`](ROADMAP.md): conda-forge feedstock, registry-pushed Docker image, bulk RNA-seq subcommand, metabolomics composite scoring, proteomics integration.

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgements

Defaults trace to the in-house R AIO pipeline (Liu lab) and clinician-bioinformatics working practice for OARSI / MSK research. Built with assistance from Claude (Anthropic).
