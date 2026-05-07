# scellrun

[![CI](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml/badge.svg)](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/scellrun.svg)](https://pypi.org/project/scellrun/)

> A scRNA-seq CLI with a decision log so two analysts — or two LLM agents — on the same data give the same answer.

[中文版 README](README_zh.md)

---

## The 30-second pitch

Two analysts on the same scRNA dataset reanalyze it and the cluster count differs by ~20% ([Lähnemann et al., 2020 / PMC9122178](https://pmc.ncbi.nlm.nih.gov/articles/PMC9122178/)). The reason is rarely the science — it's the dozens of tiny choices nobody writes down: mt% ceiling, HVG count, integration method, clustering resolution, panel pick. Six months later nobody remembers which knob got tweaked.

scellrun is what you put on top of [scanpy](https://scanpy.readthedocs.io) so every one of those choices ends up in a single grep-able file with a reason next to it. Same data → same defaults → same output, every time. If a reviewer asks "why mt% 20?" you read line 14 of `00_decisions.jsonl` to them.

It's not a workflow manager (use [nf-core](https://nf-co.re/scrnaseq) for cluster orchestration), and it doesn't replace scanpy — it calls scanpy under the hood with opinionated parameters and an audit trail.

## Quick start

```bash
conda create -n scellrun python=3.11 -y
conda activate scellrun
pip install scellrun

scellrun analyze data.h5ad --tissue "OA cartilage"
# → ./scellrun_out/run-<ts>/05_report/index.html
```

Got cellranger output instead of an `.h5ad`?

```bash
scellrun scrna convert path/to/cellranger_outs -o data.h5ad
scellrun analyze data.h5ad --tissue "OA cartilage"
```

Add `--lang zh` for a Chinese report. Add `--profile joint-disease` for cartilage / synovium / subchondral bone. Walkthrough in [`docs/quickstart.md`](docs/quickstart.md).

## What you actually get

- **Five-stage one-shot**: QC → integrate (Harmony) → markers → annotate → report. One command, one HTML you can email.
- **Decision log** at `<run>/00_decisions.jsonl` — every non-trivial choice with a one-sentence reason. Greppable; `auto`/`user`/`ai` source labels.
- **Five tissue profiles** — `default`, `joint-disease` (Fan 2024 chondrocyte panel), `tumor`, `brain`, `kidney`. One Python file each, contribute yours.
- **Self-check** — each stage detects pathologies (panel mismatch, all-fragmented clusters, single-sample-no-batch) and proposes the cheapest fix. `--auto-fix` applies it.
- **Reviewer loop** — `scellrun review <run>` runs a tiny local Flask app for cluster relabels, threshold tweaks, and notes; `analyze --apply-overrides <json>` re-runs with the human's edits as `source="user"` rows.
- **PDF export** — `scellrun export <run> --format pdf` for publication.

## Who this is for

- **The LLM agent (Claude Code, Hermes, Codex, OpenClaw) handling a clinician's request.** This is the primary user. Drop [`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md) into your agent's skill directory; the agent then knows which command maps to which user intent, how to read the decision log, when to surface a self-check trigger.
- **The clinician-bioinformatics team that wants every project to look the same.** Same QC layout, same decision table, same provenance trail across samples, students, rotations.
- **The reviewer asking "why mt% 20?"** Open `00_decisions.jsonl`; the answer is on line 14, verbatim.

[`docs/agent-demo.md`](docs/agent-demo.md) is a verbatim Claude Code transcript running scellrun end-to-end on real OA cartilage data, including the agent quoting the decision log when the user asks "why res=0.3?".

## What the decision log looks like

```jsonl
{"stage":"qc","key":"max_pct_mt","value":20.0,"default":20.0,"source":"auto",
 "rationale":"mt% ceiling 20% — joint tissue is stress-prone, the textbook 10% silently drops real chondrocytes (PI cohort 2024-2026, AIO PM=20)"}
{"stage":"analyze","key":"method_downgrade","value":"none","default":"harmony","source":"auto",
 "rationale":"no sample/batch column in obs — single-sample input; auto-downgraded --method from harmony to none"}
{"stage":"analyze","key":"chosen_resolution_for_annotate","value":0.3,"source":"auto",
 "rationale":"fewest singletons → most balanced (every resolution fragmented) — picked res=0.3: n_clusters=13, largest=31.5%, smallest=0.2%, singletons=2"}
{"stage":"analyze","key":"annotate.auto_panel","value":"celltype_broad","source":"auto",
 "rationale":"swapped to celltype_broad: chondrocyte_hits=2, broad_hits=9; required >=1.5x margin to keep chondrocyte panel"}
```

`source` is `auto` (built-in heuristic), `user` (a CLI / review override), or `ai` (an LLM call). `attempt_id` groups rows by invocation; `fix_payload` carries the structured fix on self-check `*.suggest` rows. Full schema in [`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md).

v1.3.2 onwards the `chosen_resolution_for_annotate` rationale, panel auto-pick reasoning, and self-check triggers all surface in the HTML report's "At a glance" section so a reader doesn't have to grep the jsonl to learn why a particular resolution / panel got picked.

## Profiles

A profile is community-encoded working practice for a tissue domain — defaults + marker panels in one Python file.

| profile | mt% | hb% | panels | notes |
|---|---|---|---|---|
| `default` | 20% | — | — | fresh-tissue 10x v3 baseline (OARSI ceiling) |
| `joint-disease` | 20% | tight | Fan 2024 11-subtype chondrocyte + 15-group broad | cold-validated; auto-swaps to broad on immune-rich data |
| `tumor` | 20% | — | TISCH/Sun 2021 pan-cancer TME (broad only) | not yet cold-validated |
| `brain` | 10% | — | Tasic/Hodge cortical-hippocampal (broad only) | not yet cold-validated |
| `kidney` | 15% | — | KPMP/Stewart 2019 nephron + immune (broad only) | not yet cold-validated |

```bash
scellrun profiles list
scellrun profiles show joint-disease   # thresholds + panels
```

Tissue or disease working practice diverging from the defaults? Contribute a profile — one Python file under `src/scellrun/profiles/`. See [`docs/contributing.md`](docs/contributing.md).

## Status

**v1.3 frozen surface — scRNA only.** v1.x is now in maintenance mode. Bug fixes and additive scRNA profiles only; no new public commands. Bulk RNA-seq, metabolomics, and proteomics extensions are deferred to a future v2.0; see [`ROADMAP.md`](ROADMAP.md).

CLI surface (`qc` / `integrate` / `markers` / `annotate` / `analyze` / `review` / `export` / `profiles`) is locked for the v1.x series.

## Distribution

- **PyPI**: `pip install scellrun` ([pypi.org/project/scellrun](https://pypi.org/project/scellrun/))
- **ClawHub**: `clawhub install scellrun` for agent skills ([clawhub.ai/skills/scellrun](https://clawhub.ai/skills/scellrun))
- **Docker**: `docker pull ghcr.io/drstrangerujn/scellrun:latest` (v1.0+)

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgements

Defaults trace back to the Liu-lab in-house R AIO pipeline and clinician-bioinformatics working practice for OARSI / musculoskeletal research. The Fan 2024 chondrocyte panel ships under the `joint-disease` profile.
