# scellrun — agent notes

PI: 刘希宇 (xiyu.l@icloud.com). All code authored by Claude under PI direction.

## Positioning (don't drift)

scellrun ≠ scientific skills. Skills are LLM-facing prompt scaffolds. scellrun is a CLI that produces deliverables (HTML/PDF reports) for human researchers. When in doubt: ask "would a tired postdoc on a deadline want this?" — that's the user.

## Hard rules

1. **Opinionated defaults.** Every threshold/parameter has a one-line rationale in code. If you can't justify the default in one sentence linking to a paper or PI memory, don't ship it.
2. **Report-first.** Each command produces a report artifact (HTML at minimum), not a pile of plots. Provenance trail in three tiers: data / inference / literature (matches `feedback_report_rules`).
3. **Smoke-test on real data before claiming done.** Do not commit a command that hasn't been run against an actual h5ad from `/root/oc_analysis/` or another real source.
4. **Keep v0.1 scope tight.** scRNA QC only. No Harmony, no annotation, no metabolomics until v0.2+.
5. **Domain-agnostic surface, profile-specific defaults.** Code paths must work for any human scRNA h5ad. Domain-specific knowledge lives in `scellrun.profiles.*` (currently only `default`). Don't bake "OA" or "joint" into command names, public APIs, or generic tests.

## Code conventions

- typer for CLI (subcommand groups)
- jinja2 for HTML reports
- Default output dir: `./scellrun_out/<command>/<timestamp>/`
- All numeric thresholds in `scellrun.defaults` module with rationale docstrings
- No mocking real data — tests run against tiny real h5ad fixtures (TBD)

## Memory anchors (apply in implementation)

- `feedback_neutrophil_in_snRNA` — annotation step must verify with 5-marker check
- `feedback_score_genes_z_and_mean` — score_genes outputs report z + mean
- `feedback_metaboanalyst_rsd_filter` — RSD prefilter table preserved
- `feedback_composite_score_vs_multiple_testing` — composite score, not N FDRs
- `feedback_matplotlib_cjk_hospital` — fontManager path for any CJK figure
- `feedback_no_traffic_light` — no 🔴🟡🟢 in reports
