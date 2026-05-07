# scellrun — agent notes

PI: 刘希宇 (xiyu.l@icloud.com). All code authored by Claude under PI direction.

## Positioning (don't drift)

scellrun ≠ scientific skills. Skills are LLM-facing prompt scaffolds. scellrun is a CLI that produces deliverables (HTML/PDF reports) for human researchers. When in doubt: ask "would a tired postdoc on a deadline want this?" — that's the user.

## Hard rules

1. **Opinionated defaults.** Every threshold/parameter has a one-line rationale in code. If you can't justify the default in one sentence linking to a paper or PI memory, don't ship it.
2. **Report-first.** Each command produces a report artifact (HTML at minimum), not a pile of plots. Provenance trail in three tiers: data / inference / literature (matches `feedback_report_rules`).
3. **Smoke-test on real data before claiming done.** Do not commit a command that hasn't been run against an actual h5ad from `/root/oc_analysis/` or another real source. Heavy-deps tests (scanpy/anndata stack) run on the hospital server in a conda env, never on the local German server's system Python.
4. **v1.3 frozen surface — scRNA only.** As of v1.3.1 the CLI surface is FROZEN for the v1.x series: qc, integrate, markers, annotate, analyze, review, export, profiles, scrna_full_report. No new public commands. Bug fixes and additive profiles only. Bulk RNA / metabolomics / proteomics extensions are deferred to a future v2.0 (out of scope for the v1.x cycle); see ROADMAP.md "deferred" section. Rationale: v0.1 charter asked for "fast + stable scRNA analysis"; v1.3 already met that bar, every additional surface dilutes the deliverable.
5. **Domain-agnostic surface, profile-specific defaults.** Code paths must work for any human scRNA h5ad. Domain-specific knowledge lives in `scellrun.profiles.*` (default + joint-disease + tumor + brain + kidney; cold-validated only on joint-disease). Don't bake "OA" or "joint" into command names, public APIs, or generic tests.

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
