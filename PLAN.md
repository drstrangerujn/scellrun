# scellrun roadmap to v1.0 — "no-brain" experience

This document is the standing brief for any subagent or human picking up
work toward v1.0. Each version below is self-contained: a subagent can
read its section, implement, commit, push, and verify CI without further
clarification.

## Working rules (apply to every version)

1. **Do not pip install anything to the developer's local Python.**
   All testing happens on GitHub Actions CI. Local edits → commit → push →
   wait for CI. The maintainer's local environment is intentionally empty
   of scellrun's deps.
2. **Commit hygiene:** no `Co-Authored-By` trailer, author email is the
   GitHub no-reply (`16444402+drstrangerujn@users.noreply.github.com`).
   Repo-local git config already set; do not override.
3. **CI must be green before tagging.** Tag pushes trigger PyPI release
   via `.github/workflows/release.yml`; only tag a SHA whose CI passed.
4. **No silent failures.** If a step can be skipped (LLM unavailable,
   data shape unexpected), warn loudly and surface in the report.
5. **Reports are first-class.** Every command writes an HTML report
   alongside its data artifacts. EN + ZH templates kept in sync.

## Version map

### v0.6 — `scellrun analyze` one-shot pipeline (1.5h)

**Goal:** a single command runs qc → integrate → markers → annotate → report
end to end. New users do not have to know the stage layout.

**Files to add/edit:**
- `src/scellrun/analyze.py` (new): orchestrator function. Imports
  `run_qc`, `run_integrate`, `run_markers`, `run_annotate`,
  `build_report` and chains them. Reads h5ad once at the top, passes the
  in-memory AnnData through each stage, writes intermediate h5ads at the
  canonical `<run_dir>/<NN_stage>/` paths.
- `src/scellrun/cli.py`: add `app.command("analyze")`. Surface params
  for `--profile`, `--species`, `--tissue`, `--resolutions`, `--ai/--no-ai`
  (defaults to True if `ANTHROPIC_API_KEY` is set, else False), `--lang`,
  `--run-dir`, `--force`. Stage-specific overrides via dotted form
  (`--qc.max-genes 4500` etc.) — accept them and route to the right step.
- `tests/test_analyze.py`: synthetic h5ad → pipeline runs end to end →
  all five stage subdirs exist → 05_report/index.html mentions every
  stage. Skip --ai path in tests (no API key in CI).
- ROADMAP.md: mark v0.6 shipped, document `analyze`.
- README.md: replace the multi-step quickstart with the new one-liner.

**Acceptance:** `scellrun analyze data.h5ad --tissue "OA cartilage"`
produces a populated run-dir and prints a single `index.html` link.
CI green; bump to `0.6.0`; tag `v0.6.0` after CI passes.

### v0.7 — Decision log + AI summary (2h)

**Goal:** the report's first page summarizes every non-trivial decision
the pipeline made. The user reads one screen and understands what
happened.

**Files to add/edit:**
- `src/scellrun/decisions.py` (new): `Decision` dataclass + writer that
  appends to `<run_dir>/00_decisions.jsonl`. Fields: stage, key, value,
  default, source ("user", "auto", "ai"), rationale.
- Each stage module records decisions: profile picked, mt% applied
  (default vs user override), sample-key auto-detected, resolution
  recommended by AI, panel chosen for annotate, etc.
- `src/scellrun/report.py` and templates: top section "Decision summary"
  rendering the JSONL grouped by stage. Highlight any AI calls.

**Acceptance:** A run with no overrides shows ~10-15 decisions with
rationale; the user can answer "why mt% 20?" and "why res 0.5?" without
opening any sub-report.

### v0.8 — Error self-heal + actionable suggestions (2h)

**Goal:** when something looks off (most cells failing QC, no clusters,
harmony silently degrading), the pipeline either auto-corrects or stops
with a concrete recommendation, not a stack trace.

**Files to add/edit:**
- Each stage gains a `self_check()` step that runs at the end:
  - QC: if pass-rate < 30%, suggest the smallest threshold relaxation
    that gets to ≥ 60% based on the sensitivity table.
  - Integrate: if all resolutions yield ≤ 2 clusters, suggest
    `--resolutions aio` for the wider sweep.
  - Annotate: if every cluster's `panel_margin` < 0.05, suggest a
    different `--panel` or `--profile`.
- New `--auto-fix/--no-auto-fix` flag at `analyze` level. Default off.
  When on, applies the suggestion and reruns just that stage.
- Self-check output appended to the decision log with `source="auto"`.

**Acceptance:** Synthetic adversarial inputs (e.g. all-mt > 30%, bad
profile) produce actionable suggestions; with `--auto-fix` on, the
pipeline recovers and continues.

### v0.9 — Streamlit web UI (3h)

**Goal:** a clinician with no terminal experience can drag in a `.h5ad`,
fill three fields, click run, read the report.

**Files to add/edit:**
- `src/scellrun/web/app.py` (new): Streamlit app.
  - File upload widget (or path input on the server)
  - Form: tissue, species, profile (dropdown of registered profiles),
    optional Anthropic API key (text input, masked, stored in session)
  - "Run analysis" button → runs `scellrun analyze` as a subprocess →
    streams stdout to a log panel → opens the resulting `index.html`
    inline via `streamlit.components.v1.html`.
- `pyproject.toml`: add `web` extra: `streamlit>=1.30`.
- `Dockerfile`: ships the Streamlit server on port 8501.
- README: "Web UI" section with `pip install scellrun[web]` +
  `scellrun web` (new CLI subcommand that just launches the app).

**Acceptance:** `pip install scellrun[web] && scellrun web` brings up
the UI; uploading a 10x mtx tarball or .h5ad produces the same report
the CLI does.

### v1.0 — Distribution polish (3h)

**Goal:** the project stops being a moving target. APIs frozen, install
paths broad, docs ready for a wider audience.

**Files to add/edit:**
- `conda-forge` feedstock submission: write `meta.yaml`, open the PR
  to `staged-recipes`. Document in `RELEASING.md`.
- `Dockerfile` (already created in v0.9): publish to GHCR via a new
  `release.yml` job.
- `docs/quickstart.md`: 5-minute tutorial including
  `scellrun analyze` + screenshots of the report.
- `docs/contributing.md`: how to add a profile, how to add a marker
  panel, how to add a stage.
- `pyproject.toml`: pin minimum scanpy / anndata versions explicitly to
  what we test on.
- Stable API: any function used outside `cli.py` and outside the
  `analyze` orchestrator promoted to `scellrun.api` re-exports;
  `_internal_*` for everything else.

**Acceptance:** `pip install scellrun==1.0.0` from a fresh conda env
runs the full demo; `scellrun analyze` + Streamlit UI both work; the
quickstart is readable by a non-Python user.

## Subagent dispatch model

Each version above can be handed off to a single `general-purpose`
subagent with three things in the prompt:

1. The version's spec (paste the version block from this PLAN).
2. The "Working rules" block above.
3. A line: "When done, commit + push. Do not tag — the maintainer will
   tag after reviewing CI green. Report back: list of files changed,
   commit SHA, CI run URL."

The subagent works in `isolation: "worktree"` so its branch is isolated;
the maintainer merges to main and tags.
