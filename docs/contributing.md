# Contributing to scellrun

scellrun is an opinionated tool. The point of contributing isn't "add
options"; it's "encode another working practice in code so the next
project doesn't re-litigate it." If you have a defensible default and
can write the rationale in one sentence, that's a profile / panel /
stage worth adding.

This page covers the three most common contribution shapes.

## 1. Add a profile

A profile is one Python module under `src/scellrun/profiles/`. It
exports per-assay QC threshold dataclasses and (optionally) marker
panels. It is the only place where domain-specific defaults live —
the rest of scellrun is tissue-agnostic.

Look at `src/scellrun/profiles/joint_disease.py` for a complete
example. The minimum:

```python
# src/scellrun/profiles/<your_tissue>.py
"""
<your_tissue> profile — short tissue description.

Reference: <citation or in-house pipeline this mirrors>
"""
from scellrun.defaults import ScrnaQCThresholds, ScrnaQCThresholdsSN

scrna_qc = ScrnaQCThresholds(
    species="human",
    max_pct_mt=15.0,    # or whatever your tissue tolerates
    max_pct_hb=5.0,
)
"""<one-sentence rationale that gets rendered into the QC report>"""

snrna_qc = ScrnaQCThresholdsSN(species="human")  # snRNA defaults
```

Then register the name in `src/scellrun/profiles/__init__.py`:

```python
_REGISTRY: list[str] = ["default", "joint-disease", "<your_tissue>"]
```

External names use dashes (`--profile your-tissue`); module names use
underscores (`your_tissue.py`). The loader handles the swap.

## 2. Add a marker panel

Marker panels live in the same profile module. Each panel is a `dict`
of `{label: [marker genes]}`. Panel keys with biological order
(progenitor → effector → terminal) help downstream readers; the
annotate stage doesn't enforce ordering but it shows up in the panel
table.

```python
my_panel: dict[str, list[str]] = {
    "Cell type A": ["GENE1", "GENE2", "GENE3"],
    "Cell type B": ["GENE4", "GENE5"],
    ...
}
```

Two conventions matter:

- A profile with multiple panels needs a default-pick policy. The
  `joint-disease` profile auto-picks between `chondrocyte_markers` and
  `celltype_broad` based on the data's top-marker overlap (see
  `_select_panel` in `src/scellrun/scrna/annotate.py`). If your panel
  set has a single dominant panel, set `default_panel = "<name>"` at
  module level and the loader will pick it without auto-detection.
- Gene symbols are matched case-insensitively but are otherwise
  literal. If your panel uses HGNC symbols and the data uses Ensembl
  IDs, the deterministic match will return `Unassigned` for every
  cluster — there is no symbol-to-ID resolver in scellrun. Convert
  your data on the way in.

The annotate stage's panel-margin self-check is in
`src/scellrun/self_check.py::annotate_panel_tissue_mismatch`. It
fires when every cluster's best panel margin is below `0.05`,
suggesting `--panel celltype_broad`. If you ship a profile with both
chondrocyte-style and broad panels, that suggestion path will work
out of the box.

## 3. Add a stage

A stage is a CLI subcommand under `scellrun scrna *` (or a new
top-level group), plus an entry in the `analyze` orchestrator. The
canonical example is `src/scellrun/scrna/markers.py` — ~300 lines,
one `run_markers()` function, one `write_markers_artifacts()`
function, no decision-log calls inside the stage (the orchestrator
records them).

The non-negotiable contract for any new stage:

1. Stage output goes under `<run_dir>/NN_<stage>/`. Use
   `scellrun.runlayout.stage_dir(run_dir, "<name>", force=force,
   attempt_id=attempt_id)` to compute the path; do NOT hard-code.
2. The stage writes its own `report.html` (jinja2 template under
   `src/scellrun/templates/<stage>.html.j2`).
3. Every non-trivial parameter the stage applied is recorded via
   `scellrun.decisions.record(run_dir, Decision(...))`. The orchestrator
   in `src/scellrun/analyze.py` is the canonical example of how to wrap
   a stage call with decision-logging.
4. If the stage has a self-check, add it to
   `src/scellrun/self_check.py` next to `qc_self_check` and friends —
   one `def <stage>_self_check(...) -> list[SelfCheckFinding]`. The
   orchestrator picks up self-check findings, logs every fired one as a
   `<stage>.skipped_findings` row, and applies the first actionable
   `fix_payload` once when `--auto-fix` is on.
5. Add the stage to the `analyze` pipeline in
   `src/scellrun/analyze.py::run_analyze`. Slot it between the
   existing stages with the right run-meta call and the right
   `[N/5] <stage>: <one-line summary>` print.
6. Add a smoke test under `tests/test_<stage>.py` that runs against a
   tiny synthetic h5ad fixture (see `tests/conftest.py` for the
   fixture pattern — gamma-distributed gene means dodge the seurat
   HVG NaN that bit us in v0.6).

## Local dev loop

```bash
git clone https://github.com/drstrangerujn/scellrun.git
cd scellrun
conda create -n scellrun-dev python=3.11 -y && conda activate scellrun-dev
pip install -e ".[dev]"
ruff check src tests
pytest -q
```

CI runs the same `ruff check` + `pytest -q`; local-green = CI-green.

## What we don't accept

- Stages that don't write a `report.html`. Reports are the
  deliverable; everything else is intermediate.
- Defaults you can't justify in one sentence linking to a paper or
  in-house protocol. A new threshold needs a rationale string for
  the report.
- Mocking real data in tests. Use small fixtures from `conftest.py`;
  bigger fixtures committed under `tests/data/` are fine if they
  fit in <2 MB.
- Wrapping `scanpy` / `scrublet` / `harmonypy` calls behind another
  abstraction layer. The whole point of scellrun is the thin shim
  with thick defaults — every stage is "call this scanpy function
  with these parameters, log the parameters, render the report."

## Style

- Type hints everywhere. CI runs `ruff` with `E,F,I,B,UP` selected.
- Profiles docstring with the working-practice citation: a paper, an
  in-house Rmd path, or a PI-memo identifier.
- No emojis in code, comments, or report templates.
- One-sentence rationale for every new threshold, in a docstring on
  the line after the assignment (jinja2 reads it at report time).
