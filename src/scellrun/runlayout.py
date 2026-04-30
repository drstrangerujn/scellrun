"""
Run-directory layout — single source of truth for where stage outputs go.

Every stage writes into `<run_dir>/<NN_stage>/`. This makes
`scellrun run --steps qc,integrate,markers,annotate` (v0.4) just an
orchestration shell over the same single-stage commands.

Layout:
    scellrun_out/run-20260430-094530/
        00_run.json              run metadata: profile, species, command, timestamps
        01_qc/
            report.html
            per_cell_metrics.csv
            qc.h5ad              annotated AnnData; canonical handoff to v0.2 integrate
        02_integrate/            (v0.2)
        03_markers/              (v0.3)
        04_annotate/             (v0.4)
        05_report/               (v0.5)
"""
from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path

STAGE_DIRS = {
    "qc": "01_qc",
    "integrate": "02_integrate",
    "markers": "03_markers",
    "annotate": "04_annotate",
    "report": "05_report",
}


def default_run_dir() -> Path:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    return Path("scellrun_out") / f"run-{ts}"


class StageOutputExists(FileExistsError):
    """Raised when a stage subdir already contains artifacts (and --force not set)."""


def stage_dir(run_dir: Path, stage: str, *, force: bool = False) -> Path:
    """
    Resolve and create the stage subdirectory.

    If artifacts already exist there (any non-hidden file) and `force` is
    False, raise `StageOutputExists`. This guards the silent-overwrite trap
    where re-running a stage into the same run-dir leaves the manifest
    showing two runs but only the latest payload on disk.
    """
    if stage not in STAGE_DIRS:
        raise ValueError(f"unknown stage {stage!r}; expected one of {list(STAGE_DIRS)}")
    p = run_dir / STAGE_DIRS[stage]
    if p.exists():
        existing = [f for f in p.iterdir() if not f.name.startswith(".")]
        if existing and not force:
            raise StageOutputExists(
                f"{p} already has artifacts: {[f.name for f in existing]}. "
                "Re-run with --force to overwrite, or pick a fresh --run-dir."
            )
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_run_meta(run_dir: Path, command: str, params: dict) -> None:
    """Write/append a small JSON manifest so any stage can be reconstructed later."""
    run_dir.mkdir(parents=True, exist_ok=True)
    meta_path = run_dir / "00_run.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
    else:
        meta = {"created_at": datetime.now(timezone.utc).isoformat(), "stages": []}

    safe_params = {k: _jsonable(v) for k, v in params.items()}
    meta["stages"].append({
        "command": command,
        "ran_at": datetime.now(timezone.utc).isoformat(),
        "params": safe_params,
    })
    meta_path.write_text(json.dumps(meta, indent=2, default=str))


def _jsonable(v):
    if is_dataclass(v):
        return asdict(v)
    if isinstance(v, Path):
        return str(v)
    return v
