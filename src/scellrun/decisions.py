"""
Decision log: a finer-grained sibling of `00_run.json`.

`00_run.json` (see `runlayout.write_run_meta`) captures the parameters a
stage was *called with*. The decision log captures the choices the
pipeline *made* — which thresholds were applied, what got auto-detected,
where AI was consulted — each with a human-readable one-sentence
rationale. The user reads the decision log on the report's first page
and answers "why mt% 20?" / "why res 0.5?" without opening any
sub-report.

Storage: append-only JSONL at ``<run_dir>/00_decisions.jsonl``. One line
per Decision. Stages append as they run; the report aggregator reads the
whole file and groups by stage for rendering.

Fields:
    stage:     pipeline stage that made the call (qc / integrate /
               markers / annotate / analyze)
    key:       short slug for the decision (e.g. "max_pct_mt", "sample_key")
    value:     the value chosen (anything JSON-serialisable)
    default:   the default the code would have used absent input;
               same type as `value` or None if there is no default.
               Only ``value != default`` decisions get the "user override"
               styling in the report.
    source:    one of "user" (caller passed an override),
               "auto" (deterministic logic chose it),
               "ai"   (LLM call decided / recommended it)
    rationale: one sentence explaining the choice (cite paper / heuristic /
               table entry; this is what the user reads first)
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DECISIONS_FILENAME = "00_decisions.jsonl"

VALID_SOURCES = ("user", "auto", "ai")


@dataclass
class Decision:
    """One recorded decision; serialised one-per-line into the JSONL log."""

    stage: str
    key: str
    value: Any
    default: Any = None
    source: str = "auto"
    rationale: str = ""
    ts: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def __post_init__(self) -> None:
        if self.source not in VALID_SOURCES:
            raise ValueError(
                f"Decision.source must be one of {VALID_SOURCES}, got {self.source!r}"
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "key": self.key,
            "value": _jsonable(self.value),
            "default": _jsonable(self.default),
            "source": self.source,
            "rationale": self.rationale,
            "ts": self.ts,
        }


def _jsonable(v: Any) -> Any:
    """Best-effort coerce a value into something json.dumps can handle."""
    if v is None:
        return None
    if isinstance(v, (str, int, float, bool)):
        return v
    if isinstance(v, Path):
        return str(v)
    if is_dataclass(v):
        return asdict(v)
    if isinstance(v, dict):
        return {str(k): _jsonable(val) for k, val in v.items()}
    if isinstance(v, (list, tuple, set)):
        return [_jsonable(x) for x in v]
    return str(v)


def decisions_path(run_dir: Path) -> Path:
    """Canonical path to the decision log inside a run-dir."""
    return Path(run_dir) / DECISIONS_FILENAME


def record(run_dir: Path | None, decision: Decision) -> None:
    """
    Append one Decision to ``<run_dir>/00_decisions.jsonl``.

    A None ``run_dir`` is a no-op; callers can pass it through unchanged
    when the user runs a stage without a configured run-dir (e.g. unit
    tests that don't care about the log).
    """
    if run_dir is None:
        return
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    p = decisions_path(run_dir)
    with p.open("a", encoding="utf-8") as f:
        f.write(json.dumps(decision.to_dict(), ensure_ascii=False, default=str) + "\n")


def record_many(run_dir: Path | None, decisions: list[Decision]) -> None:
    """Append a list of Decisions in one go (still one JSON object per line)."""
    if run_dir is None or not decisions:
        return
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    p = decisions_path(run_dir)
    with p.open("a", encoding="utf-8") as f:
        for d in decisions:
            f.write(json.dumps(d.to_dict(), ensure_ascii=False, default=str) + "\n")


def read_decisions(run_dir: Path) -> list[dict[str, Any]]:
    """
    Read the JSONL log and return a list of dicts in the order they were
    written. Returns an empty list if the file doesn't exist.

    Used by the report aggregator to group decisions by stage.
    """
    p = decisions_path(run_dir)
    if not p.exists():
        return []
    out: list[dict[str, Any]] = []
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                # Skip malformed lines rather than crash the report.
                continue
    return out


def group_by_stage(decisions: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """
    Bucket decisions by stage, preserving insertion order within each
    bucket. The report uses the canonical stage order (qc, integrate,
    markers, annotate, analyze) but unknown stages are kept too so a
    custom downstream stage doesn't get silently dropped.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}
    for d in decisions:
        stage = d.get("stage", "unknown")
        grouped.setdefault(stage, []).append(d)
    return grouped
