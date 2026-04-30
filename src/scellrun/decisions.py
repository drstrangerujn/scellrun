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

Concurrency: appends use ``fcntl.flock`` (LOCK_EX) on POSIX so the v0.8
auto-fix retry path and any future multi-process callers don't interleave
JSONL lines. Windows has no ``fcntl`` — the lock is a no-op there; the
single-writer assumption that held in v0.7 still holds for the Windows
path.

Fields:
    schema_version: shape version of this Decision row; bumped only on
                    incompatible changes. v0.9.1 = 1.
    stage:     pipeline stage that made the call (qc / integrate /
               markers / annotate / analyze)
    key:       short slug for the decision (e.g. "max_pct_mt", "sample_key")
    value:     the value chosen (anything JSON-serialisable)
    default:   the default the code would have used absent input. For a
               profile-influenced threshold this is the profile-applied
               baseline (e.g. joint-disease's 20% mt ceiling), not the
               library default. Pre-v0.9.1 reports treated this as the
               library default — that bug is fixed in v0.9.1.
    source:    one of "user" (caller passed an override),
               "auto" (deterministic logic chose it),
               "ai"   (LLM call decided / recommended it)
    rationale: one sentence explaining the choice (cite paper / heuristic /
               table entry; this is what the user reads first)
    fix_payload: optional structured fix dict carried over from a
                 SelfCheckFinding so downstream tools can apply the fix
                 mechanically without parsing the rationale prose.
    attempt_id: UUID4 generated once per `analyze` invocation (or per
                per-stage CLI call). Lets the report group rows by retry
                attempt when --auto-fix or a manual --force re-runs a
                stage. Empty string when unknown.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DECISIONS_FILENAME = "00_decisions.jsonl"

VALID_SOURCES = ("user", "auto", "ai")

SCHEMA_VERSION = 1


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
    schema_version: int = SCHEMA_VERSION
    attempt_id: str = ""
    fix_payload: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.source not in VALID_SOURCES:
            raise ValueError(
                f"Decision.source must be one of {VALID_SOURCES}, got {self.source!r}"
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "stage": self.stage,
            "key": self.key,
            "value": _jsonable(self.value),
            "default": _jsonable(self.default),
            "source": self.source,
            "rationale": self.rationale,
            "fix_payload": _jsonable(self.fix_payload) if self.fix_payload else None,
            "attempt_id": self.attempt_id,
            "ts": self.ts,
        }

    @classmethod
    def from_choice(
        cls,
        *,
        stage: str,
        key: str,
        value: Any,
        default: Any,
        is_user_override: bool,
        rationale: str,
        attempt_id: str = "",
        fix_payload: dict[str, Any] | None = None,
    ) -> Decision:
        """
        Build a Decision with `source` derived from explicit caller intent.

        Pre-v0.9.1 callers diffed `value` vs `default` to set source — that
        false-tagged auto-resolved values whose internal pick happened to
        differ from a hard-coded default (most visibly the orchestrator's
        chosen_resolution_for_annotate). Use this constructor when the
        call site knows whether the caller explicitly supplied that
        argument; `is_user_override=False` always becomes `source="auto"`,
        regardless of whether `value` and `default` happen to match.
        """
        return cls(
            stage=stage,
            key=key,
            value=value,
            default=default,
            source="user" if is_user_override else "auto",
            rationale=rationale,
            attempt_id=attempt_id,
            fix_payload=fix_payload,
        )


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


def _locked_append(p: Path, lines: list[str]) -> None:
    """
    Append text lines to ``p`` under an exclusive advisory lock.

    On POSIX we ``fcntl.flock`` the file for the lifetime of the write;
    a competing writer (the v0.8 auto-fix retry, a parallel CLI) waits
    on the lock and then writes its own block. On Windows fcntl doesn't
    exist — we fall back to plain append. Single-writer scellrun on
    Windows still works; the multi-writer guarantee is POSIX-only.
    """
    if os.name == "nt":
        with p.open("a", encoding="utf-8") as f:
            for line in lines:
                f.write(line)
        return

    import fcntl

    with p.open("a", encoding="utf-8") as f:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            for line in lines:
                f.write(line)
            f.flush()
        finally:
            try:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass


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
    line = json.dumps(decision.to_dict(), ensure_ascii=False, default=str) + "\n"
    _locked_append(p, [line])


def record_many(run_dir: Path | None, decisions: list[Decision]) -> None:
    """Append a list of Decisions in one go (still one JSON object per line)."""
    if run_dir is None or not decisions:
        return
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    p = decisions_path(run_dir)
    lines = [
        json.dumps(d.to_dict(), ensure_ascii=False, default=str) + "\n"
        for d in decisions
    ]
    _locked_append(p, lines)


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


def truncate_stage(run_dir: Path | None, stage: str) -> int:
    """
    Drop every JSONL row whose ``stage`` matches and rewrite the file.

    Returns the count of rows dropped. Used by stage_dir(force=True) to
    keep the decision log honest when the user re-runs a stage — the
    prior attempt's rows get cleared so the re-run's rows replace them
    (rather than accumulating duplicate `qc.profile`, `qc.max_pct_mt`,
    etc. on every retry).

    A None ``run_dir`` is a no-op (returns 0); a missing file likewise.
    """
    if run_dir is None:
        return 0
    p = decisions_path(Path(run_dir))
    if not p.exists():
        return 0

    if os.name == "nt":
        rows = p.read_text(encoding="utf-8").splitlines()
        kept: list[str] = []
        n_dropped = 0
        for line in rows:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                kept.append(line)
                continue
            if obj.get("stage") == stage:
                n_dropped += 1
                continue
            kept.append(line)
        p.write_text("\n".join(kept) + ("\n" if kept else ""), encoding="utf-8")
        return n_dropped

    import fcntl

    with p.open("r+", encoding="utf-8") as f:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            f.seek(0)
            rows = f.read().splitlines()
            kept: list[str] = []
            n_dropped = 0
            for line in rows:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    kept.append(line)
                    continue
                if obj.get("stage") == stage:
                    n_dropped += 1
                    continue
                kept.append(line)
            f.seek(0)
            f.truncate()
            if kept:
                f.write("\n".join(kept) + "\n")
            f.flush()
        finally:
            try:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
    return n_dropped


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
