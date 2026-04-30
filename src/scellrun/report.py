"""
Multi-stage report aggregator: stitch QC + integrate + markers + annotate
into one HTML deliverable with the three-tier provenance trail
(data / inference / literature).

This is intentionally an HTML aggregator, not PDF. Reasons:
- HTML preserves links and tables; PDF flattens them.
- Browsers print to PDF on demand; one less LaTeX dep in the stack.
- Future v0.6 can add WeasyPrint for true PDF if needed.

The aggregator reads each stage's report.html (and supporting CSV/PNG)
from <run_dir>/<NN_stage>/, plus the 00_run.json manifest, and emits
<run_dir>/05_report/index.html that links + summarizes.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from scellrun.decisions import group_by_stage, read_decisions
from scellrun.runlayout import STAGE_DIRS

STAGE_ORDER = ["qc", "integrate", "markers", "annotate"]
# Stage order for the decision-summary section. "analyze" goes last so the
# orchestrator's own picks (e.g. chosen_resolution_for_annotate) read after
# the per-stage decisions they depend on.
DECISION_STAGE_ORDER = ["qc", "integrate", "markers", "annotate", "analyze"]


@dataclass
class StageInfo:
    name: str
    dir_name: str
    present: bool
    report_path: Path | None
    artifacts: dict[str, Path]


def _scan_stage(run_dir: Path, stage: str) -> StageInfo:
    sub = run_dir / STAGE_DIRS[stage]
    if not sub.exists():
        return StageInfo(name=stage, dir_name=STAGE_DIRS[stage], present=False, report_path=None, artifacts={})

    artifacts: dict[str, Path] = {}
    for f in sorted(sub.iterdir()):
        if f.is_file() and not f.name.startswith("."):
            artifacts[f.name] = f
    report = sub / "report.html"
    return StageInfo(
        name=stage,
        dir_name=STAGE_DIRS[stage],
        present=True,
        report_path=report if report.exists() else None,
        artifacts=artifacts,
    )


def _read_manifest(run_dir: Path) -> dict:
    p = run_dir / "00_run.json"
    if not p.exists():
        return {"stages": []}
    return json.loads(p.read_text())


def _glance_from_decisions(decisions: list[dict]) -> dict[str, object]:
    """
    Scrape the decision log for the at-a-glance numbers the index page
    surfaces. Best-effort: if a key isn't present we leave the entry None
    rather than crashing the report build.
    """
    glance: dict[str, object] = {
        "chosen_resolution": None,
        "qc_pass_rate_pct": None,
        "qc_n_pass": None,
        "qc_n_in": None,
        "method": None,
        "panel": None,
        "self_check_codes": [],
    }
    for d in decisions:
        stage = d.get("stage")
        key = d.get("key")
        if stage == "analyze" and key == "chosen_resolution_for_annotate":
            glance["chosen_resolution"] = d.get("value")
        elif stage == "integrate" and key == "method":
            glance["method"] = d.get("value")
        elif stage == "annotate" and key == "panel":
            glance["panel"] = d.get("value")
        elif stage == "analyze" and key == "auto_fix.qc.outcome":
            try:
                v = str(d.get("value", "")).rstrip("%")
                glance["qc_pass_rate_pct"] = float(v)
            except (TypeError, ValueError):
                pass
        if isinstance(key, str) and key.startswith("self_check.") and key.endswith(".trigger"):
            code = d.get("value")
            if isinstance(code, str):
                glance["self_check_codes"].append(code)  # type: ignore[union-attr]
    return glance


def _glance_from_artifacts(run_dir: Path) -> dict[str, object]:
    """
    Scrape stage artifact files for the at-a-glance numbers.

    We deliberately prefer stage-artifact files over decision-log scrapes
    here because the artifacts are the source of truth — the decision log
    only captures parameters, not computed results like QC pass-rate or
    cluster counts. Returns None entries where the artifact is missing.
    """
    out: dict[str, object] = {
        "qc_pass_rate_pct": None,
        "qc_n_pass": None,
        "qc_n_in": None,
        "n_clusters": None,
        "top_celltypes": [],  # list of (label, n_cells)
    }

    # QC pass-rate from per_cell_metrics.csv if present.
    qc_csv = run_dir / STAGE_DIRS["qc"] / "per_cell_metrics.csv"
    if qc_csv.exists():
        try:
            with qc_csv.open() as f:
                header = f.readline().strip().split(",")
                if "scellrun_qc_pass" in header:
                    idx = header.index("scellrun_qc_pass")
                    n_in = 0
                    n_pass = 0
                    for line in f:
                        parts = line.rstrip("\n").split(",")
                        if len(parts) <= idx:
                            continue
                        n_in += 1
                        v = parts[idx].strip()
                        if v == "True" or v == "true" or v == "1":
                            n_pass += 1
                    if n_in:
                        out["qc_n_in"] = n_in
                        out["qc_n_pass"] = n_pass
                        out["qc_pass_rate_pct"] = round(100.0 * n_pass / n_in, 1)
        except Exception:
            pass

    # Annotate: n_clusters, top three labels with cell counts.
    annot_csv = run_dir / STAGE_DIRS["annotate"] / "annotations.csv"
    if annot_csv.exists():
        try:
            with annot_csv.open() as f:
                header = f.readline().strip().split(",")
                if "panel_label" in header and "cluster" in header:
                    label_idx = header.index("panel_label")
                    n_clusters = 0
                    label_counts: dict[str, int] = {}
                    for line in f:
                        parts = _split_csv_line(line.rstrip("\n"))
                        if len(parts) <= label_idx:
                            continue
                        n_clusters += 1
                        label = parts[label_idx]
                        label_counts[label] = label_counts.get(label, 0) + 1
                    out["n_clusters"] = n_clusters
                    # Top 3 labels by cluster count (cells-per-cluster not
                    # in this file; cluster count is the next-best proxy
                    # at the index level).
                    top = sorted(label_counts.items(), key=lambda kv: -kv[1])[:3]
                    out["top_celltypes"] = [
                        {"label": lbl, "n_clusters": n}
                        for lbl, n in top
                        if lbl
                    ]
        except Exception:
            pass

    return out


def _split_csv_line(line: str) -> list[str]:
    """Tiny CSV splitter that respects double-quoted fields containing commas."""
    out: list[str] = []
    cur: list[str] = []
    in_quotes = False
    for ch in line:
        if ch == '"':
            in_quotes = not in_quotes
            continue
        if ch == "," and not in_quotes:
            out.append("".join(cur))
            cur = []
            continue
        cur.append(ch)
    out.append("".join(cur))
    return out


def build_report(
    run_dir: Path,
    out_dir: Path,
    *,
    lang: str = "en",
) -> dict[str, Path]:
    """
    Build the index.html report aggregating whatever stages are present.

    Returns dict of {name: path} for created artifacts.
    """
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)

    stages = [_scan_stage(run_dir, s) for s in STAGE_ORDER]
    manifest = _read_manifest(run_dir)

    # v0.7: aggregate decisions for the "Decision summary" section.
    raw_decisions = read_decisions(run_dir)
    grouped = group_by_stage(raw_decisions)
    # Render in canonical order; tack any unknown stage names on the end.
    decisions_by_stage: list[tuple[str, list[dict]]] = []
    seen: set[str] = set()
    for s in DECISION_STAGE_ORDER:
        if s in grouped:
            decisions_by_stage.append((s, grouped[s]))
            seen.add(s)
    for s, items in grouped.items():
        if s not in seen:
            decisions_by_stage.append((s, items))
    n_decisions_total = sum(len(v) for _, v in decisions_by_stage)
    n_ai_decisions = sum(
        1 for _, items in decisions_by_stage for d in items if d.get("source") == "ai"
    )

    # v0.9.1 (#007): at-a-glance numbers from artifacts + decision log.
    glance = _glance_from_artifacts(run_dir)
    glance_decisions = _glance_from_decisions(raw_decisions)
    # Decision-log values fill in fields the artifact scrape didn't cover.
    for k in ("chosen_resolution", "method", "panel", "self_check_codes"):
        if k not in glance or not glance.get(k):
            glance[k] = glance_decisions.get(k)

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_full_report_zh.html.j2" if lang == "zh" else "scrna_full_report.html.j2"
    template = env.get_template(template_name)
    # Render only the run-dir basename in the report. Showing the absolute
    # path leaks the user's filesystem layout (e.g. /tmp/foo, /home/.../...)
    # into a deliverable that gets emailed / printed-to-PDF, and the path
    # stops being valid once the report moves. See ISSUES.md #008 from the
    # v0.7 OA dogfood.
    html = template.render(
        run_dir=run_dir.name,
        stages=stages,
        manifest=manifest,
        decisions_by_stage=decisions_by_stage,
        n_decisions_total=n_decisions_total,
        n_ai_decisions=n_ai_decisions,
        glance=glance,
    )
    index_path = out_dir / "index.html"
    index_path.write_text(html, encoding="utf-8")
    return {"index": index_path}
