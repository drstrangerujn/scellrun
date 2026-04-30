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

from scellrun.runlayout import STAGE_DIRS

STAGE_ORDER = ["qc", "integrate", "markers", "annotate"]


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

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_full_report_zh.html.j2" if lang == "zh" else "scrna_full_report.html.j2"
    template = env.get_template(template_name)
    html = template.render(
        run_dir=str(run_dir.resolve()),
        stages=stages,
        manifest=manifest,
    )
    index_path = out_dir / "index.html"
    index_path.write_text(html, encoding="utf-8")
    return {"index": index_path}
