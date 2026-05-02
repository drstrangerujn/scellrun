"""
v1.3.0 Part A: `scellrun review <run-dir>` — Flask local server + browser
UI for human-in-the-loop overrides.

Surface:

    scellrun review <run-dir>           # serves http://127.0.0.1:8788
    scellrun review <run-dir> --port 9000
    scellrun review <run-dir> --read-only
    scellrun review <run-dir> --lang en|zh

The user opens the URL, edits cluster labels / cell exclusions /
threshold tweaks / notes, hits Save. The server writes the override
file to ``<run-dir>/06_views/review_overrides.json`` and keeps running
until Ctrl-C. The orchestrator's ``analyze --apply-overrides <path>``
flag (see scellrun.analyze) consumes the file on the next run.

Override JSON schema (schema_version=1):

    {
      "schema_version": 1,
      "reviewer": "<git config user.name or empty>",
      "saved_at": "<ISO8601>",
      "run_dir": "<absolute path>",
      "cluster_label_overrides": {"5": "Pericyte", ...},
      "cell_exclusions": ["AAACCTG-1-1", ...],
      "threshold_overrides": {"max_pct_mt": 18, ...},
      "notes": "..."
    }

Hard rules:
- Bind 127.0.0.1 only (never 0.0.0.0). The CLI does not expose a host
  flag — sharing the report URL is a tunnel-level decision.
- No symlinks, no template-dir spawn — HTML/JS lives inline in this
  module so the dependency footprint stays at one (Flask).
- Single-file aggregator: keep this module under 200 LOC of business
  logic (helpers + route handlers). When in doubt, drop a feature
  rather than splitting into a package.
"""
from __future__ import annotations

import csv
import json
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scellrun.runlayout import STAGE_DIRS

OVERRIDES_FILENAME = "review_overrides.json"
OVERRIDES_SCHEMA_VERSION = 1


def overrides_path(run_dir: Path) -> Path:
    """Canonical path the review server writes to."""
    return Path(run_dir) / STAGE_DIRS["views"] / OVERRIDES_FILENAME


def _git_user_name() -> str:
    """Best-effort `git config user.name`. Returns empty on failure."""
    try:
        out = subprocess.run(
            ["git", "config", "--get", "user.name"],
            capture_output=True,
            text=True,
            timeout=2,
        )
        return out.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return ""


def _read_annotations(run_dir: Path) -> list[dict[str, str]]:
    """Read 04_annotate/annotations.csv; return [] if absent."""
    p = run_dir / STAGE_DIRS["annotate"] / "annotations.csv"
    if not p.exists():
        return []
    rows: list[dict[str, str]] = []
    with p.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(dict(r))
    return rows


def _read_qc_thresholds(run_dir: Path) -> dict[str, float | int | None]:
    """
    Pull effective QC thresholds from 00_decisions.jsonl. The QC stage
    records max_pct_mt / max_genes / min_counts on every run; we surface
    those three so the slider has a live default.
    """
    p = run_dir / "00_decisions.jsonl"
    out: dict[str, float | int | None] = {
        "max_pct_mt": None,
        "max_genes": None,
        "min_counts": None,
    }
    if not p.exists():
        return out
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("stage") != "qc":
                continue
            key = row.get("key")
            if key in out:
                v = row.get("value")
                # Most-recent wins (later QC re-runs supersede earlier).
                out[key] = v
    return out


def _read_existing_overrides(run_dir: Path) -> dict[str, Any]:
    """Pre-populate the form with whatever was last saved (if anything)."""
    p = overrides_path(run_dir)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _detect_lang_default(run_dir: Path) -> str:
    """
    Pick lang from the first stage's recorded `lang` param. Defaults to
    'en' when 00_run.json has nothing useful.
    """
    p = run_dir / "00_run.json"
    if not p.exists():
        return "en"
    try:
        meta = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "en"
    for stage in meta.get("stages", []):
        params = stage.get("params") or {}
        v = params.get("lang")
        if v in ("en", "zh"):
            return v
    return "en"


def _parse_cell_exclusions(text: str) -> list[str]:
    """Comma- or newline-separated barcodes → deduped ordered list."""
    if not text:
        return []
    raw = text.replace("\r", "\n")
    chunks: list[str] = []
    for line in raw.split("\n"):
        for piece in line.split(","):
            piece = piece.strip()
            if piece:
                chunks.append(piece)
    seen: set[str] = set()
    out: list[str] = []
    for c in chunks:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _free_port(start: int = 8788) -> int:
    """Find a free TCP port on 127.0.0.1, starting from `start`."""
    for port in range(start, start + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise RuntimeError(f"no free port found in [{start}, {start + 100})")


_HTML_EN = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>scellrun review — {run_name}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1080px; margin: 2rem auto; padding: 0 1rem; color: #222;
         line-height: 1.55; }}
  h1 {{ font-size: 1.6rem; border-bottom: 2px solid #333; padding-bottom: .3rem; }}
  h2 {{ font-size: 1.1rem; margin-top: 1.6rem; padding-bottom: .15rem;
       border-bottom: 1px solid #eee; }}
  table.plain {{ font-size: .9rem; border-collapse: collapse; margin: .4rem 0;
                width: 100%; }}
  table.plain th, table.plain td {{ border: 1px solid #ddd; padding: .25rem .5rem;
                                   text-align: left; vertical-align: top; }}
  table.plain th {{ background: #f6f6f6; }}
  input[type="text"] {{ width: 16rem; }}
  textarea {{ width: 100%; font-family: ui-monospace, Menlo, monospace;
             font-size: .9rem; }}
  input[type="range"] {{ width: 18rem; vertical-align: middle; }}
  button {{ padding: .5rem 1rem; font-size: 1rem; cursor: pointer; }}
  .status {{ font-size: .9rem; color: #1f4e8a; margin-left: 1rem; }}
  .footer {{ font-size: .85rem; color: #666; margin-top: 1.5rem; }}
  .preview {{ font-size: .85rem; color: #555; }}
</style>
</head>
<body>
<h1>scellrun review &mdash; {run_name}</h1>
<p>Reviewer: <strong>{reviewer}</strong>{readonly_label}</p>
<p>Edit overrides below; press Save when done. Saved file:
<code>{overrides_rel}</code>.</p>

<h2>Cluster labels</h2>
{cluster_table_html}

<h2>Cell exclusions</h2>
<p>Comma- or newline-separated barcodes. <span class="preview" id="excl-preview">{excl_preview}</span></p>
<textarea id="cell_exclusions" rows="6">{cell_exclusions_text}</textarea>

<h2>Threshold overrides</h2>
<p>Sliders set intent only; nothing re-runs until you invoke
<code>scellrun analyze --apply-overrides &lt;file&gt;</code>.</p>
<table class="plain">
<tbody>
<tr><th>max_pct_mt</th>
<td><input type="range" id="max_pct_mt" min="0" max="50" step="1" value="{max_pct_mt}">
<output id="max_pct_mt_out">{max_pct_mt}</output></td></tr>
<tr><th>max_genes</th>
<td><input type="range" id="max_genes" min="500" max="10000" step="100" value="{max_genes}">
<output id="max_genes_out">{max_genes}</output></td></tr>
<tr><th>min_counts</th>
<td><input type="range" id="min_counts" min="0" max="5000" step="50" value="{min_counts}">
<output id="min_counts_out">{min_counts}</output></td></tr>
</tbody>
</table>

<h2>Notes</h2>
<textarea id="notes" rows="4" placeholder="Reviewer notes (free text)">{notes}</textarea>

<p><button id="save_btn"{readonly_attr}>Save</button>
<span class="status" id="status"></span></p>

<p class="footer">scellrun review server &middot; bound to 127.0.0.1:{port} &middot;
Ctrl-C in the terminal to stop.</p>

<script>
const readOnly = {read_only_js};
const slider_ids = ["max_pct_mt", "max_genes", "min_counts"];
slider_ids.forEach(id => {{
  const el = document.getElementById(id);
  const out = document.getElementById(id + "_out");
  el.addEventListener("input", () => {{ out.textContent = el.value; }});
}});
const exclEl = document.getElementById("cell_exclusions");
const exclPreview = document.getElementById("excl-preview");
function updatePreview() {{
  const text = exclEl.value;
  const items = text.split(/[,\\n]/).map(s => s.trim()).filter(Boolean);
  let preview = "count: " + items.length;
  if (items.length > 0) {{
    const head = items.slice(0, 5).join(", ");
    const tail = items.length > 5 ? " ... " + items.slice(-5).join(", ") : "";
    preview += " | " + head + tail;
  }}
  exclPreview.textContent = preview;
}}
exclEl.addEventListener("input", updatePreview);
updatePreview();

if (!readOnly) {{
  document.getElementById("save_btn").addEventListener("click", async () => {{
    const labels = {{}};
    document.querySelectorAll('input[data-cluster]').forEach(inp => {{
      const c = inp.getAttribute("data-cluster");
      const v = inp.value.trim();
      if (v) labels[c] = v;
    }});
    const payload = {{
      cluster_label_overrides: labels,
      cell_exclusions: exclEl.value,
      threshold_overrides: {{
        max_pct_mt: parseFloat(document.getElementById("max_pct_mt").value),
        max_genes: parseInt(document.getElementById("max_genes").value, 10),
        min_counts: parseInt(document.getElementById("min_counts").value, 10),
      }},
      notes: document.getElementById("notes").value,
    }};
    const r = await fetch("/save", {{
      method: "POST",
      headers: {{"Content-Type": "application/json"}},
      body: JSON.stringify(payload),
    }});
    const data = await r.json();
    const status = document.getElementById("status");
    if (data.ok) {{
      const t = new Date().toLocaleTimeString();
      status.textContent = "saved at " + t;
    }} else {{
      status.textContent = "save failed: " + (data.error || "unknown");
    }}
  }});
}}
</script>
</body>
</html>
"""

_HTML_ZH = """<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8">
<title>scellrun 复审 — {run_name}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "PingFang SC", sans-serif;
         max-width: 1080px; margin: 2rem auto; padding: 0 1rem; color: #222;
         line-height: 1.55; }}
  h1 {{ font-size: 1.6rem; border-bottom: 2px solid #333; padding-bottom: .3rem; }}
  h2 {{ font-size: 1.1rem; margin-top: 1.6rem; padding-bottom: .15rem;
       border-bottom: 1px solid #eee; }}
  table.plain {{ font-size: .9rem; border-collapse: collapse; margin: .4rem 0;
                width: 100%; }}
  table.plain th, table.plain td {{ border: 1px solid #ddd; padding: .25rem .5rem;
                                   text-align: left; vertical-align: top; }}
  table.plain th {{ background: #f6f6f6; }}
  input[type="text"] {{ width: 16rem; }}
  textarea {{ width: 100%; font-family: ui-monospace, Menlo, monospace;
             font-size: .9rem; }}
  input[type="range"] {{ width: 18rem; vertical-align: middle; }}
  button {{ padding: .5rem 1rem; font-size: 1rem; cursor: pointer; }}
  .status {{ font-size: .9rem; color: #1f4e8a; margin-left: 1rem; }}
  .footer {{ font-size: .85rem; color: #666; margin-top: 1.5rem; }}
  .preview {{ font-size: .85rem; color: #555; }}
</style>
</head>
<body>
<h1>scellrun 复审 &mdash; {run_name}</h1>
<p>审阅者：<strong>{reviewer}</strong>{readonly_label}</p>
<p>下面调整后点保存。覆写文件：<code>{overrides_rel}</code>。</p>

<h2>聚类标签</h2>
{cluster_table_html}

<h2>排除细胞</h2>
<p>逗号或换行分隔条码。<span class="preview" id="excl-preview">{excl_preview}</span></p>
<textarea id="cell_exclusions" rows="6">{cell_exclusions_text}</textarea>

<h2>阈值覆写</h2>
<p>滑块只记录意图；要重跑须执行
<code>scellrun analyze --apply-overrides &lt;file&gt;</code>。</p>
<table class="plain">
<tbody>
<tr><th>max_pct_mt</th>
<td><input type="range" id="max_pct_mt" min="0" max="50" step="1" value="{max_pct_mt}">
<output id="max_pct_mt_out">{max_pct_mt}</output></td></tr>
<tr><th>max_genes</th>
<td><input type="range" id="max_genes" min="500" max="10000" step="100" value="{max_genes}">
<output id="max_genes_out">{max_genes}</output></td></tr>
<tr><th>min_counts</th>
<td><input type="range" id="min_counts" min="0" max="5000" step="50" value="{min_counts}">
<output id="min_counts_out">{min_counts}</output></td></tr>
</tbody>
</table>

<h2>备注</h2>
<textarea id="notes" rows="4" placeholder="审阅者备注">{notes}</textarea>

<p><button id="save_btn"{readonly_attr}>保存</button>
<span class="status" id="status"></span></p>

<p class="footer">scellrun review 服务 &middot; 绑定 127.0.0.1:{port} &middot;
终端 Ctrl-C 退出。</p>

<script>
const readOnly = {read_only_js};
const slider_ids = ["max_pct_mt", "max_genes", "min_counts"];
slider_ids.forEach(id => {{
  const el = document.getElementById(id);
  const out = document.getElementById(id + "_out");
  el.addEventListener("input", () => {{ out.textContent = el.value; }});
}});
const exclEl = document.getElementById("cell_exclusions");
const exclPreview = document.getElementById("excl-preview");
function updatePreview() {{
  const text = exclEl.value;
  const items = text.split(/[,\\n]/).map(s => s.trim()).filter(Boolean);
  let preview = "数量: " + items.length;
  if (items.length > 0) {{
    const head = items.slice(0, 5).join(", ");
    const tail = items.length > 5 ? " ... " + items.slice(-5).join(", ") : "";
    preview += " | " + head + tail;
  }}
  exclPreview.textContent = preview;
}}
exclEl.addEventListener("input", updatePreview);
updatePreview();

if (!readOnly) {{
  document.getElementById("save_btn").addEventListener("click", async () => {{
    const labels = {{}};
    document.querySelectorAll('input[data-cluster]').forEach(inp => {{
      const c = inp.getAttribute("data-cluster");
      const v = inp.value.trim();
      if (v) labels[c] = v;
    }});
    const payload = {{
      cluster_label_overrides: labels,
      cell_exclusions: exclEl.value,
      threshold_overrides: {{
        max_pct_mt: parseFloat(document.getElementById("max_pct_mt").value),
        max_genes: parseInt(document.getElementById("max_genes").value, 10),
        min_counts: parseInt(document.getElementById("min_counts").value, 10),
      }},
      notes: document.getElementById("notes").value,
    }};
    const r = await fetch("/save", {{
      method: "POST",
      headers: {{"Content-Type": "application/json"}},
      body: JSON.stringify(payload),
    }});
    const data = await r.json();
    const status = document.getElementById("status");
    if (data.ok) {{
      const t = new Date().toLocaleTimeString();
      status.textContent = "已保存 " + t;
    }} else {{
      status.textContent = "保存失败: " + (data.error || "unknown");
    }}
  }});
}}
</script>
</body>
</html>
"""


def _render_cluster_table(
    annotations: list[dict[str, str]],
    saved_labels: dict[str, str],
    *,
    lang: str,
    read_only: bool,
) -> str:
    from html import escape as _e

    if not annotations:
        if lang == "zh":
            return "<p><em>没找到 04_annotate/annotations.csv，先跑一次 analyze。</em></p>"
        return "<p><em>No 04_annotate/annotations.csv found; run analyze first.</em></p>"

    if lang == "zh":
        head = "<thead><tr><th>聚类</th><th>当前标签</th><th>新标签</th></tr></thead>"
    else:
        head = "<thead><tr><th>cluster</th><th>current label</th><th>new label</th></tr></thead>"

    rows: list[str] = []
    disabled = " disabled" if read_only else ""
    for r in annotations:
        cid = r.get("cluster", "")
        cur = r.get("panel_label", "")
        prev = saved_labels.get(cid, "")
        rows.append(
            "<tr>"
            f"<td>{_e(cid)}</td>"
            f"<td>{_e(cur)}</td>"
            f'<td><input type="text" data-cluster="{_e(cid)}" '
            f'value="{_e(prev)}"{disabled}></td>'
            "</tr>"
        )
    return f'<table class="plain">{head}<tbody>{"".join(rows)}</tbody></table>'


def _excl_preview(items: list[str], lang: str) -> str:
    if lang == "zh":
        prefix = "数量: "
    else:
        prefix = "count: "
    if not items:
        return prefix + "0"
    head = ", ".join(items[:5])
    tail = ""
    if len(items) > 5:
        tail = " ... " + ", ".join(items[-5:])
    return f"{prefix}{len(items)} | {head}{tail}"


def create_app(
    run_dir: Path,
    *,
    read_only: bool = False,
    lang: str = "en",
):
    """
    Build a Flask app for the given run-dir.

    Lazy import of Flask: importing review.py at CLI startup mustn't pay
    the Flask import cost on commands that never use the review server.
    """
    from flask import Flask, jsonify, request

    run_dir = Path(run_dir).resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"run-dir does not exist: {run_dir}")

    app = Flask(__name__)
    # Disable Flask's own log noise on stdout — the CLI prints its own
    # "Serving at..." message and we don't want werkzeug duplicating.
    import logging as _logging

    _logging.getLogger("werkzeug").setLevel(_logging.ERROR)

    template = _HTML_ZH if lang == "zh" else _HTML_EN

    @app.route("/", methods=["GET"])
    def index() -> str:
        annotations = _read_annotations(run_dir)
        thresholds = _read_qc_thresholds(run_dir)
        existing = _read_existing_overrides(run_dir)
        saved_labels = existing.get("cluster_label_overrides") or {}
        existing_excl: list[str] = list(existing.get("cell_exclusions") or [])
        existing_thresh = existing.get("threshold_overrides") or {}
        notes = existing.get("notes") or ""

        cluster_html = _render_cluster_table(
            annotations, saved_labels, lang=lang, read_only=read_only
        )

        from html import escape as _e

        # Threshold initial value: prefer existing override, fall back to
        # what the QC stage actually applied, and finally to a sane const.
        max_pct_mt = existing_thresh.get("max_pct_mt", thresholds.get("max_pct_mt") or 20)
        max_genes = existing_thresh.get("max_genes", thresholds.get("max_genes") or 4000)
        min_counts = existing_thresh.get("min_counts", thresholds.get("min_counts") or 500)

        readonly_label = (
            (" — 只读" if lang == "zh" else " — read-only")
            if read_only
            else ""
        )
        readonly_attr = " disabled" if read_only else ""
        return template.format(
            run_name=_e(run_dir.name),
            reviewer=_e(_git_user_name() or ("(未设置)" if lang == "zh" else "(unset)")),
            readonly_label=readonly_label,
            overrides_rel=_e(str(overrides_path(run_dir).relative_to(run_dir))),
            cluster_table_html=cluster_html,
            cell_exclusions_text=_e("\n".join(existing_excl)),
            excl_preview=_excl_preview(existing_excl, lang),
            max_pct_mt=_e(str(max_pct_mt)),
            max_genes=_e(str(max_genes)),
            min_counts=_e(str(min_counts)),
            notes=_e(str(notes)),
            readonly_attr=readonly_attr,
            read_only_js="true" if read_only else "false",
            port=app.config.get("SCELLRUN_PORT", 8788),
        )

    @app.route("/save", methods=["POST"])
    def save() -> Any:
        if read_only:
            return jsonify({"ok": False, "error": "server is read-only"}), 403
        payload = request.get_json(silent=True) or {}
        labels = payload.get("cluster_label_overrides") or {}
        if not isinstance(labels, dict):
            return jsonify({"ok": False, "error": "cluster_label_overrides must be a dict"}), 400
        excl_raw = payload.get("cell_exclusions") or ""
        if isinstance(excl_raw, list):
            excl_list = [str(x) for x in excl_raw if str(x).strip()]
        else:
            excl_list = _parse_cell_exclusions(str(excl_raw))
        thresholds = payload.get("threshold_overrides") or {}
        notes = str(payload.get("notes") or "")

        out = {
            "schema_version": OVERRIDES_SCHEMA_VERSION,
            "reviewer": _git_user_name(),
            "saved_at": datetime.now(timezone.utc).isoformat(),
            "run_dir": str(run_dir),
            "cluster_label_overrides": {str(k): str(v) for k, v in labels.items() if str(v).strip()},
            "cell_exclusions": excl_list,
            "threshold_overrides": {
                k: thresholds[k]
                for k in ("max_pct_mt", "max_genes", "min_counts")
                if k in thresholds and thresholds[k] is not None
            },
            "notes": notes,
        }
        p = overrides_path(run_dir)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        return jsonify({"ok": True, "path": str(p), "saved_at": out["saved_at"]})

    return app


def serve(
    run_dir: Path,
    *,
    port: int = 8788,
    read_only: bool = False,
    lang: str | None = None,
) -> None:
    """
    Start the review server (blocking). Imported by the CLI command.

    Picks a free port if `port` is taken (scans 100 above). Lang
    defaults to whatever the run-dir's first stage used, falling back
    to 'en'.
    """
    if lang is None:
        lang = _detect_lang_default(run_dir)
    actual_port = port
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", port))
    except OSError:
        actual_port = _free_port(start=port + 1)

    app = create_app(run_dir, read_only=read_only, lang=lang)
    app.config["SCELLRUN_PORT"] = actual_port
    print(f"scellrun review serving at http://127.0.0.1:{actual_port}/  (Ctrl-C to stop)")
    app.run(host="127.0.0.1", port=actual_port, debug=False, use_reloader=False)


def make_review_overrides_json(
    run_dir: Path,
    *,
    cluster_label_overrides: dict[str, str] | None = None,
    cell_exclusions: list[str] | None = None,
    threshold_overrides: dict[str, float | int] | None = None,
    notes: str = "",
    reviewer: str = "",
    out_path: Path | None = None,
) -> Path:
    """
    Test helper: write a schema-conformant override file without spinning
    up the Flask app. Tests + the CLI's `--apply-overrides` consumer
    both round-trip through this shape.
    """
    out = {
        "schema_version": OVERRIDES_SCHEMA_VERSION,
        "reviewer": reviewer,
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(Path(run_dir).resolve()),
        "cluster_label_overrides": dict(cluster_label_overrides or {}),
        "cell_exclusions": list(cell_exclusions or []),
        "threshold_overrides": dict(threshold_overrides or {}),
        "notes": notes,
    }
    if out_path is None:
        out_path = overrides_path(Path(run_dir))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    return out_path


__all__ = [
    "OVERRIDES_FILENAME",
    "OVERRIDES_SCHEMA_VERSION",
    "create_app",
    "make_review_overrides_json",
    "overrides_path",
    "serve",
]
