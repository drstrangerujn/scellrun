"""
Views layer (v1.2.0): pure-HTML cross-sections of a finished run-dir.

`05_report/index.html` aggregates the run end-to-end. The views layer
slices the same artifacts three other ways:

- ``06_views/by_resolution/<R>/index.html`` — one page per resolution,
  surfacing per-cluster top markers + (when applicable) annotations.
- ``06_views/by_cluster/cluster_<C>.html`` — one page per cluster at
  the resolution that was annotated, with that cluster's top markers,
  panel + AI rationale, and sample distribution.
- ``06_views/by_decision_source/index.html`` — three lists (ai / user /
  auto-with-non-trivial-rationale) parsed from ``00_decisions.jsonl``.

The views are pure HTML referencing stage artifacts via relative paths
(``../../../02_integrate/umap_grid.png`` etc.). No symlinks, no copies,
no hardlinks — the view directory is cheap and disposable, the source
of truth still lives in the per-stage subdirs.

Resolution naming: the URL segment matches the leiden_res_X column
suffix (`0.5` → `by_resolution/0.5/`). Cluster IDs are taken verbatim
from ``annotations.csv`` (already strings on disk).
"""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from scellrun.report import _split_csv_line  # tiny CSV splitter we reuse
from scellrun.runlayout import STAGE_DIRS

# Auto-rationale keys that aren't worth showing on the by_decision_source
# page. These are the "boring" knob defaults that don't tell a reviewer
# anything new — every run has them, and they bury the interesting picks.
_BORING_AUTO_KEYS: frozenset[str] = frozenset({
    "min_genes",
    "max_genes",
    "min_counts",
    "min_cells_per_gene",
    "max_pct_ribo",
    "flag_doublets",
    "logfc_threshold",
    "pct_min",
    "only_positive",
    "top_n",
    "n_pcs",
})


def _read_decisions_jsonl(run_dir: Path) -> list[dict]:
    """Parse 00_decisions.jsonl; return [] if absent or unreadable."""
    p = run_dir / "00_decisions.jsonl"
    if not p.exists():
        return []
    out: list[dict] = []
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def _resolutions_present(run_dir: Path) -> list[float]:
    """
    Discover resolutions that have a markers CSV on disk.
    File pattern is ``markers_res_<R>.csv`` per markers stage.
    """
    markers_dir = run_dir / STAGE_DIRS["markers"]
    if not markers_dir.exists():
        return []
    out: list[float] = []
    for p in sorted(markers_dir.glob("markers_res_*.csv")):
        stem = p.stem.removeprefix("markers_res_")
        try:
            out.append(float(stem))
        except ValueError:
            continue
    return out


def _read_markers_csv(p: Path) -> list[dict[str, str]]:
    """Read markers CSV into a list of dicts. Returns [] on missing/empty."""
    if not p.exists():
        return []
    rows: list[dict[str, str]] = []
    with p.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(dict(r))
    return rows


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


def _annotate_resolution(decisions: list[dict]) -> float | None:
    """
    Return the resolution `analyze` chose for the annotate stage, by
    scanning the decision log for the orchestrator's
    `analyze.chosen_resolution_for_annotate` row. Returns None if not
    present (e.g. the user ran per-stage commands without analyze).
    """
    for d in decisions:
        if d.get("stage") == "analyze" and d.get("key") == "chosen_resolution_for_annotate":
            v = d.get("value")
            try:
                return float(v)
            except (TypeError, ValueError):
                continue
    return None


def _e(s: object) -> str:
    """HTML-escape any value, coercing to str. Handles None as empty string."""
    if s is None:
        return ""
    return html.escape(str(s), quote=True)


# ---------- by_resolution ---------------------------------------------------


def _format_resolution_segment(res: float) -> str:
    """Render a resolution float into the URL path segment used on disk."""
    return f"{res:g}"


def _top_markers_per_cluster(
    rows: list[dict[str, str]], top_n: int = 5
) -> dict[str, list[dict[str, str]]]:
    """
    Group markers CSV rows by cluster and keep the first top_n rows
    per cluster (CSVs are written in score order by markers stage).
    """
    out: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        c = str(r.get("cluster", ""))
        if not c:
            continue
        if c not in out:
            out[c] = []
        if len(out[c]) < top_n:
            out[c].append(r)
    return out


def _markers_table_html(
    grouped: dict[str, list[dict[str, str]]],
    columns: tuple[str, ...] = ("cluster", "gene", "log2fc"),
) -> str:
    """Render a small inline table for top-N markers."""
    head = "".join(f"<th>{_e(c)}</th>" for c in columns)
    body_rows: list[str] = []
    # Sort clusters by integer if all parseable, else lexicographically.
    try:
        clusters_sorted = sorted(grouped.keys(), key=lambda c: int(c))
    except ValueError:
        clusters_sorted = sorted(grouped.keys())
    for c in clusters_sorted:
        for r in grouped[c]:
            cells = []
            for col in columns:
                # markers CSV uses 'logfoldchanges' (scanpy default) but
                # may also write 'log2fc'; surface whichever is present.
                if col == "log2fc":
                    v = r.get("log2fc") or r.get("logfoldchanges") or ""
                else:
                    v = r.get(col, "")
                cells.append(f"<td>{_e(v)}</td>")
            body_rows.append(f"<tr>{''.join(cells)}</tr>")
    if not body_rows:
        return "<p><em>No markers recorded.</em></p>"
    return (
        "<table class=\"plain\">"
        f"<thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table>"
    )


def _annotations_table_html(annotations: list[dict[str, str]]) -> str:
    """Inline table of cluster | panel_label | panel_score | margin."""
    if not annotations:
        return "<p><em>No annotations recorded.</em></p>"
    rows: list[str] = []
    for r in annotations:
        rows.append(
            "<tr>"
            f"<td>{_e(r.get('cluster', ''))}</td>"
            f"<td>{_e(r.get('panel_label', ''))}</td>"
            f"<td>{_e(r.get('panel_score', ''))}</td>"
            f"<td>{_e(r.get('panel_margin', ''))}</td>"
            "</tr>"
        )
    return (
        "<table class=\"plain\">"
        "<thead><tr>"
        "<th>cluster</th><th>panel_label</th><th>panel_score</th><th>margin</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody>"
        "</table>"
    )


def _by_resolution_html(
    res: float,
    markers_rows: list[dict[str, str]],
    annotations: list[dict[str, str]] | None,
) -> str:
    seg = _format_resolution_segment(res)
    top_grouped = _top_markers_per_cluster(markers_rows, top_n=5)
    annotations_section = ""
    if annotations:
        annotations_section = (
            "<h2>Annotations</h2>"
            f"{_annotations_table_html(annotations)}"
        )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>scellrun — resolution {_e(seg)}</title>
{_VIEWS_CSS}
</head>
<body>
<p><a href="../../index.html">&larr; views index</a></p>
<h1>Resolution {_e(seg)}</h1>

<h2>UMAP</h2>
<p>The combined per-resolution UMAP grid lives in the integrate stage:</p>
<p><img src="../../../{STAGE_DIRS['integrate']}/umap_grid.png" alt="UMAP grid (all resolutions)"></p>

<h2>Cluster sizes</h2>
<p><img src="../../../{STAGE_DIRS['integrate']}/cluster_sizes.png" alt="Cluster sizes per resolution"></p>

<h2>Top markers (top 5 per cluster)</h2>
{_markers_table_html(top_grouped)}

{annotations_section}

<p class="footer">
  <a href="../../../{STAGE_DIRS['markers']}/markers_res_{_e(seg)}.csv">Full markers CSV</a>
</p>
</body>
</html>
"""


# ---------- by_cluster ------------------------------------------------------


def _read_cluster_by_sample(run_dir: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    """
    Parse ``02_integrate/cluster_by_sample.csv`` if present.
    Returns (sample_columns, mapping[cluster -> {sample: count}]).
    """
    p = run_dir / STAGE_DIRS["integrate"] / "cluster_by_sample.csv"
    if not p.exists():
        return [], {}
    samples: list[str] = []
    out: dict[str, dict[str, str]] = {}
    with p.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n")
        parts = _split_csv_line(header)
        # First column is the cluster id (no name in pandas crosstab).
        samples = parts[1:]
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cells = _split_csv_line(line)
            if len(cells) < 2:
                continue
            cluster = cells[0]
            out[cluster] = {samples[i]: cells[i + 1] for i in range(len(samples)) if i + 1 < len(cells)}
    return samples, out


def _by_cluster_html(
    cluster_id: str,
    annotation: dict[str, str] | None,
    markers_for_cluster: list[dict[str, str]],
    sample_dist: tuple[list[str], dict[str, str]] | None,
    used_ai: bool,
) -> str:
    label = annotation.get("panel_label", "") if annotation else ""

    # Top 10 markers for this cluster
    if markers_for_cluster:
        rows = []
        for r in markers_for_cluster[:10]:
            rows.append(
                "<tr>"
                f"<td>{_e(r.get('gene', ''))}</td>"
                f"<td>{_e(r.get('log2fc') or r.get('logfoldchanges') or '')}</td>"
                f"<td>{_e(r.get('pct_in') or r.get('pct.1') or r.get('pct.in') or '')}</td>"
                f"<td>{_e(r.get('pct_out') or r.get('pct.2') or r.get('pct.out') or '')}</td>"
                "</tr>"
            )
        markers_table = (
            "<table class=\"plain\">"
            "<thead><tr><th>gene</th><th>log2fc</th><th>pct in</th><th>pct out</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody>"
            "</table>"
        )
    else:
        markers_table = "<p><em>No markers recorded.</em></p>"

    if annotation:
        why_rows = [
            ("panel_score", annotation.get("panel_score", "")),
            ("panel_margin", annotation.get("panel_margin", "")),
        ]
        if annotation.get("panel_rationale"):
            why_rows.append(("panel_rationale", annotation.get("panel_rationale", "")))
        if used_ai or annotation.get("ai_label"):
            why_rows.append(("ai_label", annotation.get("ai_label", "")))
            why_rows.append(("ai_rationale", annotation.get("ai_rationale", "")))
        why_table = (
            "<table class=\"plain\">"
            "<tbody>"
            + "".join(
                f"<tr><th>{_e(k)}</th><td>{_e(v)}</td></tr>"
                for k, v in why_rows
            )
            + "</tbody></table>"
        )
    else:
        why_table = "<p><em>No annotation row for this cluster.</em></p>"

    if sample_dist is not None and sample_dist[1]:
        samples, dist = sample_dist
        cells_total = 0
        for v in dist.values():
            try:
                cells_total += int(v)
            except (TypeError, ValueError):
                continue
        rows_dist = "".join(
            f"<tr><td>{_e(s)}</td><td>{_e(dist.get(s, ''))}</td></tr>"
            for s in samples
        )
        dist_section = (
            f"<p>Total cells in cluster: <strong>{cells_total}</strong></p>"
            "<table class=\"plain\">"
            "<thead><tr><th>sample</th><th>n_cells</th></tr></thead>"
            f"<tbody>{rows_dist}</tbody>"
            "</table>"
        )
    else:
        dist_section = "<p><em>No cluster_by_sample.csv on disk (single-sample run, or integrate skipped the contingency table).</em></p>"

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>scellrun — cluster {_e(cluster_id)}{(" — " + _e(label)) if label else ""}</title>
{_VIEWS_CSS}
</head>
<body>
<p><a href="../index.html">&larr; views index</a></p>
<h1>Cluster {_e(cluster_id)}{(" &mdash; " + _e(label)) if label else ""}</h1>

<h2>Top 10 markers</h2>
{markers_table}

<h2>Why this label</h2>
{why_table}

<h2>Cell count + sample distribution</h2>
{dist_section}

<p class="footer">
  <a href="cluster_{_e(cluster_id)}_cells.csv">Cell barcodes (CSV)</a>
</p>
</body>
</html>
"""


# ---------- by_decision_source ----------------------------------------------


def _is_boring(d: dict) -> bool:
    """
    Return True if this auto row is a "boring" knob default — same
    threshold every run, no value to a reviewer. Skipped on the
    by_decision_source page to keep the auto list useful.
    """
    if d.get("source") != "auto":
        return False
    key = str(d.get("key", ""))
    if key in _BORING_AUTO_KEYS:
        return True
    rationale = (d.get("rationale") or "").strip()
    return not rationale


def _decisions_table_html(rows: list[dict]) -> str:
    if not rows:
        return "<p><em>No rows.</em></p>"
    body = []
    for d in rows:
        body.append(
            "<tr>"
            f"<td><code>{_e(d.get('stage', ''))}.{_e(d.get('key', ''))}</code></td>"
            f"<td>{_e(d.get('value', ''))}</td>"
            f"<td>{_e(d.get('rationale', ''))}</td>"
            "</tr>"
        )
    return (
        "<table class=\"plain\">"
        "<thead><tr><th>key</th><th>value</th><th>rationale</th></tr></thead>"
        f"<tbody>{''.join(body)}</tbody>"
        "</table>"
    )


def _by_decision_source_html(decisions: list[dict]) -> str:
    ai_rows = [d for d in decisions if d.get("source") == "ai"]
    user_rows = [d for d in decisions if d.get("source") == "user"]
    auto_rows = [d for d in decisions if d.get("source") == "auto" and not _is_boring(d)]

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>scellrun — decisions by source</title>
{_VIEWS_CSS}
</head>
<body>
<p><a href="../index.html">&larr; views index</a></p>
<h1>Decisions by source</h1>

<h2>AI-driven ({len(ai_rows)})</h2>
{_decisions_table_html(ai_rows)}

<h2>User overrides ({len(user_rows)})</h2>
{_decisions_table_html(user_rows)}

<h2>Auto picks with rationale ({len(auto_rows)})</h2>
<p>Boring knob defaults (min_genes, max_genes, etc.) are filtered out.</p>
{_decisions_table_html(auto_rows)}
</body>
</html>
"""


# ---------- index page ------------------------------------------------------


def _index_html(
    resolutions: list[float],
    annotations: list[dict[str, str]],
    annotate_res: float | None,
) -> str:
    res_links = "".join(
        f'<li><a href="by_resolution/{_e(_format_resolution_segment(r))}/index.html">'
        f"resolution {_e(_format_resolution_segment(r))}</a></li>"
        for r in resolutions
    )

    annot_table = ""
    if annotations:
        rows = []
        for a in annotations:
            cluster = a.get("cluster", "")
            label = a.get("panel_label", "")
            rows.append(
                "<tr>"
                f'<td><a href="by_cluster/cluster_{_e(cluster)}.html">cluster {_e(cluster)}</a></td>'
                f"<td>{_e(label)}</td>"
                "</tr>"
            )
        annot_table = (
            "<table class=\"plain\">"
            "<thead><tr><th>cluster</th><th>panel_label</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody>"
            "</table>"
        )
    else:
        annot_table = "<p><em>No annotations recorded.</em></p>"

    annotate_res_label = (
        f" (at the chosen annotate resolution {_format_resolution_segment(annotate_res)})"
        if annotate_res is not None
        else ""
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>scellrun — views</title>
{_VIEWS_CSS}
</head>
<body>
<p><a href="../{STAGE_DIRS['report']}/index.html">&larr; main report</a></p>
<h1>Views</h1>
<p>Cross-sections of this run, derived from the same artifacts the main report references.</p>

<h2>By resolution</h2>
{('<ul>' + res_links + '</ul>') if res_links else '<p><em>No resolutions on disk.</em></p>'}

<h2>By cluster{_e(annotate_res_label)}</h2>
{annot_table}

<h2>By decision source</h2>
<p><a href="by_decision_source/index.html">decisions grouped by source (ai / user / auto)</a></p>
</body>
</html>
"""


# ---------- shared CSS ------------------------------------------------------

_VIEWS_CSS = """<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1080px; margin: 2rem auto; padding: 0 1rem;
         color: #222; line-height: 1.55; }
  h1 { font-size: 1.6rem; border-bottom: 2px solid #333; padding-bottom: .3rem; }
  h2 { font-size: 1.1rem; margin-top: 1.6rem; padding-bottom: .15rem;
       border-bottom: 1px solid #eee; }
  a { color: #1f4e8a; }
  table.plain { font-size: .9rem; border-collapse: collapse; margin: .4rem 0; }
  table.plain th, table.plain td { border: 1px solid #ddd; padding: .25rem .5rem;
                                   text-align: left; vertical-align: top; }
  table.plain th { background: #f6f6f6; }
  img { max-width: 100%; height: auto; }
  code { background: #f0f0f0; padding: .05rem .3rem; border-radius: 3px; font-size: .9em; }
  .footer { font-size: .85rem; color: #666; margin-top: 1.5rem; }
</style>"""


# ---------- public entrypoint -----------------------------------------------


def build_views(run_dir: Path) -> dict[str, Path]:
    """
    Generate ``06_views/`` HTML index files. References stage artifacts via
    relative paths (``../02_integrate/...``); no file copying, no symlinks.

    Returns a mapping of view-name → absolute path of the file written.
    Best-effort: if a stage's artifacts are missing the corresponding view
    section degrades to a placeholder; the function never raises on
    missing inputs (it's an aggregator, not a stage).
    """
    run_dir = Path(run_dir)
    views_root = run_dir / STAGE_DIRS["views"]
    views_root.mkdir(parents=True, exist_ok=True)

    out: dict[str, Path] = {}

    decisions = _read_decisions_jsonl(run_dir)
    resolutions = _resolutions_present(run_dir)
    annotations = _read_annotations(run_dir)
    annotate_res = _annotate_resolution(decisions)
    samples, cluster_by_sample = _read_cluster_by_sample(run_dir)

    # by_resolution/<R>/index.html
    by_res_root = views_root / "by_resolution"
    by_res_root.mkdir(parents=True, exist_ok=True)
    for r in resolutions:
        res_dir = by_res_root / _format_resolution_segment(r)
        res_dir.mkdir(parents=True, exist_ok=True)
        markers_p = run_dir / STAGE_DIRS["markers"] / f"markers_res_{r:g}.csv"
        markers_rows = _read_markers_csv(markers_p)
        # Annotations only attach to the resolution that was annotated.
        annots_for_res = annotations if (annotate_res is not None and r == annotate_res) else None
        page = res_dir / "index.html"
        page.write_text(
            _by_resolution_html(r, markers_rows, annots_for_res),
            encoding="utf-8",
        )
        out[f"by_resolution_{r:g}"] = page

    # by_cluster/cluster_<C>.html
    by_clu_root = views_root / "by_cluster"
    by_clu_root.mkdir(parents=True, exist_ok=True)
    used_ai = any(
        d.get("stage") == "annotate" and d.get("key") == "use_ai" and d.get("value") is True
        for d in decisions
    )
    # Markers for the annotate resolution (if known); fall back to the
    # markers CSV with the most rows so the by-cluster page still
    # populates when analyze didn't record chosen_resolution_for_annotate.
    annot_markers: list[dict[str, str]] = []
    if annotate_res is not None:
        annot_markers = _read_markers_csv(
            run_dir / STAGE_DIRS["markers"] / f"markers_res_{annotate_res:g}.csv"
        )
    elif resolutions:
        # Best-effort fallback: pick the first resolution on disk.
        annot_markers = _read_markers_csv(
            run_dir / STAGE_DIRS["markers"] / f"markers_res_{resolutions[0]:g}.csv"
        )

    annot_by_cluster: dict[str, dict[str, str]] = {
        a.get("cluster", ""): a for a in annotations
    }
    markers_by_cluster: dict[str, list[dict[str, str]]] = {}
    for r in annot_markers:
        c = str(r.get("cluster", ""))
        markers_by_cluster.setdefault(c, []).append(r)

    cluster_ids = list(annot_by_cluster.keys())
    if not cluster_ids:
        # No annotations.csv on disk; fall back to the markers' clusters.
        cluster_ids = list(markers_by_cluster.keys())

    for cluster_id in cluster_ids:
        if not cluster_id:
            continue
        page = by_clu_root / f"cluster_{cluster_id}.html"
        sample_dist: tuple[list[str], dict[str, str]] | None = None
        if cluster_id in cluster_by_sample:
            sample_dist = (samples, cluster_by_sample[cluster_id])
        page.write_text(
            _by_cluster_html(
                cluster_id=cluster_id,
                annotation=annot_by_cluster.get(cluster_id),
                markers_for_cluster=markers_by_cluster.get(cluster_id, []),
                sample_dist=sample_dist,
                used_ai=used_ai,
            ),
            encoding="utf-8",
        )
        out[f"by_cluster_{cluster_id}"] = page

        # Per-cluster cell barcode CSV — empty by default; analyze-time
        # callers can populate the cells from adata.obs if they want a
        # full barcode dump. Writing the file unconditionally so the
        # footer link in the by_cluster page never 404s on disk.
        cells_path = by_clu_root / f"cluster_{cluster_id}_cells.csv"
        if not cells_path.exists():
            cells_path.write_text(
                "barcode\n",
                encoding="utf-8",
            )

    # by_decision_source/index.html
    by_src_root = views_root / "by_decision_source"
    by_src_root.mkdir(parents=True, exist_ok=True)
    src_page = by_src_root / "index.html"
    src_page.write_text(_by_decision_source_html(decisions), encoding="utf-8")
    out["by_decision_source"] = src_page

    # 06_views/index.html landing page
    landing = views_root / "index.html"
    landing.write_text(
        _index_html(resolutions, annotations, annotate_res),
        encoding="utf-8",
    )
    out["index"] = landing

    return out
