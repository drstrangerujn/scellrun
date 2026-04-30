"""
scrna annotate stage: assign cell-type labels to clusters.

Two-tier strategy:
- Deterministic panel match (always runs): for every cluster, compute a
  per-panel-group score from the cluster's top markers vs. the profile's
  marker dict. Best-scoring panel = label.
- Optional AI enhance (`--ai`): send (top markers + tissue + profile context)
  to Anthropic API; the model returns a celltype + reasoning. The result
  is reported alongside the deterministic label so the user can compare,
  not delegate.

PubMed evidence: for each cluster's top markers, fetch top recent papers
referencing both the gene and the tissue keyword. Embedded into the report
as the third "literature" tier of the data / inference / literature
provenance trail (`feedback_report_rules`).

Maps to Rmd § 9-10. Profiles ship the panels; this module is the matching
+ explanation engine, not a panel definition.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import urlopen

import anndata as ad
import pandas as pd

from scellrun.decisions import Decision, record_many
from scellrun.self_check import SelfCheckFinding, annotate_self_check, record_findings

DEFAULT_TOP_N_MARKERS = 30  # how many top markers per cluster to consider for matching
DEFAULT_PUBMED_PER_GENE = 3
PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class ClusterAnnotation:
    cluster: str
    panel_label: str
    panel_score: float
    panel_margin: float  # gap to runner-up; a small margin = ambiguous
    panel_rationale: str
    ai_label: str | None = None
    ai_rationale: str | None = None
    top_markers: list[str] = field(default_factory=list)
    pubmed: dict[str, list[dict]] = field(default_factory=dict)


@dataclass
class AnnotateResult:
    resolution: float
    panel_name: str
    annotations: list[ClusterAnnotation]
    used_ai: bool
    used_pubmed: bool
    tissue: str | None
    findings: list[SelfCheckFinding] = field(default_factory=list)
    """v0.8 self-check findings (panel ambiguity / tissue mismatch)."""


def _resolution_keys(adata: ad.AnnData) -> dict[float, str]:
    out: dict[float, str] = {}
    for col in adata.obs.columns:
        if not col.startswith("leiden_res_"):
            continue
        try:
            out[float(col.removeprefix("leiden_res_").replace("_", "."))] = col
        except ValueError:
            continue
    return out


def _top_markers_per_cluster(
    adata: ad.AnnData,
    cluster_key: str,
    *,
    top_n: int = DEFAULT_TOP_N_MARKERS,
    use_raw: bool | None = None,
) -> dict[str, list[str]]:
    """Run rank_genes_groups (or reuse if already cached) and pull top genes per cluster."""
    import scanpy as sc

    if use_raw is None:
        use_raw = adata.raw is not None

    rank_key = f"rank_{cluster_key}"
    if rank_key not in adata.uns:
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            pts=True,
            key_added=rank_key,
        )

    rg = adata.uns[rank_key]
    cluster_names = list(rg["names"].dtype.names)

    out: dict[str, list[str]] = {}
    for c in cluster_names:
        names = list(rg["names"][c][:top_n])
        out[c] = [str(g) for g in names]
    return out


def _score_cluster_against_panel(
    cluster_top_markers: list[str],
    panel: dict[str, list[str]],
) -> dict[str, float]:
    """
    Score: for each panel group, how many of its marker genes appear in the
    cluster's top-N markers, divided by the panel size. 0..1 range.

    A more sophisticated version would use scanpy.tl.score_genes, but for
    label assignment this overlap-fraction is robust and dropout-resistant.
    """
    scores: dict[str, float] = {}
    cluster_set = set(cluster_top_markers)
    for label, genes in panel.items():
        if not genes:
            scores[label] = 0.0
            continue
        hits = sum(1 for g in genes if g in cluster_set)
        scores[label] = hits / len(genes)
    return scores


def _best_panel_match(scores: dict[str, float]) -> tuple[str, float, float]:
    """Return (best_label, best_score, margin_to_runner_up).

    When every panel group scores 0 (none of the panel's genes appear in the
    cluster's top markers), return "Unassigned" rather than the first panel
    label by dict-iteration order. The previous behavior silently labeled
    every score-0 cluster with the first key in the panel, so for the
    joint-disease chondrocyte_markers panel a non-chondrocyte cluster (e.g.
    NK / T / mast) would come back as "ProC" with score 0.0 — see
    ISSUES.md #003 from the v0.7 OA dogfood.
    """
    if not scores:
        return ("Unassigned", 0.0, 0.0)
    sorted_labels = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    best_label, best_score = sorted_labels[0]
    runner_up = sorted_labels[1][1] if len(sorted_labels) > 1 else 0.0
    if best_score == 0.0:
        return ("Unassigned", 0.0, 0.0)
    return (best_label, best_score, best_score - runner_up)


def _select_panel(profile_module: Any, panel_name: str | None) -> tuple[str, dict[str, list[str]]]:
    """
    Pick the panel to match against. If panel_name is None, prefer
    chondrocyte_markers (fine subtype) when available, else celltype_broad.
    """
    if panel_name is not None:
        if not hasattr(profile_module, panel_name):
            raise ValueError(
                f"profile {profile_module.__name__!r} doesn't define panel {panel_name!r}"
            )
        return panel_name, getattr(profile_module, panel_name)

    for candidate in ("chondrocyte_markers", "celltype_broad"):
        if hasattr(profile_module, candidate):
            return candidate, getattr(profile_module, candidate)
    raise ValueError(
        f"profile {profile_module.__name__!r} has no marker panels — "
        "annotate needs at least one of: chondrocyte_markers, celltype_broad."
    )


def _fetch_pubmed_for_marker(
    gene: str,
    tissue: str | None,
    *,
    n: int = DEFAULT_PUBMED_PER_GENE,
    timeout: int = 5,
) -> list[dict]:
    """
    Best-effort PubMed lookup. Failures (network, rate limit) are non-fatal —
    we return [] and the report omits the evidence section for that gene.
    """
    if tissue:
        term = f"{gene}[Title/Abstract] AND {tissue}[Title/Abstract]"
    else:
        term = f"{gene}[Title/Abstract]"
    try:
        # esearch
        params = urlencode({
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "retmax": n,
            "sort": "date",
        })
        with urlopen(f"{PUBMED_EUTILS}/esearch.fcgi?{params}", timeout=timeout) as r:
            import json
            data = json.loads(r.read().decode())
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        # esummary
        params = urlencode({
            "db": "pubmed",
            "id": ",".join(ids),
            "retmode": "json",
        })
        with urlopen(f"{PUBMED_EUTILS}/esummary.fcgi?{params}", timeout=timeout) as r:
            summ = json.loads(r.read().decode())

        out = []
        for pmid in ids:
            doc = summ.get("result", {}).get(pmid, {})
            if not doc:
                continue
            out.append({
                "pmid": pmid,
                "title": doc.get("title", ""),
                "year": (doc.get("pubdate", "") or "").split(" ")[0],
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })
        return out
    except Exception:
        return []


def _ai_annotate(
    cluster_id: str,
    top_markers: list[str],
    panel_label: str,
    panel_score: float,
    profile_name: str,
    tissue: str | None,
    *,
    model: str = "claude-haiku-4-5-20251001",
) -> tuple[str | None, str | None]:
    """
    Optional LLM second-opinion. Returns (label, rationale). Failures are
    non-fatal — caller sees None, None and the report omits the AI column.

    Requires ANTHROPIC_API_KEY env var. If not present, returns (None, None)
    with a note in the rationale slot.
    """
    import os

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return (None, "ANTHROPIC_API_KEY not set; AI step skipped")

    try:
        from anthropic import Anthropic
    except ImportError:
        return (None, "anthropic package not installed; AI step skipped (pip install anthropic)")

    tissue_str = f" The tissue context is **{tissue}**." if tissue else ""
    prompt = f"""You are reviewing a single-cell RNA-seq cluster annotation.

Profile in use: {profile_name}.{tissue_str}
Cluster ID: {cluster_id}
Top markers (in order): {", ".join(top_markers[:20])}

The deterministic panel-match step assigned this cluster the label "{panel_label}" with a score of {panel_score:.2f} (overlap fraction with that panel's gene list).

Question: based on the marker list, do you agree with "{panel_label}"? If not, what's a better cell-type label given the profile context?

Respond in this exact form:
LABEL: <one or two words, e.g. "Macrophages" or "preInfC chondrocytes">
RATIONALE: <one or two sentences explaining the markers that drove the call>

Be conservative: if the markers don't clearly support any cell type, label it "Ambiguous" and explain why.
"""

    try:
        client = Anthropic(api_key=api_key)
        msg = client.messages.create(
            model=model,
            max_tokens=300,
            messages=[{"role": "user", "content": prompt}],
        )
        content = msg.content[0].text if msg.content else ""
    except Exception as e:
        return (None, f"AI call failed: {type(e).__name__}: {e}")

    label = None
    rationale = None
    for line in content.splitlines():
        line = line.strip()
        if line.startswith("LABEL:"):
            label = line.removeprefix("LABEL:").strip()
        elif line.startswith("RATIONALE:"):
            rationale = line.removeprefix("RATIONALE:").strip()
    return (label, rationale)


def run_annotate(
    adata: ad.AnnData,
    profile_module: Any,
    *,
    resolution: float,
    panel_name: str | None = None,
    use_ai: bool = False,
    ai_model: str = "claude-haiku-4-5-20251001",
    use_pubmed: bool = False,
    tissue: str | None = None,
    top_n_markers: int = DEFAULT_TOP_N_MARKERS,
    pubmed_per_gene: int = DEFAULT_PUBMED_PER_GENE,
    pubmed_top_genes: int = 3,
    run_dir: Path | None = None,
) -> AnnotateResult:
    """
    Annotate every cluster at the given resolution.

    `profile_module` is the loaded scellrun.profiles.* submodule. Must define
    at least one marker panel (chondrocyte_markers or celltype_broad).
    """
    res_to_key = _resolution_keys(adata)
    if resolution not in res_to_key:
        raise ValueError(
            f"resolution {resolution} not present in adata.obs. "
            f"Available: {sorted(res_to_key)}"
        )

    panel_actual_name, panel = _select_panel(profile_module, panel_name)

    cluster_key = res_to_key[resolution]
    top_per_cluster = _top_markers_per_cluster(adata, cluster_key, top_n=top_n_markers)

    annotations: list[ClusterAnnotation] = []
    for cluster_id, top_markers in sorted(top_per_cluster.items(), key=lambda kv: kv[0]):
        scores = _score_cluster_against_panel(top_markers, panel)
        best_label, best_score, margin = _best_panel_match(scores)

        rationale_bits = []
        for label, score in sorted(scores.items(), key=lambda kv: kv[1], reverse=True)[:3]:
            matched_genes = [g for g in panel.get(label, []) if g in set(top_markers)]
            if matched_genes:
                rationale_bits.append(f"{label} ({score:.2f}): {', '.join(matched_genes)}")
            else:
                rationale_bits.append(f"{label} ({score:.2f}): —")
        panel_rationale = " | ".join(rationale_bits)

        ai_label, ai_rationale = (None, None)
        if use_ai:
            ai_label, ai_rationale = _ai_annotate(
                cluster_id, top_markers, best_label, best_score,
                profile_module.__name__.split(".")[-1], tissue, model=ai_model,
            )

        pubmed: dict[str, list[dict]] = {}
        if use_pubmed:
            for gene in top_markers[:pubmed_top_genes]:
                pubmed[gene] = _fetch_pubmed_for_marker(gene, tissue, n=pubmed_per_gene)
                # NCBI rate limit: 3/sec without API key
                time.sleep(0.4)

        annotations.append(ClusterAnnotation(
            cluster=cluster_id,
            panel_label=best_label,
            panel_score=best_score,
            panel_margin=margin,
            panel_rationale=panel_rationale,
            ai_label=ai_label,
            ai_rationale=ai_rationale,
            top_markers=top_markers,
            pubmed=pubmed,
        ))

    # Write back to adata.obs as a new category column
    cluster_col = res_to_key[resolution]
    label_map_panel = {a.cluster: a.panel_label for a in annotations}
    adata.obs["scellrun_celltype_panel"] = (
        adata.obs[cluster_col].map(label_map_panel).astype("category")
    )
    if use_ai:
        label_map_ai = {a.cluster: (a.ai_label or a.panel_label) for a in annotations}
        adata.obs["scellrun_celltype_ai"] = (
            adata.obs[cluster_col].map(label_map_ai).astype("category")
        )

    # v0.8 self-check: panel-margin and panel-tissue-mismatch guards.
    findings = annotate_self_check(
        annotations=annotations,
        panel_name=panel_actual_name,
        profile_module=profile_module,
    )

    if run_dir is not None:
        profile_short = profile_module.__name__.split(".")[-1]
        decisions: list[Decision] = [
            Decision(
                stage="annotate",
                key="profile",
                value=profile_short,
                default="default",
                source="user" if profile_short != "default" else "auto",
                rationale="profile selects which marker panels are available for matching",
            ),
            Decision(
                stage="annotate",
                key="panel",
                value=panel_actual_name,
                default=None,
                source="user" if panel_name is not None else "auto",
                rationale=(
                    f"--panel {panel_name!r} forced"
                    if panel_name is not None
                    else (
                        "auto-picked chondrocyte_markers (fine subtype) — preferred when available"
                        if panel_actual_name == "chondrocyte_markers"
                        else f"auto-picked {panel_actual_name!r} — first panel the profile defines"
                    )
                ),
            ),
            Decision(
                stage="annotate",
                key="resolution",
                value=resolution,
                default=None,
                source="user",
                rationale=(
                    f"matching against leiden_res_{resolution:g} clusters; the orchestrator picks "
                    "this from integrate's quality table when called via `scellrun analyze`"
                ),
            ),
            Decision(
                stage="annotate",
                key="use_ai",
                value=use_ai,
                default=False,
                source="user" if use_ai else "auto",
                rationale=(
                    f"AI second-opinion enabled (model={ai_model}); deterministic panel call still authoritative"
                    if use_ai
                    else "AI second-opinion off — deterministic panel match only"
                ),
            ),
            Decision(
                stage="annotate",
                key="use_pubmed",
                value=use_pubmed,
                default=False,
                source="user" if use_pubmed else "auto",
                rationale=(
                    f"PubMed evidence column on; tissue={tissue!r}; up to "
                    f"{pubmed_per_gene} papers per top gene"
                    if use_pubmed
                    else "PubMed lookup off — turn on with --pubmed for the literature evidence column"
                ),
            ),
            Decision(
                stage="annotate",
                key="tissue",
                value=tissue,
                default=None,
                source="user" if tissue else "auto",
                rationale=(
                    f"tissue context {tissue!r} drives PubMed scoping and AI prompt"
                    if tissue
                    else "no tissue context supplied; PubMed scoping and AI prompt run un-anchored"
                ),
            ),
        ]
        record_many(run_dir, decisions)
        record_findings(run_dir, findings)

    return AnnotateResult(
        resolution=resolution,
        panel_name=panel_actual_name,
        annotations=annotations,
        used_ai=use_ai,
        used_pubmed=use_pubmed,
        tissue=tissue,
        findings=findings,
    )


def write_artifacts(
    result: AnnotateResult,
    adata: ad.AnnData,
    out_dir: Path,
    *,
    write_h5ad: bool = True,
    lang: str = "en",
) -> dict[str, Path]:
    """Write annotations.csv + annotated h5ad + HTML report."""
    from jinja2 import Environment, PackageLoader, select_autoescape

    out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, Path] = {}

    rows = []
    for a in result.annotations:
        rows.append({
            "cluster": a.cluster,
            "panel_label": a.panel_label,
            "panel_score": a.panel_score,
            "panel_margin": a.panel_margin,
            "ai_label": a.ai_label or "",
            "ai_rationale": a.ai_rationale or "",
            "top_markers": ", ".join(a.top_markers[:10]),
        })
    csv_path = out_dir / "annotations.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    artifacts["annotations"] = csv_path

    if write_h5ad:
        h5_path = out_dir / "annotated.h5ad"
        adata.write_h5ad(h5_path)
        artifacts["annotated_h5ad"] = h5_path

    env = Environment(
        loader=PackageLoader("scellrun", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    template_name = "scrna_annotate_zh.html.j2" if lang == "zh" else "scrna_annotate.html.j2"
    template = env.get_template(template_name)
    html = template.render(result=result)
    html_path = out_dir / "report.html"
    html_path.write_text(html, encoding="utf-8")
    artifacts["report"] = html_path
    return artifacts
