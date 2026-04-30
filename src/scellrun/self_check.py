"""
Stage self-check primitives (v0.8).

Each stage runs a `self_check()` at the end of its `run_*()` to spot common
failure modes (most cells failing QC, no clusters, panel doesn't match data)
and surface a concrete recommendation rather than a stack trace.

A self-check produces zero-or-more `SelfCheckFinding`s. Each finding is
appended to the decision log as two `source="auto"` rows: one for the
trigger (what looked wrong) and one for the suggestion (what to try). The
orchestrator (`scellrun analyze`) can additionally consume the structured
finding to apply a fix and re-run the stage when `--auto-fix` is on; cap is
1 retry per stage to avoid loops.

Findings are intentionally a small dataclass (not free text) so the
auto-fix path can act on them without re-parsing rationale prose.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from scellrun.decisions import Decision, record_many


@dataclass
class SelfCheckFinding:
    """One thing the stage flagged about its own output."""

    stage: str
    """Pipeline stage that ran the check (qc / integrate / annotate)."""

    code: str
    """Short slug naming the issue (e.g. 'qc_low_pass_rate', 'integrate_too_few_clusters')."""

    trigger: str
    """One-sentence description of what looked wrong (rendered as the trigger decision rationale)."""

    suggestion: str
    """One-sentence description of what to try (rendered as the suggestion decision rationale)."""

    fix: dict[str, Any] = field(default_factory=dict)
    """
    Structured params the orchestrator can apply when --auto-fix is on.
    Empty dict means the suggestion is human-only (no automatic recovery
    available). Concrete keys per stage:
      qc:        {"max_pct_mt": 25} or {"max_genes": 6000}
      integrate: {"resolutions": "aio"} or {"regress_cell_cycle": True}
      annotate:  {"panel_name": "celltype_broad"}
    """


def record_findings(
    run_dir: Path | None,
    findings: list[SelfCheckFinding],
    *,
    attempt_id: str = "",
) -> None:
    """
    Append every finding to the decision log as two rows: trigger + suggestion.

    Both rows carry `source="auto"` so the report's decision summary
    surfaces them under the stage that ran the check, ordered after the
    stage's regular decisions (since the check runs at the end).

    The suggestion row carries the structured `fix` dict in the
    Decision's `fix_payload` field so downstream tools can mechanically
    apply the fix without having to parse the rationale prose.
    """
    if not findings:
        return
    decisions: list[Decision] = []
    for f in findings:
        decisions.append(Decision(
            stage=f.stage,
            key=f"self_check.{f.code}.trigger",
            value=f.code,
            default=None,
            source="auto",
            attempt_id=attempt_id,
            rationale=f.trigger,
        ))
        decisions.append(Decision(
            stage=f.stage,
            key=f"self_check.{f.code}.suggest",
            value=_summarize_fix(f),
            default=None,
            source="auto",
            attempt_id=attempt_id,
            rationale=f.suggestion,
            fix_payload=dict(f.fix) if f.fix else None,
        ))
    record_many(run_dir, decisions)


def _summarize_fix(f: SelfCheckFinding) -> str:
    """
    Render a fix dict as the human-readable string stored in Decision.value.

    Examples:
        {"max_pct_mt": 25}            -> "raise --max-pct-mt to 25"
        {"max_genes": 6000}           -> "raise --max-genes to 6000"
        {"resolutions": "aio"}        -> "switch to --resolutions aio (13-step sweep)"
        {"regress_cell_cycle": True}  -> "turn on --regress-cell-cycle"
        {"panel_name": "celltype_broad"} -> "switch --panel to celltype_broad"
        {} (no auto-fix path)         -> the suggestion text itself, verbatim
    """
    if not f.fix:
        return f.suggestion
    if "max_pct_mt" in f.fix:
        return f"raise --max-pct-mt to {f.fix['max_pct_mt']}"
    if "max_genes" in f.fix:
        return f"raise --max-genes to {f.fix['max_genes']}"
    if "min_genes" in f.fix:
        return f"lower --min-genes to {f.fix['min_genes']}"
    if "min_counts" in f.fix:
        return f"lower --min-counts to {f.fix['min_counts']}"
    if "resolutions" in f.fix:
        return f"switch to --resolutions {f.fix['resolutions']} (wider sweep)"
    if "regress_cell_cycle" in f.fix and f.fix["regress_cell_cycle"]:
        return "turn on --regress-cell-cycle"
    if "panel_name" in f.fix:
        return f"switch --panel to {f.fix['panel_name']}"
    if "profile" in f.fix:
        return f"switch --profile to {f.fix['profile']}"
    # Fallback: dump the dict
    return ", ".join(f"{k}={v}" for k, v in f.fix.items())


# ---------------------------------------------------------------------------
# QC self-check
# ---------------------------------------------------------------------------

QC_PASS_RATE_TRIGGER = 60.0  # below this %, fire the check (v0.9.1: was 30%)
QC_PASS_RATE_TARGET = 60.0   # find the cheapest single threshold change to reach this %


def qc_self_check(
    n_cells_in: int,
    n_cells_pass: int,
    sensitivity: dict[str, list[dict]],
    thresholds: Any,
) -> list[SelfCheckFinding]:
    """
    Walk QCResult.sensitivity for the smallest single threshold change that
    pushes pass-rate to QC_PASS_RATE_TARGET (60%) when the current rate is
    below QC_PASS_RATE_TRIGGER (60%, lowered from 30% in v0.9.1).

    The 60% trigger ceiling catches the degraded-prep range (30-60% pass)
    that the v0.8 30% trigger silently waved through. Real OA samples
    routinely sit at 40-55% on default thresholds because joint tissue
    is stress-prone; that's prime suggestion territory.

    Strategy: for each candidate knob (max_pct_mt, max_genes), find the
    smallest threshold in the sensitivity table where the marginal pass-rate
    on JUST that axis reaches QC_PASS_RATE_TARGET. Pick the suggestion
    whose threshold is the smallest absolute relaxation from the current
    setting (cheapest single change).

    Note: sensitivity rows count cells passing JUST that one threshold, not
    the AND of every threshold. So a row at 60% pass on max_pct_mt may not
    yield 60% overall pass after the other knobs still apply. This is OK
    for v0.8 — the suggestion is the cheapest *single* change; the user
    re-runs and sees the new pass-rate. The auto-fix loop then fires once
    and accepts whatever pass-rate the relaxation produces.
    """
    findings: list[SelfCheckFinding] = []
    if n_cells_in == 0:
        return findings

    pass_rate = 100.0 * n_cells_pass / n_cells_in
    if pass_rate >= QC_PASS_RATE_TRIGGER:
        return findings

    # Candidate relaxations: (knob name, current value, sensitivity rows, comparison direction)
    # comparison: "ge" means we want threshold > current (relax UPWARDS).
    #             "le" means we want threshold < current (relax DOWNWARDS).
    candidates: list[tuple[str, float, list[dict], str]] = []
    if "max_pct_mt" in sensitivity:
        candidates.append((
            "max_pct_mt",
            float(thresholds.max_pct_mt),
            sensitivity["max_pct_mt"],
            "ge",
        ))
    if "max_genes" in sensitivity:
        candidates.append((
            "max_genes",
            float(thresholds.max_genes),
            sensitivity["max_genes"],
            "ge",
        ))

    best: tuple[str, float, float] | None = None  # (knob, suggested_value, marginal_pct)

    for knob, current, rows, direction in candidates:
        for row in rows:
            t = float(row["threshold"])
            pct = float(row["pct_pass"])
            # Only consider relaxations (raise the ceiling), not tightenings.
            if direction == "ge" and t <= current:
                continue
            if direction == "le" and t >= current:
                continue
            if pct >= QC_PASS_RATE_TARGET:
                # Found the smallest relaxation for this knob (rows are sorted asc).
                if best is None or (t - current) < (best[1] - current):
                    best = (knob, t, pct)
                break  # first row hitting target on this knob is the smallest

    if best is None:
        # No single relaxation gets us to 60%; surface a human-only suggestion.
        findings.append(SelfCheckFinding(
            stage="qc",
            code="qc_low_pass_rate_no_easy_fix",
            trigger=(
                f"only {pass_rate:.1f}% of cells passed QC ({n_cells_pass}/{n_cells_in}); "
                "the sensitivity sweep shows no single threshold relaxation reaches "
                f"{QC_PASS_RATE_TARGET:.0f}%"
            ),
            suggestion=(
                "review the per-cell-metrics CSV — multiple QC axes (mt%, n_genes, n_counts) "
                "are likely failing together; consider re-dissociating the sample or accepting "
                "the lower pass-rate explicitly via --max-pct-mt / --max-genes overrides"
            ),
            fix={},
        ))
        return findings

    knob, sugg_value, marginal_pct = best
    fix_value: int | float = int(sugg_value) if knob == "max_genes" else sugg_value
    findings.append(SelfCheckFinding(
        stage="qc",
        code="qc_low_pass_rate",
        trigger=(
            f"only {pass_rate:.1f}% of cells passed QC ({n_cells_pass}/{n_cells_in}); "
            f"below the {QC_PASS_RATE_TRIGGER:.0f}% trigger threshold"
        ),
        suggestion=(
            f"raising --{knob.replace('_', '-')} from {thresholds_current_str(knob, thresholds)} "
            f"to {fix_value} would let {marginal_pct:.0f}% of cells pass that single test "
            f"(sensitivity sweep, smallest relaxation reaching {QC_PASS_RATE_TARGET:.0f}%)"
        ),
        fix={knob: fix_value},
    ))
    return findings


def thresholds_current_str(knob: str, thresholds: Any) -> str:
    """Format the current value of a threshold knob for the suggestion text."""
    v = getattr(thresholds, knob)
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


# ---------------------------------------------------------------------------
# Integrate self-check
# ---------------------------------------------------------------------------

INTEGRATE_TOO_FEW_CLUSTERS = 2  # if every resolution yields <= this, fire
INTEGRATE_LARGEST_PCT_TRIGGER = 50.0  # if largest cluster >this % at every res, fire


def integrate_self_check(
    cluster_counts: dict[float, int],
    quality: dict[float, dict[str, Any]] | None,
    resolutions_source: str,
    regress_cell_cycle_already_on: bool,
    n_cells: int = 0,
) -> list[SelfCheckFinding]:
    """
    Two checks:
      1. If every resolution yields <= 2 clusters, suggest --resolutions aio
         (13-step wider sweep). Skip if the user is already on the aio sweep
         — there's no wider one to suggest. Skip when n_cells < 500: a small
         sample with ≤2 populations may genuinely be biologically
         homogeneous, and a wider sweep won't summon clusters that aren't
         in the data.
      2. If the largest cluster at every resolution exceeds 50%, suggest
         --regress-cell-cycle, since the overrepresented cluster is often
         cycling cells. Skip if cell-cycle regression is already on.
    """
    findings: list[SelfCheckFinding] = []

    if cluster_counts:
        max_clusters = max(cluster_counts.values())
        # v0.9.1 (B2): skip the suggestion on small datasets where ≤2 clusters
        # may reflect actual biology (single compartment, low input). The
        # 500-cell threshold mirrors the standard "you need this many to even
        # reasonably cluster" rule of thumb in the scanpy tutorials.
        small_homogeneous = (n_cells > 0 and n_cells < 500)
        if (
            max_clusters <= INTEGRATE_TOO_FEW_CLUSTERS
            and resolutions_source != "aio"
            and not small_homogeneous
        ):
            findings.append(SelfCheckFinding(
                stage="integrate",
                code="integrate_too_few_clusters",
                trigger=(
                    f"every requested resolution yielded <= {INTEGRATE_TOO_FEW_CLUSTERS} clusters "
                    f"(max={max_clusters})"
                ),
                suggestion=(
                    "the resolution sweep tops out at 1.0 — if you expect more substructure, "
                    "switch to --resolutions aio (13-step sweep from 0.01 to 2.0); otherwise "
                    "this dataset may genuinely have ≤2 cell populations"
                ),
                fix={"resolutions": "aio"},
            ))

    if quality and not regress_cell_cycle_already_on:
        largest_pcts = [
            float(m.get("largest_pct", 0.0)) for m in quality.values()
        ]
        if largest_pcts and all(p > INTEGRATE_LARGEST_PCT_TRIGGER for p in largest_pcts):
            min_largest = min(largest_pcts)
            findings.append(SelfCheckFinding(
                stage="integrate",
                code="integrate_dominant_cluster",
                trigger=(
                    f"the largest cluster exceeded {INTEGRATE_LARGEST_PCT_TRIGGER:.0f}% at every "
                    f"resolution (min largest_pct={min_largest:.1f}%); cycling cells often dominate "
                    "this way"
                ),
                suggestion=(
                    "turn on --regress-cell-cycle (Tirosh S/G2M regression) — "
                    "if the dominant cluster is a cell-cycle artifact, regressing it out usually "
                    "uncovers the real biological substructure underneath"
                ),
                fix={"regress_cell_cycle": True},
            ))

    return findings


# ---------------------------------------------------------------------------
# Annotate self-check
# ---------------------------------------------------------------------------

ANNOTATE_AMBIGUOUS_MARGIN = 0.05  # if every cluster's margin is below this, fire
ANNOTATE_IMMUNE_HEURISTIC_FRACTION = 0.4  # v1.1.0: >=this fraction of clusters looking immune fires the panel-swap (was 0.5)
ANNOTATE_IMMUNE_HITS_PER_CLUSTER = 1  # v1.1.0: >=this many immune-marker hits in a cluster's top markers = "looks immune" (was 2)


# Rough immune-marker shortlist for the chondrocyte_markers misuse heuristic.
# Hits any of these in a cluster's top markers is a soft "this looks immune"
# vote; if >=40% of clusters look immune AND we're matching against the
# chondrocyte panel, switch to celltype_broad (which has macrophage / T-cell
# / B-cell groups that will actually score).
#
# v1.1.0 broadens the shortlist (cold-validation gap 2 + ISSUES.md #012):
# the v0.8 list was a textbook PBMC marker set and missed the markers that
# actually show up in subchondral bone / synovium scRNA — HLA class II,
# plasma-cell IGs, mast-cell tryptases, B-cell CD37, dendritic CD1C.
# Concretely, on the BML_1 cold-validation sample the old list found
# immune signal in only 1/13 clusters (NK only); the broadened list catches
# macrophage (HLA-DRA/CD74), plasma (MZB1/JCHAIN/IGHG1), mast (CPA3/TPSAB1),
# B (CD37), and dendritic (CD1C/FCER1A) clusters too. Combined with the
# >=1 hit and >=40% thresholds (lowered from >=2 / >50%), the trigger
# fires comfortably on a 9/13 = 69% immune dataset.
IMMUNE_MARKER_HINTS: tuple[str, ...] = (
    # T cell
    "CD3D", "CD3E", "CD3G", "CD8A", "CD4", "TRAC", "TRBC1", "TRBC2",
    # B cell
    "CD19", "MS4A1", "CD79A", "CD79B", "CD37",
    # Plasma cell
    "MZB1", "JCHAIN", "XBP1", "IGHG1", "IGHA1", "IGHM", "IGKC", "IGLC2", "DERL3",
    # Myeloid (macrophage / monocyte)
    "CD68", "CD163", "LYZ", "CD14", "FCN1", "S100A8", "S100A9",
    "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G",
    # MHC class II (myeloid + DC + B; absent on chondrocytes / fibroblasts)
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "CD74",
    # Dendritic cell (incl. plasmacytoid)
    "CD1C", "FCER1A", "CLEC9A", "IRF7", "IRF8", "IL3RA", "PLD4",
    # Mast cell
    "CPA3", "TPSAB1", "TPSB2", "MS4A2", "CTSG", "KIT", "HPGDS",
    # NK
    "NKG7", "GNLY", "KLRD1", "GZMA", "GZMB", "KLRB1",
    # Pan-leukocyte
    "PTPRC",
)


def annotate_self_check(
    annotations: list[Any],
    panel_name: str,
    profile_module: Any,
) -> list[SelfCheckFinding]:
    """
    Two checks:
      1. If every cluster's panel_margin < 0.05, the chosen panel looks
         under-resolved for this data. Suggest a different --panel (if the
         profile has another) or a different --profile.
      2. If panel_name == "chondrocyte_markers" and most clusters' top
         markers are dominated by immune-cell genes, suggest switching to
         celltype_broad (which has the immune groups defined). This is a
         common misconfiguration: someone runs a synovium / joint-fluid
         dataset against the chondrocyte panel.
    """
    findings: list[SelfCheckFinding] = []
    if not annotations:
        return findings

    margins = [float(getattr(a, "panel_margin", 0.0)) for a in annotations]
    if all(m < ANNOTATE_AMBIGUOUS_MARGIN for m in margins):
        max_margin = max(margins)
        # Find an alternative panel on the profile.
        alt_panel: str | None = None
        for candidate in ("celltype_broad", "chondrocyte_markers"):
            if candidate != panel_name and hasattr(profile_module, candidate):
                alt_panel = candidate
                break

        trigger = (
            f"every cluster's panel_margin < {ANNOTATE_AMBIGUOUS_MARGIN:g} "
            f"(max={max_margin:.3f}); the chosen panel {panel_name!r} is not separating "
            "cell types — labels are ambiguous"
        )
        if alt_panel is not None:
            findings.append(SelfCheckFinding(
                stage="annotate",
                code="annotate_ambiguous_panel",
                trigger=trigger,
                suggestion=(
                    f"switch --panel to {alt_panel!r} — the current panel may be too "
                    "narrow / wrong-tissue for this dataset"
                ),
                fix={"panel_name": alt_panel},
            ))
        else:
            findings.append(SelfCheckFinding(
                stage="annotate",
                code="annotate_ambiguous_no_alt_panel",
                trigger=trigger,
                suggestion=(
                    "switch to a different --profile that ships a panel matching this tissue; "
                    "the current profile has no alternate panel to fall back to"
                ),
                fix={},
            ))

    # Heuristic: chondrocyte panel + dataset that looks immune-dominated.
    if panel_name == "chondrocyte_markers" and hasattr(profile_module, "celltype_broad"):
        n_immune_clusters = 0
        for a in annotations:
            top_markers = [str(g).upper() for g in getattr(a, "top_markers", [])[:30]]
            hits = sum(1 for g in IMMUNE_MARKER_HINTS if g in top_markers)
            # v1.1.0: ANNOTATE_IMMUNE_HITS_PER_CLUSTER (>=1 hit) marks a
            # cluster as "looks immune". Was >=2 in v0.8 — too strict on
            # subchondral bone / synovium where each cluster's top 30
            # often only has 1-3 of the curated markers visible above
            # housekeeping noise. See ISSUES.md #012.
            if hits >= ANNOTATE_IMMUNE_HITS_PER_CLUSTER:
                n_immune_clusters += 1
        immune_fraction = n_immune_clusters / len(annotations)
        # v1.1.0: >=40% (was >50%). On BML_1 cold-validation 9/13 = 69%
        # comfortably clears either bar; the lower threshold catches the
        # boundary case (synovium with mixed chondrocyte + immune lineages)
        # where the old >50% bar silently passed.
        if immune_fraction >= ANNOTATE_IMMUNE_HEURISTIC_FRACTION:
            # Don't double-suggest if we already fired the ambiguous-panel finding
            # against the same alt panel.
            already_suggested_broad = any(
                f.fix.get("panel_name") == "celltype_broad" for f in findings
            )
            if not already_suggested_broad:
                findings.append(SelfCheckFinding(
                    stage="annotate",
                    code="annotate_panel_tissue_mismatch",
                    trigger=(
                        f"{n_immune_clusters}/{len(annotations)} clusters have >=2 immune-cell "
                        "marker hints in their top markers, but the chondrocyte_markers panel "
                        "doesn't define immune groups"
                    ),
                    suggestion=(
                        "switch --panel to celltype_broad — the dataset looks immune-dominated "
                        "(synovium / fluid?) and the chondrocyte panel will assign every cluster "
                        "to a near-zero-score chondrocyte subtype"
                    ),
                    fix={"panel_name": "celltype_broad"},
                ))

    return findings
