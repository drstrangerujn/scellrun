"""
One-shot orchestrator: qc → integrate → markers → annotate → report.

This is the v0.6 entrypoint that chains every single-stage command into
one call. The user gives an .h5ad and a tissue keyword and gets a single
`05_report/index.html` deliverable. New users do not have to know the
stage layout exists.

Design notes:
- The orchestrator reads the input h5ad once, then writes intermediate
  h5ads at each canonical stage path (`<run_dir>/01_qc/qc.h5ad`,
  `02_integrate/integrated.h5ad`, etc.) and re-reads the next stage's
  input from disk. This keeps each stage self-contained and the run-dir
  fully reproducible after the fact (any stage can be re-run from its
  upstream artifact).
- Each stage appends to `00_run.json` via `write_run_meta` with
  `command="analyze:<stage>"` so the manifest captures the orchestrated
  path without losing the per-stage params.
- Errors in one stage abort the pipeline, but stages already completed
  remain on disk; the caller surfaces what got done before the failure.
"""
from __future__ import annotations

import dataclasses
import shutil
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from scellrun.decisions import Decision, record
from scellrun.runlayout import (
    STAGE_DIRS,
    StageOutputExists,
    default_run_dir,
    stage_dir,
    write_run_meta,
)


@dataclass
class AnalyzeResult:
    """Summary of one orchestrated run, returned for tests + CLI prints."""
    run_dir: Path
    stages_completed: list[str]
    qc_summary: str
    integrate_summary: str
    markers_summary: str
    annotate_summary: str
    report_index: Path | None
    chosen_resolution: float
    attempt_id: str = ""


class StageFailure(RuntimeError):
    """Raised when a downstream stage hits an unrecoverable error."""

    def __init__(self, stage: str, original: BaseException):
        super().__init__(f"stage {stage!r} failed: {original}")
        self.stage = stage
        self.original = original


def _pick_first_actionable(findings: list[Any]) -> dict[str, Any] | None:
    """
    Walk a list of SelfCheckFinding and return the first non-empty `fix` dict.

    Findings without a fix dict (the human-only suggestions) are skipped;
    they're already in the decision log for the user to read but the
    orchestrator can't act on them.
    """
    for f in findings:
        fix = getattr(f, "fix", None) or {}
        if fix:
            return dict(fix)
    return None


def _all_finding_codes(findings: list[Any]) -> list[str]:
    """Collect every finding's code in the order seen."""
    return [getattr(f, "code", "") for f in findings if getattr(f, "code", "")]


def _summarize_fix_dict(fix: dict[str, Any]) -> str:
    """Compact stringification of a fix dict for decision-log values."""
    return ", ".join(f"{k}={v}" for k, v in fix.items())


def _move_to_failed(run_dir: Path, stage: str) -> Path | None:
    """
    Rename ``<run_dir>/<NN_stage>/`` to ``<NN_stage>.failed-1/`` and return
    the new path. Used by the --auto-fix retry path so the failed first-pass
    artifacts are preserved next to the retry's fresh stage dir.

    If a `.failed-1` already exists (defensive: somebody re-ran twice in a
    weird state), bump to `.failed-2`, `.failed-3`, … so we never clobber
    older debug payloads.
    """
    sub = STAGE_DIRS.get(stage)
    if sub is None:
        return None
    src = run_dir / sub
    if not src.exists():
        return None
    n = 1
    while True:
        dst = run_dir / f"{sub}.failed-{n}"
        if not dst.exists():
            break
        n += 1
    shutil.move(str(src), str(dst))
    return dst


def _pick_best_resolution(
    quality: dict[float, dict[str, float | int]] | None,
    fallback: float = 0.5,
) -> tuple[float, str]:
    """
    Pick the "best" Leiden resolution from `integrate`'s quality table.

    Returns (resolution, criterion_str) — the criterion string explains
    which branch of the heuristic was used, suitable for the decision-log
    rationale.

    Heuristic:
      1. Prefer "non-fragmented" resolutions (≤1 singleton, ≥2 clusters).
         Among those, take the largest n_clusters; ties broken by smaller
         resolution number (more conservative).
      2. If every resolution is fragmented, prefer the resolution with the
         FEWEST singletons (least fragmented). Tie-break by lowest
         largest_pct (most balanced). This dodges the v0.7 trap where
         every resolution had >1 singleton and the picker fell through to
         max-n_clusters, which on real OA data picked the most fragmented
         resolution as the annotate target. See ISSUES.md #005.
      3. Falls back to `fallback` if quality is empty or no resolution has
         ≥2 clusters.
    """
    if not quality:
        return fallback, "fallback (no quality table)"

    candidates = [
        (res, m["n_clusters"])
        for res, m in quality.items()
        if int(m.get("n_singletons", 0)) <= 1 and int(m.get("n_clusters", 0)) >= 2
    ]
    if candidates:
        max_n = max(int(n) for _, n in candidates)
        best = min(res for res, n in candidates if int(n) == max_n)
        return float(best), "largest n_clusters among non-fragmented (≤1 singleton)"

    # Every resolution is fragmented OR has <2 clusters. Among those with
    # ≥2 clusters, prefer fewest singletons → smallest largest_pct. This
    # gives the user the LEAST fragmented + most balanced resolution
    # available, instead of the most fragmented (= max n_clusters) one.
    any_two_plus = [
        (res, m) for res, m in quality.items() if int(m.get("n_clusters", 0)) >= 2
    ]
    if not any_two_plus:
        return fallback, "fallback (no resolution had ≥2 clusters)"

    def _key(item: tuple[float, dict[str, Any]]) -> tuple[int, float, float]:
        _, m = item
        return (
            int(m.get("n_singletons", 0)),
            float(m.get("largest_pct", 0.0)),
            float(item[0]),
        )

    best_res, _ = min(any_two_plus, key=_key)
    return float(best_res), "fewest singletons → most balanced (every resolution fragmented)"


PANEL_AUTOPICK_CHONDRO_MARGIN = 1.5
"""
v1.1.0: chondrocyte panel must beat broad panel by 1.5x in cluster-level hits
to be kept. A tie or anything weaker swaps to celltype_broad. The pre-1.1.0
"only count broad-only clusters and fire on >50%" heuristic missed the BML_1
case where the chondrocyte panel got 1-2 weak cluster hits and the broad panel
got 7+ cluster hits — that's a tie under the old rule (no broad-only clusters
because each had at least one weak chondrocyte gene), so the chondrocyte panel
won by default. The 1.5x margin requirement matches the cold-validation
journal's "fine-subtype panel only when chondrocyte hits dominate" wording
that the rationale string promised but the code didn't enforce.
"""


def _autopick_panel_for_data(
    profile_module: Any,
    integrated_adata: Any,
    chosen_res: float,
) -> tuple[str, str]:
    """
    Decide which panel to use when both `chondrocyte_markers` and
    `celltype_broad` are defined on the profile.

    v1.1.0 tie-break: count clusters where the chondrocyte panel hits
    (any chondrocyte_markers gene in cluster top markers) and clusters
    where the broad panel hits. Keep chondrocyte_markers ONLY when
    chondrocyte hits >= 1.5x broad hits. Tie or smaller margin → swap
    to celltype_broad. Pre-1.1.0 logic only counted "broad-only" clusters
    and fired on >50%, which silently kept the chondrocyte panel on
    immune-rich subchondral-bone data where every immune cluster picked
    up at least one weak chondrocyte hit (cold-validation gap 1).

    Returns (panel_name, rationale_str).
    """
    has_chondro = hasattr(profile_module, "chondrocyte_markers")
    has_broad = hasattr(profile_module, "celltype_broad")
    if not has_chondro:
        return ("celltype_broad" if has_broad else "", "only celltype_broad available")
    if not has_broad:
        return ("chondrocyte_markers", "only chondrocyte_markers available")

    chondro_panel: dict[str, list[str]] = profile_module.chondrocyte_markers
    broad_panel: dict[str, list[str]] = profile_module.celltype_broad
    chondro_genes = {g.upper() for genes in chondro_panel.values() for g in genes}
    broad_genes = {g.upper() for genes in broad_panel.values() for g in genes}

    # Look up per-cluster top markers via rank_genes_groups on the chosen res.
    # Use scanpy's existing facility; if any of the inputs aren't ready we
    # fall back to the chondrocyte-first preference (current default).
    res_key = f"leiden_res_{chosen_res:g}".replace(".", "_")
    if res_key not in integrated_adata.obs.columns:
        return ("chondrocyte_markers", "default (no leiden col for chosen res)")

    try:
        import scanpy as sc

        rank_key = f"rank_panel_pick_{res_key}"
        if rank_key not in integrated_adata.uns:
            sc.tl.rank_genes_groups(
                integrated_adata,
                groupby=res_key,
                method="wilcoxon",
                use_raw=integrated_adata.raw is not None,
                pts=False,
                key_added=rank_key,
            )
        rg = integrated_adata.uns[rank_key]
        names = rg["names"]
        cluster_names = list(names.dtype.names)
    except Exception:
        return ("chondrocyte_markers", "default (rank_genes_groups failed)")

    n_clusters = len(cluster_names)
    if n_clusters == 0:
        return ("chondrocyte_markers", "default (no clusters)")

    chondro_hits = 0
    broad_hits = 0
    for c in cluster_names:
        top = [str(g).upper() for g in names[c][:30]]
        if any(g in chondro_genes for g in top):
            chondro_hits += 1
        if any(g in broad_genes for g in top):
            broad_hits += 1

    # v1.1.0: chondrocyte panel must clear a 1.5x margin over broad panel
    # to be kept. Tie or smaller margin → broad panel wins. Required margin
    # is documented at PANEL_AUTOPICK_CHONDRO_MARGIN above.
    keeps_chondro = chondro_hits >= PANEL_AUTOPICK_CHONDRO_MARGIN * broad_hits
    if not keeps_chondro:
        return (
            "celltype_broad",
            f"swapped to celltype_broad: chondrocyte_hits={chondro_hits}, "
            f"broad_hits={broad_hits}; required >={PANEL_AUTOPICK_CHONDRO_MARGIN:g}x margin "
            "to keep chondrocyte panel.",
        )
    return (
        "chondrocyte_markers",
        f"kept chondrocyte_markers: chondrocyte_hits={chondro_hits}, "
        f"broad_hits={broad_hits}; cleared the >={PANEL_AUTOPICK_CHONDRO_MARGIN:g}x margin "
        "required to keep chondrocyte panel.",
    )


def _shrink_h5ad(adata: Any) -> None:
    """
    In-place: shrink the integrated h5ad before writing.

    Cast `.X` to float32 (from the float64 that scale() produces) to halve
    on-disk size; if already sparse, leave it alone.
    """
    try:
        import scipy.sparse as sp
    except ImportError:
        sp = None  # type: ignore[assignment]
    X = adata.X
    if sp is not None and sp.issparse(X):
        # Sparse already-compact; convert dtype only if needed.
        if X.dtype != np.float32:
            adata.X = X.astype(np.float32)
        return
    if hasattr(X, "dtype") and X.dtype != np.float32:
        adata.X = X.astype(np.float32, copy=False)


def _shrink_h5ad_for_annotate(adata: Any) -> Any:
    """
    Return a copy of `adata` whose `.X` is the log-normalized matrix
    (from `.raw`) instead of the dense scaled float64 matrix that
    integrate produces. The annotated h5ad is the user-facing artifact;
    they want the log-normalized matrix downstream, not the scaled one.

    If `.raw` isn't set, fall back to shrinking dtype only.
    """
    if adata.raw is not None:
        try:
            new = adata.raw.to_adata()
            # Carry over obs / obsm / obsp / uns from the integrated adata
            # so cluster columns, UMAP, neighbors stay attached.
            for col in adata.obs.columns:
                if col not in new.obs.columns:
                    new.obs[col] = adata.obs[col].values
            for k, v in adata.obsm.items():
                new.obsm[k] = v
            for k, v in adata.obsp.items():
                new.obsp[k] = v
            for k, v in adata.uns.items():
                new.uns[k] = v
            _shrink_h5ad(new)
            return new
        except Exception:
            pass
    _shrink_h5ad(adata)
    return adata


# Marker file lives in the per-user cache dir; first-run hint logic.
_INSTALL_MARKER_PATH = Path.home() / ".cache" / "scellrun" / "installed.touch"


def first_run_hint() -> str | None:
    """
    Return a one-line hint to print on the first invocation of `scellrun`
    in this user's environment. Returns None on subsequent invocations.

    The hint surfaces install heaviness so a non-engineer PI doesn't think
    the install hung. The marker file is `~/.cache/scellrun/installed.touch`
    — once present, the hint goes silent for that user.
    """
    try:
        if _INSTALL_MARKER_PATH.exists():
            return None
        _INSTALL_MARKER_PATH.parent.mkdir(parents=True, exist_ok=True)
        _INSTALL_MARKER_PATH.write_text(
            "scellrun first-run marker; safe to delete to re-trigger the heaviness hint\n",
            encoding="utf-8",
        )
        return (
            "scellrun's deps include scanpy + scrublet + harmonypy + leidenalg; "
            "first install ~3-5 min."
        )
    except OSError:
        # Read-only home / permission issue: skip the hint silently.
        return None


def _load_overrides(path: Path | None) -> dict[str, Any]:
    """Read a review_overrides.json; return {} if no path."""
    if path is None:
        return {}
    import json

    with Path(path).open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"overrides file {path} must contain a JSON object at top level")
    return data


def run_analyze(
    h5ad: Path,
    *,
    profile: str = "default",
    species: str = "human",
    tissue: str | None = None,
    resolutions: tuple[float, ...] | str = "default",
    use_ai: bool = False,
    ai_model: str = "claude-haiku-4-5-20251001",
    lang: str = "en",
    run_dir: Path | None = None,
    force: bool = False,
    max_genes: int | None = None,
    method: str = "harmony",
    method_user_supplied: bool = False,
    regress_cell_cycle: bool = False,
    use_pubmed: bool = False,
    write_h5ad: bool = True,
    auto_fix: bool = False,
    apply_overrides: Path | None = None,
    on_progress: Any = None,
    attempt_id: str | None = None,
) -> AnalyzeResult:
    """
    Run the full QC → integrate → markers → annotate → report pipeline.

    `on_progress`, if provided, is a callable receiving a single string
    summary line per stage; the CLI passes a console.print so the user
    sees per-stage output as it happens.
    """
    import anndata as ad

    from scellrun.profiles import load as load_profile
    from scellrun.report import build_report
    from scellrun.scrna.annotate import run_annotate
    from scellrun.scrna.annotate import write_artifacts as write_annotate_artifacts
    from scellrun.scrna.integrate import (
        AIO_FULL_RESOLUTIONS,
        DEFAULT_RESOLUTIONS,
        run_integrate,
    )
    from scellrun.scrna.integrate import write_artifacts as write_integrate_artifacts
    from scellrun.scrna.markers import run_markers
    from scellrun.scrna.markers import write_artifacts as write_markers_artifacts
    from scellrun.scrna.qc import run_qc
    from scellrun.scrna.qc import write_artifacts as write_qc_artifacts
    from scellrun.views import build_views

    def _say(msg: str) -> None:
        if on_progress is not None:
            on_progress(msg)

    if attempt_id is None:
        attempt_id = uuid.uuid4().hex

    if run_dir is None:
        run_dir = default_run_dir()
    run_dir = Path(run_dir)

    # v1.3.0: load review-server overrides up front so QC threshold
    # tweaks land before the QC stage runs.
    overrides_data = _load_overrides(apply_overrides)
    threshold_overrides_user = (
        overrides_data.get("threshold_overrides") or {}
    ) if overrides_data else {}
    cluster_label_overrides_user = (
        overrides_data.get("cluster_label_overrides") or {}
    ) if overrides_data else {}
    cell_exclusions_user: list[str] = list(
        overrides_data.get("cell_exclusions") or []
    ) if overrides_data else []
    review_notes = str(overrides_data.get("notes") or "") if overrides_data else ""

    if overrides_data:
        record(
            run_dir,
            Decision(
                stage="analyze",
                key="applied_overrides",
                value=str(apply_overrides),
                default=None,
                source="user",
                attempt_id=attempt_id,
                rationale=(
                    f"reviewer-supplied overrides loaded from {apply_overrides}; "
                    f"thresholds={list(threshold_overrides_user.keys()) or 'none'}, "
                    f"label_overrides={len(cluster_label_overrides_user)}, "
                    f"cell_exclusions={len(cell_exclusions_user)}"
                ),
            ),
        )

    # Resolve the resolution sweep up front so we can pass it to integrate.
    if isinstance(resolutions, str):
        if resolutions == "aio":
            res_tuple: tuple[float, ...] = AIO_FULL_RESOLUTIONS
        elif resolutions == "default":
            res_tuple = DEFAULT_RESOLUTIONS
        else:
            res_tuple = tuple(float(x) for x in resolutions.split(",") if x.strip())
    else:
        res_tuple = tuple(resolutions)
    if not res_tuple:
        raise ValueError("resolutions must list at least one value")

    profile_module = load_profile(profile)

    # Early profile sanity check: the annotate stage needs a marker panel.
    # Bail before doing any expensive work if the chosen profile can't
    # satisfy the full chain. Default profile (no tissue) currently has
    # no panel; users wanting the one-shot pipeline should pick a
    # tissue-specific profile (e.g. joint-disease) or supply their own.
    if not (
        hasattr(profile_module, "chondrocyte_markers")
        or hasattr(profile_module, "celltype_broad")
    ):
        raise ValueError(
            f"profile {profile!r} has no marker panels — "
            "scellrun analyze needs a panel at the annotate stage. "
            "Pick a tissue-specific profile (e.g. --profile joint-disease) "
            "or run the per-stage commands and skip annotate."
        )

    # Issue #001: when the default --method is harmony but the input has no
    # sample/batch column to integrate over, downgrade to "none" rather
    # than failing 1m+ into the run. Only auto-downgrade when the caller
    # didn't explicitly pass --method.
    if method == "harmony" and not method_user_supplied:
        try:
            quick = ad.read_h5ad(h5ad, backed="r")
            obs_cols = list(quick.obs.columns)
            quick.file.close()
        except Exception:
            obs_cols = []
        sample_candidates = [
            c for c in ("orig.ident", "sample", "batch", "donor") if c in obs_cols
        ]
        if not sample_candidates:
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="method_downgrade",
                    value="none",
                    default="harmony",
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "no sample/batch column (orig.ident/sample/batch/donor) in obs — "
                        "single-sample input; auto-downgraded --method from harmony to none. "
                        "Pass --method harmony explicitly to force the original behavior."
                    ),
                ),
            )
            method = "none"

    stages_completed: list[str] = []

    # ---- Stage 1: QC -------------------------------------------------------
    qc_out = _resolve_stage_dir(run_dir, "qc", force=force, attempt_id=attempt_id)

    base_thresholds = profile_module.scrna_qc
    overrides: dict[str, Any] = {"species": species}
    if max_genes is not None:
        overrides["max_genes"] = max_genes
    # v1.3.0: review-server threshold overrides land here so the user's
    # slider edits become source="user" rows in the decision log.
    if threshold_overrides_user:
        for k in ("max_pct_mt", "max_genes", "min_counts"):
            if k in threshold_overrides_user and threshold_overrides_user[k] is not None:
                overrides[k] = threshold_overrides_user[k]
    qc_thresholds = dataclasses.replace(base_thresholds, **overrides)
    user_overrides_for_log = {k: v for k, v in overrides.items() if k != "species"}

    try:
        adata = ad.read_h5ad(h5ad)
        qc_result = run_qc(
            adata,
            assay="scrna",
            thresholds=qc_thresholds,
            run_dir=run_dir,
            profile=profile,
            user_thresholds_overrides=user_overrides_for_log,
            profile_applied_thresholds=base_thresholds,
            profile_user_supplied=profile != "default",
            assay_user_supplied=False,
            species_user_supplied=species != "human",
            lang_user_supplied=lang != "en",
            flag_doublets_user_supplied=False,
            attempt_id=attempt_id,
            lang=lang,
        )
        qc_artifacts = write_qc_artifacts(
            qc_result, adata, qc_out, write_h5ad=True, lang=lang
        )
    except Exception as e:
        raise StageFailure("qc", e) from e

    # v0.8 auto-fix: if QC self-check fired with an actionable fix, re-run
    # the QC stage once with the relaxed threshold applied. Cap = 1 retry.
    # v0.9.1: log every fired finding (not just the one we acted on), and
    # preserve the failed-first-pass dir at NN_qc.failed-N/ so the user
    # can compare original vs retry artifacts.
    if auto_fix and qc_result.findings:
        all_codes = _all_finding_codes(qc_result.findings)
        applied = _pick_first_actionable(qc_result.findings)
        skipped = [c for c in all_codes if applied is None or c != all_codes[0]]
        if skipped:
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="qc.skipped_findings",
                    value=skipped,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "self-check fired multiple findings; --auto-fix only acts on "
                        "the first actionable one. The remaining codes are still in "
                        "the trigger/suggest rows above for the user to review."
                    ),
                ),
            )
        if applied is not None:
            new_overrides = dict(overrides)
            new_overrides.update(applied)
            qc_thresholds = dataclasses.replace(base_thresholds, **new_overrides)
            user_overrides_for_log = {
                k: v for k, v in new_overrides.items() if k != "species"
            }
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.qc.applied",
                    value=_summarize_fix_dict(applied),
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    fix_payload=dict(applied),
                    rationale=(
                        f"--auto-fix on; QC self-check suggestion {applied} applied; "
                        "re-running stage once. Failed first-pass artifacts preserved at "
                        "01_qc.failed-1/."
                    ),
                ),
            )
            _say(f"[1/5] qc: --auto-fix applied {applied}; re-running QC")
            failed_dir = _move_to_failed(run_dir, "qc")
            if failed_dir is not None:
                record(
                    run_dir,
                    Decision(
                        stage="analyze",
                        key="auto_fix.qc.failed_first_pass",
                        value=failed_dir.name,
                        default=None,
                        source="auto",
                        attempt_id=attempt_id,
                        rationale=(
                            "preserved the failed first-pass QC artifacts so the user "
                            "can diff original vs retry"
                        ),
                    ),
                )
            try:
                # Auto-fix retry preserves the first-pass self-check rows;
                # truncate_decisions=False keeps them visible in the log.
                qc_out = _resolve_stage_dir(
                    run_dir,
                    "qc",
                    force=True,
                    attempt_id=attempt_id,
                    truncate_decisions=False,
                )
                adata = ad.read_h5ad(h5ad)
                qc_result = run_qc(
                    adata,
                    assay="scrna",
                    thresholds=qc_thresholds,
                    run_dir=run_dir,
                    profile=profile,
                    user_thresholds_overrides=user_overrides_for_log,
                    profile_applied_thresholds=base_thresholds,
                    profile_user_supplied=profile != "default",
                    assay_user_supplied=False,
                    species_user_supplied=species != "human",
                    lang_user_supplied=lang != "en",
                    flag_doublets_user_supplied=False,
                    attempt_id=attempt_id,
                    lang=lang,
                )
                qc_artifacts = write_qc_artifacts(
                    qc_result, adata, qc_out, write_h5ad=True, lang=lang
                )
            except Exception as e:
                raise StageFailure("qc", e) from e
            new_pct = 100.0 * qc_result.n_cells_pass / max(qc_result.n_cells_in, 1)
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.qc.outcome",
                    value=f"{new_pct:.1f}%",
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        f"after auto-fix retry, QC pass-rate is {new_pct:.1f}% "
                        f"({qc_result.n_cells_pass}/{qc_result.n_cells_in}); "
                        + (
                            "retry rescued the stage"
                            if new_pct >= 30.0
                            else "retry did not improve to a passing rate; "
                            "review per_cell_metrics.csv before continuing"
                        )
                    ),
                ),
            )

    write_run_meta(
        run_dir,
        command="analyze:qc",
        params={
            "h5ad": h5ad,
            "out_dir": qc_out,
            "profile": profile,
            "species": species,
            "thresholds": qc_thresholds,
            "force": force,
            "lang": lang,
            "attempt_id": attempt_id,
        },
    )
    qc_pct = 100.0 * qc_result.n_cells_pass / max(qc_result.n_cells_in, 1)
    qc_summary = (
        f"[1/5] qc: {qc_result.n_cells_pass:,} / {qc_result.n_cells_in:,} cells passed "
        f"({qc_pct:.1f}%)"
    )
    _say(qc_summary)
    stages_completed.append("qc")
    qc_h5ad_path = qc_artifacts["qc_h5ad"]

    # ---- Stage 2: integrate ------------------------------------------------
    integ_out = _resolve_stage_dir(
        run_dir, "integrate", force=force, attempt_id=attempt_id
    )

    integrate_res_tuple: tuple[float, ...] = res_tuple
    integrate_regress_cc: bool = regress_cell_cycle

    try:
        qc_adata = ad.read_h5ad(qc_h5ad_path)
        # The orchestrator passes resolutions through verbatim; mark them as
        # "user" if the caller specified anything other than the literal default.
        resolutions_source = (
            "user" if tuple(integrate_res_tuple) != DEFAULT_RESOLUTIONS else "auto"
        )
        if tuple(integrate_res_tuple) == AIO_FULL_RESOLUTIONS:
            resolutions_source = "aio"
        integ_result, integrated = run_integrate(
            qc_adata,
            method=method,
            method_user_supplied=method_user_supplied,
            n_pcs=30,
            resolutions=integrate_res_tuple,
            regress_cell_cycle=integrate_regress_cc,
            regress_cell_cycle_user_supplied=regress_cell_cycle,
            species=species,
            drop_qc_fail=True,
            use_ai=use_ai,
            ai_model=ai_model,
            tissue=tissue,
            run_dir=run_dir,
            sample_key_user_supplied=False,
            resolutions_source=resolutions_source,
            attempt_id=attempt_id,
        )
        # Issue #006: shrink the integrated h5ad before writing.
        _shrink_h5ad(integrated)
        integ_artifacts = write_integrate_artifacts(
            integ_result, integrated, integ_out, write_h5ad=True, lang=lang
        )
    except Exception as e:
        raise StageFailure("integrate", e) from e

    # v0.8 auto-fix: if integrate self-check fired with an actionable fix,
    # re-run the integrate stage once with the fix applied. Cap = 1 retry.
    if auto_fix and integ_result.findings:
        all_codes = _all_finding_codes(integ_result.findings)
        applied = _pick_first_actionable(integ_result.findings)
        skipped = [c for c in all_codes if applied is None or c != all_codes[0]]
        if skipped:
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="integrate.skipped_findings",
                    value=skipped,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "self-check fired multiple findings; --auto-fix only acts on "
                        "the first actionable one. The remaining codes are in the "
                        "trigger/suggest rows above for the user to review."
                    ),
                ),
            )
        if applied is not None:
            if "resolutions" in applied and applied["resolutions"] == "aio":
                integrate_res_tuple = AIO_FULL_RESOLUTIONS
            if "regress_cell_cycle" in applied:
                integrate_regress_cc = bool(applied["regress_cell_cycle"])
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.integrate.applied",
                    value=_summarize_fix_dict(applied),
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    fix_payload=dict(applied),
                    rationale=(
                        f"--auto-fix on; integrate self-check suggestion {applied} applied; "
                        "re-running stage once. Failed first-pass artifacts preserved at "
                        "02_integrate.failed-1/."
                    ),
                ),
            )
            _say(f"[2/5] integrate: --auto-fix applied {applied}; re-running integrate")
            failed_dir = _move_to_failed(run_dir, "integrate")
            if failed_dir is not None:
                record(
                    run_dir,
                    Decision(
                        stage="analyze",
                        key="auto_fix.integrate.failed_first_pass",
                        value=failed_dir.name,
                        default=None,
                        source="auto",
                        attempt_id=attempt_id,
                        rationale=(
                            "preserved the failed first-pass integrate artifacts "
                            "so the user can diff original vs retry"
                        ),
                    ),
                )
            try:
                integ_out = _resolve_stage_dir(
                    run_dir,
                    "integrate",
                    force=True,
                    attempt_id=attempt_id,
                    truncate_decisions=False,
                )
                qc_adata = ad.read_h5ad(qc_h5ad_path)
                resolutions_source = (
                    "aio"
                    if tuple(integrate_res_tuple) == AIO_FULL_RESOLUTIONS
                    else (
                        "user"
                        if tuple(integrate_res_tuple) != DEFAULT_RESOLUTIONS
                        else "auto"
                    )
                )
                integ_result, integrated = run_integrate(
                    qc_adata,
                    method=method,
                    method_user_supplied=method_user_supplied,
                    n_pcs=30,
                    resolutions=integrate_res_tuple,
                    regress_cell_cycle=integrate_regress_cc,
                    regress_cell_cycle_user_supplied=regress_cell_cycle,
                    species=species,
                    drop_qc_fail=True,
                    use_ai=use_ai,
                    ai_model=ai_model,
                    tissue=tissue,
                    run_dir=run_dir,
                    sample_key_user_supplied=False,
                    resolutions_source=resolutions_source,
                    attempt_id=attempt_id,
                )
                _shrink_h5ad(integrated)
                integ_artifacts = write_integrate_artifacts(
                    integ_result, integrated, integ_out, write_h5ad=True, lang=lang
                )
            except Exception as e:
                raise StageFailure("integrate", e) from e
            outcome_text = (
                "retry rescued the stage"
                if (
                    integ_result.cluster_counts
                    and max(integ_result.cluster_counts.values()) > 2
                )
                else "retry did not improve cluster counts; review the cluster sweep"
            )
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.integrate.outcome",
                    value=str(dict(integ_result.cluster_counts)),
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "after auto-fix retry, cluster counts per resolution "
                        f"{dict(integ_result.cluster_counts)} — {outcome_text}"
                    ),
                ),
            )
            res_tuple = integrate_res_tuple
            regress_cell_cycle = integrate_regress_cc

    write_run_meta(
        run_dir,
        command="analyze:integrate",
        params={
            "h5ad": qc_h5ad_path,
            "out_dir": integ_out,
            "method": method,
            "sample_key": integ_result.sample_key,
            "n_pcs": 30,
            "resolutions": list(res_tuple),
            "regress_cell_cycle": regress_cell_cycle,
            "species": species,
            "use_ai": use_ai,
            "tissue": tissue,
            "lang": lang,
            "attempt_id": attempt_id,
        },
    )
    integ_summary = (
        f"[2/5] integrate: {integ_result.n_cells_used:,} cells, method={integ_result.method}, "
        f"clusters per res={dict(integ_result.cluster_counts)}"
    )
    _say(integ_summary)
    stages_completed.append("integrate")
    integrated_h5ad_path = integ_artifacts["integrated_h5ad"]

    # Decide which resolution to feed annotate. Prefer integrate's quality;
    # if AI gave a recommendation and it's parseable, take that.
    chosen_res, criterion = _pick_best_resolution(integ_result.quality)
    chosen_res_source = "auto"
    chosen_res_rationale = criterion
    if integ_result.quality and chosen_res in integ_result.quality:
        m = integ_result.quality[chosen_res]
        chosen_res_rationale += (
            f" — picked res={chosen_res:g}: n_clusters={int(m.get('n_clusters', 0))}, "
            f"largest={float(m.get('largest_pct', 0.0)):.1f}%, "
            f"smallest={float(m.get('smallest_pct', 0.0)):.1f}%, "
            f"singletons={int(m.get('n_singletons', 0))}"
        )
    if integ_result.ai_recommendation is not None:
        try:
            ai_res = float(integ_result.ai_recommendation.get("recommended_resolution", ""))
            if ai_res in res_tuple:
                chosen_res = ai_res
                chosen_res_source = "ai"
                ai_why = integ_result.ai_recommendation.get("rationale") or ""
                chosen_res_rationale = (
                    f"AI recommendation overrode the deterministic pick — {ai_why}"
                ).strip()
        except (TypeError, ValueError):
            pass

    # A2: this decision is never `source="user"` regardless of whether the
    # value happens to differ from a hard-coded default. Use Decision
    # directly with explicit `source` so the AI branch can override it.
    record(
        run_dir,
        Decision(
            stage="analyze",
            key="chosen_resolution_for_annotate",
            value=chosen_res,
            default=None,
            source=chosen_res_source,
            attempt_id=attempt_id,
            rationale=chosen_res_rationale,
        ),
    )

    # ---- Stage 3: markers --------------------------------------------------
    markers_out = _resolve_stage_dir(
        run_dir, "markers", force=force, attempt_id=attempt_id
    )

    try:
        markers_adata = ad.read_h5ad(integrated_h5ad_path)
        markers_result, per_res_df = run_markers(
            markers_adata, run_dir=run_dir, attempt_id=attempt_id
        )
        write_markers_artifacts(
            markers_result, per_res_df, markers_out, lang=lang
        )
    except Exception as e:
        raise StageFailure("markers", e) from e

    write_run_meta(
        run_dir,
        command="analyze:markers",
        params={
            "h5ad": integrated_h5ad_path,
            "out_dir": markers_out,
            "resolutions": list(markers_result.resolutions),
            "logfc_threshold": markers_result.logfc_threshold,
            "pct_min": markers_result.pct_min,
            "only_positive": markers_result.only_positive,
            "top_n": markers_result.top_n,
            "lang": lang,
            "attempt_id": attempt_id,
        },
    )
    n_markers_total = sum(markers_result.n_markers_per_resolution.values())
    markers_summary = (
        f"[3/5] markers: {n_markers_total:,} markers across "
        f"{len(markers_result.resolutions)} resolutions"
    )
    _say(markers_summary)
    stages_completed.append("markers")

    # ---- Stage 4: annotate -------------------------------------------------
    annot_out = _resolve_stage_dir(
        run_dir, "annotate", force=force, attempt_id=attempt_id
    )

    # If chosen_res is not actually in the integrated h5ad, fall back to 0.5
    # (or the first available resolution).
    available = list(markers_result.resolutions)
    if chosen_res not in available:
        if 0.5 in available:
            chosen_res = 0.5
        elif available:
            chosen_res = available[0]

    # Issue #003 part 2: when the profile defines both panels, peek at the
    # data's top markers to decide which panel to use, instead of always
    # preferring chondrocyte_markers. Falls back to the legacy preference
    # if the data doesn't make the choice clear.
    annot_panel_name: str | None = None
    try:
        peek_adata = ad.read_h5ad(integrated_h5ad_path)
        picked_panel, panel_pick_rationale = _autopick_panel_for_data(
            profile_module, peek_adata, chosen_res
        )
        if picked_panel:
            annot_panel_name = picked_panel
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="annotate.auto_panel",
                    value=annot_panel_name,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=panel_pick_rationale,
                ),
            )
    except Exception:
        # Auto-pick is best-effort; on failure fall through to annotate's
        # built-in _select_panel which prefers chondrocyte_markers.
        annot_panel_name = None

    try:
        annot_adata = ad.read_h5ad(integrated_h5ad_path)
        annot_result = run_annotate(
            annot_adata,
            profile_module,
            resolution=chosen_res,
            panel_name=annot_panel_name,
            use_ai=use_ai,
            ai_model=ai_model,
            use_pubmed=use_pubmed,
            tissue=tissue,
            run_dir=run_dir,
            profile_user_supplied=profile != "default",
            resolution_user_supplied=False,  # orchestrator picked it
            # v1.1.0: orchestrator-injected panel name is auto, not user.
            # Even when annot_panel_name is a non-None string, it came
            # from the auto-pick step or a self-check fix payload, NOT
            # from the user passing --panel on the CLI. Cold-validation gap 3.
            panel_name_user_supplied=False,
            use_ai_user_supplied=False,
            use_pubmed_user_supplied=use_pubmed,
            tissue_user_supplied=tissue is not None,
            attempt_id=attempt_id,
        )
        annot_adata = _shrink_h5ad_for_annotate(annot_adata)
        write_annotate_artifacts(
            annot_result, annot_adata, annot_out, write_h5ad=write_h5ad, lang=lang
        )
    except Exception as e:
        raise StageFailure("annotate", e) from e

    # v0.8 auto-fix: if annotate self-check fired with an actionable fix
    # (e.g. switch panel), re-run annotate once with the new panel.
    if auto_fix and annot_result.findings:
        all_codes = _all_finding_codes(annot_result.findings)
        applied = _pick_first_actionable(annot_result.findings)
        skipped = [c for c in all_codes if applied is None or c != all_codes[0]]
        if skipped:
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="annotate.skipped_findings",
                    value=skipped,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "self-check fired multiple findings; --auto-fix only acts on "
                        "the first actionable one. The remaining codes are in the "
                        "trigger/suggest rows above for the user to review."
                    ),
                ),
            )
        if applied is not None and "panel_name" in applied:
            annot_panel_name = str(applied["panel_name"])
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.annotate.applied",
                    value=_summarize_fix_dict(applied),
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    fix_payload=dict(applied),
                    rationale=(
                        f"--auto-fix on; annotate self-check suggestion {applied} applied; "
                        "re-running stage once. Failed first-pass artifacts preserved at "
                        "04_annotate.failed-1/."
                    ),
                ),
            )
            _say(
                f"[4/5] annotate: --auto-fix applied {applied}; re-running annotate "
                f"with panel={annot_panel_name}"
            )
            failed_dir = _move_to_failed(run_dir, "annotate")
            if failed_dir is not None:
                record(
                    run_dir,
                    Decision(
                        stage="analyze",
                        key="auto_fix.annotate.failed_first_pass",
                        value=failed_dir.name,
                        default=None,
                        source="auto",
                        attempt_id=attempt_id,
                        rationale=(
                            "preserved the failed first-pass annotate artifacts "
                            "so the user can diff original vs retry"
                        ),
                    ),
                )
            try:
                annot_out = _resolve_stage_dir(
                    run_dir,
                    "annotate",
                    force=True,
                    attempt_id=attempt_id,
                    truncate_decisions=False,
                )
                annot_adata = ad.read_h5ad(integrated_h5ad_path)
                annot_result = run_annotate(
                    annot_adata,
                    profile_module,
                    resolution=chosen_res,
                    panel_name=annot_panel_name,
                    use_ai=use_ai,
                    ai_model=ai_model,
                    use_pubmed=use_pubmed,
                    tissue=tissue,
                    run_dir=run_dir,
                    profile_user_supplied=profile != "default",
                    resolution_user_supplied=False,
                    # v1.1.0: see note on first run_annotate call above —
                    # this retry's panel comes from a self-check fix payload,
                    # still orchestrator-driven, not user-driven.
                    panel_name_user_supplied=False,
                    use_ai_user_supplied=False,
                    use_pubmed_user_supplied=use_pubmed,
                    tissue_user_supplied=tissue is not None,
                    attempt_id=attempt_id,
                )
                annot_adata = _shrink_h5ad_for_annotate(annot_adata)
                write_annotate_artifacts(
                    annot_result, annot_adata, annot_out, write_h5ad=write_h5ad, lang=lang
                )
            except Exception as e:
                raise StageFailure("annotate", e) from e
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.annotate.outcome",
                    value=annot_result.panel_name,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        "after auto-fix retry, annotate ran with panel "
                        f"{annot_result.panel_name!r}"
                    ),
                ),
            )

    write_run_meta(
        run_dir,
        command="analyze:annotate",
        params={
            "h5ad": integrated_h5ad_path,
            "out_dir": annot_out,
            "profile": profile,
            "panel_name": annot_result.panel_name,
            "resolution": chosen_res,
            "use_ai": use_ai,
            "use_pubmed": use_pubmed,
            "tissue": tissue,
            "lang": lang,
            "attempt_id": attempt_id,
        },
    )
    annotate_summary = (
        f"[4/5] annotate: {len(annot_result.annotations)} clusters labeled at "
        f"res={chosen_res:g} via panel {annot_result.panel_name}"
    )
    _say(annotate_summary)
    stages_completed.append("annotate")

    # ---- v1.3.0: apply review overrides -----------------------------------
    # Cluster-label and cell-exclusion overrides land AFTER the annotate
    # stage so the deterministic call still ran (and is in the log) and
    # the user's edit becomes a source="user" override on top.
    if cluster_label_overrides_user or cell_exclusions_user:
        _apply_post_annotate_overrides(
            run_dir=run_dir,
            annot_out=annot_out,
            integrated_h5ad_path=integrated_h5ad_path,
            cluster_label_overrides=cluster_label_overrides_user,
            cell_exclusions=cell_exclusions_user,
            attempt_id=attempt_id,
        )

    # ---- Stage 5: report ---------------------------------------------------
    report_out = _resolve_stage_dir(
        run_dir, "report", force=force, attempt_id=attempt_id
    )

    try:
        report_artifacts = build_report(run_dir, report_out, lang=lang)
    except Exception as e:
        raise StageFailure("report", e) from e

    write_run_meta(
        run_dir,
        command="analyze:report",
        params={"out_dir": report_out, "lang": lang, "attempt_id": attempt_id},
    )
    report_index = report_artifacts["index"]
    report_summary = f"[5/5] report: {report_index}"
    _say(report_summary)
    stages_completed.append("report")

    # ---- v1.2.0: views layer (best-effort; pure HTML index over the
    # existing stage artifacts). Failures here don't break the pipeline
    # because the canonical 05_report/index.html is already on disk.
    try:
        views_artifacts = build_views(run_dir)
    except Exception as e:
        views_artifacts = {}
        _say(f"[5/5] views: skipped — {type(e).__name__}: {e}")
    else:
        write_run_meta(
            run_dir,
            command="analyze:views",
            params={
                "out_dir": run_dir / STAGE_DIRS["views"],
                "n_views": len(views_artifacts),
                "attempt_id": attempt_id,
            },
        )
        if "index" in views_artifacts:
            _say(f"[5/5] views: {views_artifacts['index']}")

    if review_notes:
        # Top-level "review_notes" field on the run manifest so any reader
        # of 00_run.json sees the reviewer's reasoning without parsing
        # the per-stage params.
        import json as _json

        meta_path = run_dir / "00_run.json"
        if meta_path.exists():
            try:
                meta = _json.loads(meta_path.read_text(encoding="utf-8"))
                meta["review_notes"] = review_notes
                meta_path.write_text(
                    _json.dumps(meta, indent=2, default=str), encoding="utf-8"
                )
            except (OSError, _json.JSONDecodeError):
                pass

    return AnalyzeResult(
        run_dir=run_dir,
        stages_completed=stages_completed,
        qc_summary=qc_summary,
        integrate_summary=integ_summary,
        markers_summary=markers_summary,
        annotate_summary=annotate_summary,
        report_index=report_index,
        chosen_resolution=chosen_res,
        attempt_id=attempt_id,
    )


def _resolve_stage_dir(
    run_dir: Path,
    stage: str,
    *,
    force: bool,
    attempt_id: str,
    truncate_decisions: bool | None = None,
) -> Path:
    """
    Stage-dir resolution that handles auto-resume on incomplete prior runs.

    Issue #010: if a previous run aborted mid-pipeline, the stage dir
    exists but the stage hasn't completed (no `report.html` for stages
    that produce one, or no `index.html` for the report stage). We
    overwrite that incomplete dir without making the user pass --force.
    A `<stage>.auto_resume` decision row records the implicit overwrite.

    For complete prior runs (report.html present), behavior is unchanged:
    StageOutputExists fires and the caller wraps it in StageFailure.

    `truncate_decisions` defaults to `force` — pass False explicitly when
    an internal retry (auto-fix) wants to preserve the trigger/suggest
    rows from the first pass.
    """
    if truncate_decisions is None:
        truncate_decisions = force
    sub = run_dir / STAGE_DIRS[stage]
    if not force and sub.exists():
        sentinel = sub / ("index.html" if stage == "report" else "report.html")
        if sub.exists() and not sentinel.exists():
            # Incomplete prior run; auto-resume into this stage dir.
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key=f"{stage}.auto_resume",
                    value=sub.name,
                    default=None,
                    source="auto",
                    attempt_id=attempt_id,
                    rationale=(
                        f"prior run left {sub.name}/ without a completed report; "
                        "overwriting in place rather than erroring on StageOutputExists"
                    ),
                ),
            )
            try:
                return stage_dir(
                    run_dir, stage, force=True, truncate_decisions=truncate_decisions
                )
            except StageOutputExists as e:
                raise StageFailure(stage, e) from None
    try:
        return stage_dir(
            run_dir, stage, force=force, truncate_decisions=truncate_decisions
        )
    except StageOutputExists as e:
        raise StageFailure(stage, e) from None


def _apply_post_annotate_overrides(
    *,
    run_dir: Path,
    annot_out: Path,
    integrated_h5ad_path: Path,
    cluster_label_overrides: dict[str, str],
    cell_exclusions: list[str],
    attempt_id: str,
) -> None:
    """
    v1.3.0: rewrite 04_annotate/annotations.csv with the reviewer's
    label overrides applied (deterministic call still recorded under
    `panel_label`; new `final_label` column carries the user's pick),
    and patch 04_annotate/annotated.h5ad's `obs.scellrun_qc_pass` to
    flag out user-excluded barcodes.

    Both edits are recorded as `source="user"` rows in the decision log
    so the report carries the override trail. Excluded cells stay in
    the .h5ad — they're flagged out, not removed (the QC stage's
    flag-don't-drop policy extends here too).
    """
    import csv as _csv

    annotations_csv = annot_out / "annotations.csv"
    if cluster_label_overrides and annotations_csv.exists():
        rows: list[dict[str, str]] = []
        with annotations_csv.open("r", encoding="utf-8", newline="") as f:
            reader = _csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
            for r in reader:
                rows.append(dict(r))
        if "final_label" not in fieldnames:
            fieldnames.append("final_label")
        if "label_source" not in fieldnames:
            fieldnames.append("label_source")
        for r in rows:
            cid = str(r.get("cluster", ""))
            if cid in cluster_label_overrides:
                r["final_label"] = cluster_label_overrides[cid]
                r["label_source"] = "user"
            else:
                r["final_label"] = r.get("panel_label", "")
                r["label_source"] = "auto"
        with annotations_csv.open("w", encoding="utf-8", newline="") as f:
            writer = _csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                writer.writerow(r)

        for cid, label in cluster_label_overrides.items():
            record(
                run_dir,
                Decision(
                    stage="annotate",
                    key=f"label_override.cluster_{cid}",
                    value=label,
                    default=None,
                    source="user",
                    attempt_id=attempt_id,
                    rationale=(
                        f"reviewer overrode cluster {cid} label to {label!r} via "
                        "scellrun review (final_label column)"
                    ),
                ),
            )

    if cell_exclusions:
        # Patch the annotated.h5ad obs in place — the deterministic
        # qc_pass row still exists in 01_qc/qc.h5ad; the annotate-stage
        # h5ad gets the user-flagged-out cells.
        h5_path = annot_out / "annotated.h5ad"
        if h5_path.exists():
            try:
                import anndata as ad

                a = ad.read_h5ad(h5_path)
                if "scellrun_qc_pass" not in a.obs.columns:
                    a.obs["scellrun_qc_pass"] = True
                excl_set = set(cell_exclusions)
                # Match by exact obs_names; barcodes the user supplied
                # but not present in the h5ad are silently no-ops (we
                # log the count below).
                mask = a.obs_names.isin(excl_set)
                n_hit = int(mask.sum())
                a.obs.loc[mask, "scellrun_qc_pass"] = False
                a.write_h5ad(h5_path)
            except Exception:
                n_hit = 0
        else:
            n_hit = 0
        record(
            run_dir,
            Decision(
                stage="annotate",
                key="cell_exclusions",
                value=len(cell_exclusions),
                default=None,
                source="user",
                attempt_id=attempt_id,
                rationale=(
                    f"reviewer marked {len(cell_exclusions)} barcodes for exclusion; "
                    f"matched {n_hit} in annotated.h5ad — flagged scellrun_qc_pass=False"
                ),
            ),
        )


# Re-export for tests / callers that imported os from this module historically.
__all__ = [
    "AnalyzeResult",
    "StageFailure",
    "_apply_post_annotate_overrides",
    "_load_overrides",
    "_pick_best_resolution",
    "first_run_hint",
    "run_analyze",
]
