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
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from scellrun.decisions import Decision, record
from scellrun.runlayout import (
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


def _summarize_fix_dict(fix: dict[str, Any]) -> str:
    """Compact stringification of a fix dict for decision-log values."""
    return ", ".join(f"{k}={v}" for k, v in fix.items())


def _pick_best_resolution(
    quality: dict[float, dict[str, float | int]] | None,
    fallback: float = 0.5,
) -> float:
    """
    Pick the "best" Leiden resolution from `integrate`'s quality table.

    Heuristic: pick the largest n_clusters among non-fragmented resolutions,
    where "non-fragmented" means at most one singleton cluster (clusters
    holding < 2% of cells). On ties, prefer the smaller resolution number
    (more conservative). Falls back to `fallback` if quality is empty or
    no resolution is non-fragmented.
    """
    if not quality:
        return fallback

    candidates = [
        (res, m["n_clusters"])
        for res, m in quality.items()
        if int(m.get("n_singletons", 0)) <= 1 and int(m.get("n_clusters", 0)) >= 2
    ]
    if not candidates:
        # Every resolution is fragmented. Try "any with >=2 clusters", else fallback.
        any_two_plus = [
            (res, m["n_clusters"])
            for res, m in quality.items()
            if int(m.get("n_clusters", 0)) >= 2
        ]
        if not any_two_plus:
            return fallback
        candidates = any_two_plus

    # max n_clusters, ties broken by smallest resolution number
    max_n = max(int(n) for _, n in candidates)
    best = min(res for res, n in candidates if int(n) == max_n)
    return float(best)


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
    regress_cell_cycle: bool = False,
    use_pubmed: bool = False,
    write_h5ad: bool = True,
    auto_fix: bool = False,
    on_progress: Any = None,
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

    def _say(msg: str) -> None:
        if on_progress is not None:
            on_progress(msg)

    if run_dir is None:
        run_dir = default_run_dir()
    run_dir = Path(run_dir)

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

    stages_completed: list[str] = []

    # ---- Stage 1: QC -------------------------------------------------------
    try:
        qc_out = stage_dir(run_dir, "qc", force=force)
    except StageOutputExists as e:
        raise StageFailure("qc", e) from None

    base_thresholds = profile_module.scrna_qc
    overrides: dict[str, Any] = {"species": species}
    if max_genes is not None:
        overrides["max_genes"] = max_genes
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
            lang=lang,
        )
        qc_artifacts = write_qc_artifacts(
            qc_result, adata, qc_out, write_h5ad=True, lang=lang
        )
    except Exception as e:
        raise StageFailure("qc", e) from e

    # v0.8 auto-fix: if QC self-check fired with an actionable fix, re-run
    # the QC stage once with the relaxed threshold applied. Retry capped at
    # 1 to avoid loops; if the second pass-rate still fails, log and continue.
    if auto_fix and qc_result.findings:
        applied = _pick_first_actionable(qc_result.findings)
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
                    rationale=(
                        f"--auto-fix on; QC self-check suggestion {applied} applied; "
                        "re-running stage once"
                    ),
                ),
            )
            _say(f"[1/5] qc: --auto-fix applied {applied}; re-running QC")
            try:
                qc_out = stage_dir(run_dir, "qc", force=True)
                adata = ad.read_h5ad(h5ad)
                qc_result = run_qc(
                    adata,
                    assay="scrna",
                    thresholds=qc_thresholds,
                    run_dir=run_dir,
                    profile=profile,
                    user_thresholds_overrides=user_overrides_for_log,
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
                    rationale=(
                        f"after auto-fix retry, QC pass-rate is {new_pct:.1f}% "
                        f"({qc_result.n_cells_pass}/{qc_result.n_cells_in})"
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
    try:
        integ_out = stage_dir(run_dir, "integrate", force=force)
    except StageOutputExists as e:
        raise StageFailure("integrate", e) from None

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
            n_pcs=30,
            resolutions=integrate_res_tuple,
            regress_cell_cycle=integrate_regress_cc,
            species=species,
            drop_qc_fail=True,
            use_ai=use_ai,
            ai_model=ai_model,
            tissue=tissue,
            run_dir=run_dir,
            sample_key_user_supplied=False,
            resolutions_source=resolutions_source,
        )
        integ_artifacts = write_integrate_artifacts(
            integ_result, integrated, integ_out, write_h5ad=True, lang=lang
        )
    except Exception as e:
        raise StageFailure("integrate", e) from e

    # v0.8 auto-fix: if integrate self-check fired with an actionable fix,
    # re-run the integrate stage once with the fix applied. Cap = 1 retry.
    if auto_fix and integ_result.findings:
        applied = _pick_first_actionable(integ_result.findings)
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
                    rationale=(
                        f"--auto-fix on; integrate self-check suggestion {applied} applied; "
                        "re-running stage once"
                    ),
                ),
            )
            _say(f"[2/5] integrate: --auto-fix applied {applied}; re-running integrate")
            try:
                integ_out = stage_dir(run_dir, "integrate", force=True)
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
                    n_pcs=30,
                    resolutions=integrate_res_tuple,
                    regress_cell_cycle=integrate_regress_cc,
                    species=species,
                    drop_qc_fail=True,
                    use_ai=use_ai,
                    ai_model=ai_model,
                    tissue=tissue,
                    run_dir=run_dir,
                    sample_key_user_supplied=False,
                    resolutions_source=resolutions_source,
                )
                integ_artifacts = write_integrate_artifacts(
                    integ_result, integrated, integ_out, write_h5ad=True, lang=lang
                )
            except Exception as e:
                raise StageFailure("integrate", e) from e
            record(
                run_dir,
                Decision(
                    stage="analyze",
                    key="auto_fix.integrate.outcome",
                    value=str(dict(integ_result.cluster_counts)),
                    default=None,
                    source="auto",
                    rationale=(
                        "after auto-fix retry, cluster counts per resolution "
                        f"{dict(integ_result.cluster_counts)}"
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
    chosen_res = _pick_best_resolution(integ_result.quality)
    chosen_res_source = "auto"
    chosen_res_rationale = (
        "largest n_clusters among non-fragmented resolutions "
        "(<=1 singleton, >=2 clusters); ties broken by smallest resolution"
    )
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

    record(
        run_dir,
        Decision(
            stage="analyze",
            key="chosen_resolution_for_annotate",
            value=chosen_res,
            default=None,
            source=chosen_res_source,
            rationale=chosen_res_rationale,
        ),
    )

    # ---- Stage 3: markers --------------------------------------------------
    try:
        markers_out = stage_dir(run_dir, "markers", force=force)
    except StageOutputExists as e:
        raise StageFailure("markers", e) from None

    try:
        markers_adata = ad.read_h5ad(integrated_h5ad_path)
        markers_result, per_res_df = run_markers(markers_adata, run_dir=run_dir)
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
    try:
        annot_out = stage_dir(run_dir, "annotate", force=force)
    except StageOutputExists as e:
        raise StageFailure("annotate", e) from None

    # If chosen_res is not actually in the integrated h5ad, fall back to 0.5
    # (or the first available resolution).
    available = list(markers_result.resolutions)
    if chosen_res not in available:
        if 0.5 in available:
            chosen_res = 0.5
        elif available:
            chosen_res = available[0]

    annot_panel_name: str | None = None

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
        )
        write_annotate_artifacts(
            annot_result, annot_adata, annot_out, write_h5ad=write_h5ad, lang=lang
        )
    except Exception as e:
        raise StageFailure("annotate", e) from e

    # v0.8 auto-fix: if annotate self-check fired with an actionable fix
    # (e.g. switch panel), re-run annotate once with the new panel.
    if auto_fix and annot_result.findings:
        applied = _pick_first_actionable(annot_result.findings)
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
                    rationale=(
                        f"--auto-fix on; annotate self-check suggestion {applied} applied; "
                        "re-running stage once"
                    ),
                ),
            )
            _say(
                f"[4/5] annotate: --auto-fix applied {applied}; re-running annotate "
                f"with panel={annot_panel_name}"
            )
            try:
                annot_out = stage_dir(run_dir, "annotate", force=True)
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
                )
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
        },
    )
    annotate_summary = (
        f"[4/5] annotate: {len(annot_result.annotations)} clusters labeled at "
        f"res={chosen_res:g} via panel {annot_result.panel_name}"
    )
    _say(annotate_summary)
    stages_completed.append("annotate")

    # ---- Stage 5: report ---------------------------------------------------
    try:
        report_out = stage_dir(run_dir, "report", force=force)
    except StageOutputExists as e:
        raise StageFailure("report", e) from None

    try:
        report_artifacts = build_report(run_dir, report_out, lang=lang)
    except Exception as e:
        raise StageFailure("report", e) from e

    write_run_meta(
        run_dir,
        command="analyze:report",
        params={"out_dir": report_out, "lang": lang},
    )
    report_index = report_artifacts["index"]
    report_summary = f"[5/5] report: {report_index}"
    _say(report_summary)
    stages_completed.append("report")

    return AnalyzeResult(
        run_dir=run_dir,
        stages_completed=stages_completed,
        qc_summary=qc_summary,
        integrate_summary=integ_summary,
        markers_summary=markers_summary,
        annotate_summary=annotate_summary,
        report_index=report_index,
        chosen_resolution=chosen_res,
    )
