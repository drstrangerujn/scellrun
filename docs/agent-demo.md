# scellrun + agent: end-to-end demo on real OA cartilage data

This is what an LLM agent (Claude Code / Hermes / Codex) does when a
clinician asks for help with single-cell data. The conversation is
real; the data is the PI's actual OA cohort sample (BML_1, the same
sample dogfooded against v0.7); the agent reasoning quotes scellrun's
decision log verbatim.

The full run-dir lives on the hospital server at
`/tmp/scellrun_v1demo/scellrun_out/run-20260430-151851/`. Its
`05_report/index.html` is committed under
[`v1demo/index.html`](v1demo/index.html); the decision log under
[`v1demo/decisions.jsonl`](v1demo/decisions.jsonl); the run manifest
under [`v1demo/run.json`](v1demo/run.json).

## The conversation

> **User:** Look at my OA cartilage scRNA data on the lab server. The
> cellranger output is at
> `~/Gtest/R_data/南方医科大软骨下骨/seurat/BML_1/`.

> **Agent:** I'll set up a clean conda env on the hospital server, install
> scellrun, convert the cellranger directory to an h5ad, and run
> `scellrun analyze` end-to-end with the joint-disease profile (because
> "OA cartilage" → joint tissue → tighter hb% ceiling and the Fan 2024
> chondrocyte panel come pre-loaded). I'll keep `--no-ai` since
> there's no Anthropic key on this server, and `--auto-fix` so the
> orchestrator can apply the cheapest stage-self-check suggestion if
> something looks off. Output will land in `/tmp/scellrun_v1demo/`.

```bash
ssh hospital 'bash -lc "
  source ~/miniconda3/etc/profile.d/conda.sh
  conda create -n scellrun-v1demo python=3.11 -y
  conda activate scellrun-v1demo
  pip install --upgrade scellrun
  scellrun --version    # 0.9.1
"'

ssh hospital 'bash -lc "
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate scellrun-v1demo
  rm -rf /tmp/scellrun_v1demo && mkdir -p /tmp/scellrun_v1demo
  cd /tmp/scellrun_v1demo
  scellrun scrna convert ~/Gtest/R_data/南方医科大软骨下骨/seurat/BML_1 -o data.h5ad
"'
# loaded: 12,451 cells × 33,538 genes
# wrote: data.h5ad

ssh hospital 'bash -lc "
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate scellrun-v1demo
  cd /tmp/scellrun_v1demo
  scellrun analyze data.h5ad \
      --profile joint-disease \
      --tissue \"OA cartilage\" \
      --no-ai \
      --auto-fix \
      --lang en 2>&1 | tee v1demo.log
"'
```

> **Agent:** Pipeline finished in ~6 minutes. Here's what scellrun saw:
>
> **QC** (stage 1/5) flagged 10,835 of 12,451 cells as passing
> (87.0%). Doublet flag count was 0 — known scrublet auto-threshold
> issue logged as ISSUES.md #004 from the v0.7 dogfood; the score
> distribution itself was sane. Pass-rate is well above the 60%
> self-check trigger so no mt% relaxation was needed.
>
> **Integrate** (stage 2/5) ran with `method=none`. The orchestrator
> auto-downgraded from the default `harmony` because the obs has no
> sample/batch column — single-sample input. This is logged as a
> `source="auto"` decision row with the rationale verbatim:
> *"no sample/batch column (orig.ident/sample/batch/donor) in obs —
> single-sample input; auto-downgraded --method from harmony to none.
> Pass --method harmony explicitly to force the original behavior."*
> Cluster counts swept across resolutions: `{0.1: 10, 0.3: 13, 0.5:
> 19, 0.8: 23, 1.0: 25}`.
>
> **Markers** (stage 3/5) found 48,790 differentially expressed genes
> across the five resolutions (positive markers, log2fc≥1, pct≥0.25).
>
> **Annotate** (stage 4/5) labeled 13 clusters at res=0.3 against the
> `chondrocyte_markers` panel. Six clusters scored 0 against every
> panel group and were honestly labeled `Unassigned` (the v0.9.1 fix
> for ISSUES.md #003). Top labeled clusters: cluster 4 — `FC` (fibro
> chondrocyte) at score 1.0, cluster 1 — `InfC` at 0.67, cluster 0 —
> `preInfC` at 0.5.
>
> **Report** (stage 5/5) at
> [`v1demo/index.html`](v1demo/index.html). The first block is "At a
> glance"; below it, a per-stage decision table with `auto` / `user`
> badges and rationales. 27 decision rows total.
>
> Anything you want to drill into?

> **User:** Why did you pick resolution 0.3? The integrate output had
> resolutions 0.1 through 1.0 and you picked the second-lowest.

> **Agent:** Decision-log row for `chosen_resolution_for_annotate`,
> verbatim:
>
> ```
> stage:       analyze
> key:         chosen_resolution_for_annotate
> value:       0.3
> source:      auto
> rationale:   fewest singletons → most balanced (every resolution
>              fragmented) — picked res=0.3: n_clusters=13,
>              largest=31.5%, smallest=0.2%, singletons=2
> ```
>
> The picker has two branches. Preferred branch: pick the resolution
> with the largest n_clusters that still has ≤1 singleton cluster
> (defined as <0.5% of cells). On this sample every resolution was
> "fragmented" — even res=0.1 had 2 singletons — so the picker fell
> through to its v0.9.1 fallback: among resolutions with ≥2 clusters,
> pick the one with the fewest singletons; tie-break on smallest
> largest_pct. res=0.3 had 2 singletons (lowest, tied with 0.1); 0.3
> had 13 clusters vs 10 at res=0.1, so 0.3 wins. The picker explicitly
> avoids the legacy "max n_clusters" rule that made the v0.7 run pick
> res=1.0 with 7 singletons — see ISSUES.md #005.

> **User:** What's cluster 5? It's labeled "ProC" in the annotate
> table — that doesn't look right for an OA dataset.

> **Agent:** Pulling the annotate row for cluster 5 verbatim from
> [`v1demo/index.html`](v1demo/index.html) → `04_annotate/annotations.csv`:
>
> ```
> cluster:       5
> panel_label:   ProC
> panel_score:   0.333
> panel_margin:  0.333
> ai_label:      (AI not used; ANTHROPIC_API_KEY not set on hospital server)
> top_markers:   IGFBP7, CALD1, NDUFA4L2, ADIRF, RGS5, TPM2, MYL9,
>                TAGLN, NOTCH3, MCAM
> ```
>
> The deterministic call is wrong here. The top markers are textbook
> pericyte / mural-cell signature: RGS5, MCAM, NOTCH3, MYL9, TAGLN,
> TPM2 — that's the Bondjers / pericyte panel almost gene-for-gene.
> The chondrocyte_markers panel matched 1/3 of the `ProC` group's
> three genes (`HMGA1` — generic high-mobility-group; not specific),
> margin 0.33 — non-zero but uninformative.
>
> The right call is to switch to the `celltype_broad` panel, which has
> a `Pericytes: [RGS5, MCAM, ACTA2]` group, or to use `--ai` (the LLM
> reads the top markers without the panel-fit constraint and
> typically picks pericyte directly). Three of the 13 clusters look
> like this — chondrocyte_markers labels them with low non-zero
> scores when celltype_broad would label them honestly:
>
> | cluster | top markers (truncated) | panel_label | what it actually is |
> | --- | --- | --- | --- |
> | 5 | RGS5, MCAM, NOTCH3, MYL9, TAGLN, TPM2 | ProC (0.33) | Pericyte |
> | 7 | GZMB, PLD4, JCHAIN, IRF8, IL3RA, IRF7 | InfC (0.33) | Plasmacytoid DC |
> | 9 | ACP5, MMP9, CTSK, RGS10, TCIRG1, SPP1, CD68 | HTC (0.33) | Osteoclast |
>
> The other six clusters honestly say `Unassigned` (score 0) because
> the chondrocyte panel has no group that matches NK / T / B / plasma
> / mast / monocyte signatures.
>
> The v0.8 self-check that would catch this — `annotate_panel_tissue_mismatch`
> — counts immune-marker hints in each cluster's top genes against a
> shortlist (CD3D, NKG7, CD79A, LYZ, …). On this sample only 1 of 13
> clusters cleared the ≥2-hits-per-cluster bar, so the heuristic
> didn't fire even though most clusters are clearly non-chondrocyte.
> The shortlist is too narrow — it misses HLA-DR\*, MZB1, CPA3, CD37
> — so the auto-panel-pick stayed on chondrocyte_markers. Logged as
> ISSUES.md #012.

> **User:** Can you re-run with the celltype_broad panel? I want to
> see the immune labels.

> **Agent:** Two ways. The cheapest is to re-run only the annotate
> stage against the existing run-dir, which preserves QC + integrate +
> markers and only re-labels:
>
> ```bash
> ssh hospital 'bash -lc "
>   source ~/miniconda3/etc/profile.d/conda.sh
>   conda activate scellrun-v1demo
>   cd /tmp/scellrun_v1demo
>   scellrun scrna annotate \
>       scellrun_out/run-20260430-151851/02_integrate/integrated.h5ad \
>       --profile joint-disease \
>       --panel celltype_broad \
>       --resolution 0.3 \
>       --tissue \"OA cartilage\" \
>       --run-dir scellrun_out/run-20260430-151851 \
>       --force
> "'
> ```
>
> This rewrites `04_annotate/annotations.csv` with celltype_broad
> (Chondrocytes / Fibroblasts / Endothelial cells / Pericytes /
> Macrophages / Monocytes / T cells / B cells / Plasma cells / NK
> cells / DC / Osteoblasts / Osteoclasts / SMC / MSC). The decision
> log appends a new `panel = celltype_broad (user)` row alongside the
> earlier auto-picked one — both stay in the file, both get an
> `attempt_id` so the report renderer knows which annotate-pass each
> row belongs to. Then `scellrun report scellrun_out/run-20260430-151851`
> regenerates `index.html`.
>
> Alternative — re-run the full pipeline with `--panel celltype_broad`
> from the start. Costs the integrate step again (~3 min), so I
> wouldn't unless you also want to bump resolution or change QC. Say
> the word and I'll do it.

## What this proves

1. The auto-degrade for single-sample input lands cleanly (no
   harmony-without-sample-key crash); it's logged with rationale, not
   silently swallowed.
2. The marker tables show real signal at the top of every cluster
   (CLEC14A/RAMP3/EMCN for endothelial, RGS5/MCAM/NOTCH3 for
   pericyte) — no log2fc=21 ribosomal noise. The v0.9.1 fix to sort
   by `pval_adj` then `log2fc` is the difference.
3. The annotate stage labels honestly: clusters with no panel match
   say `Unassigned` rather than "ProC at score 0.0". That alone
   makes the report safe to hand a clinician.
4. The decision log + at-a-glance block on the index page give the
   user one screen to read and the answer to "why mt% 20?" / "why
   res 0.3?" in plain text — without the user having to read scanpy
   docs or the source.
5. Where the system is still wrong (cluster 5 / 7 / 9 mislabeled by
   chondrocyte panel; immune-heuristic shortlist too narrow), the
   wrongness is visible in the artifacts and recoverable in one
   re-run command. The agent does not have to lie about what it sees.
