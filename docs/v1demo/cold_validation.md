# Cold-agent SKILL.md validation — v1.0.1

Goal: a fresh agent uses ONLY SKILL.md to handle the user's prompt.
Did SKILL.md provide enough information? Where did it fall short?

User prompt:
> Look at my OA cartilage scRNA data on the lab server. The cellranger
> output is at `~/Gtest/R_data/南方医科大软骨下骨/seurat/BML_1/`.

Notable wrinkle the prompt adds beyond the canonical example:
the path literally contains the substring "软骨下骨" (subchondral bone),
even though the user said "OA cartilage". SKILL.md addresses this case
(joint-disease panel auto-pick can swap chondrocyte → broad), but a
fresh agent has to notice the mismatch.

## My read of SKILL.md (1-2 sentences per section)

- **Verification status** — pinned to v1.0.0, with `docs/agent-demo.md`
  as the canonical reference; if my run disagrees on the same input, I
  should suspect SKILL.md is stale, not myself.
- **Decision tree** — `analyze` is default; for a cellranger / mtx dir,
  the path is `scrna convert` → `analyze`. So this prompt = convert + analyze.
- **End-to-end one-liner** — `analyze h5ad --profile joint-disease
  --tissue "OA cartilage" --ai --auto-fix`; the report at
  `<run>/05_report/index.html` is the canonical handoff.
- **`--method` and single-sample auto-degrade** — harmony auto-falls
  back to none if no sample/batch column; recorded as `auto` /
  `method_downgrade` row. Don't second-guess, just read the log.
- **Reading the decision log** — `00_decisions.jsonl`, schema_version=1,
  group by stage, translate `source` (auto/user/ai), highlight
  `value != default`, quote `rationale` verbatim.
- **Self-check findings** — paired `.trigger` / `.suggest` rows;
  always surface the trigger before any data findings; `fix_payload`
  is structured form for tools, `value` is human prose; ask before
  re-running unless user pre-authorized.
- **`.failed-N/` directories** — old artifacts preserved next to retry's
  output; check `auto_fix.<stage>.outcome` row before claiming retry
  rescued anything.
- **Trigger codes table** — qc_low_pass_rate, integrate_too_few_clusters,
  integrate_dominant_cluster, annotate_ambiguous_panel,
  annotate_panel_tissue_mismatch, etc., each with structured fix.
- **Profile selection** — `joint-disease` for any joint tissue keyword
  (cartilage, synovium, subchondral bone, OA, RA). Don't default to
  `default` if a tissue-specific profile exists.
- **Panel auto-selection inside joint-disease** — chondrocyte_markers
  vs celltype_broad auto-picked from cluster top markers; subchondral-
  bone / synovium with mostly immune clusters → celltype_broad. Safe
  to use joint-disease even on subchondral bone.
- **Finding the data** — use the path the user gave; multi-sample is
  not auto-merged; `analyze` takes one h5ad.
- **Environment hygiene** — per-user conda env, python 3.11, `pip
  install scellrun`, confirm `scellrun --version`.
- **Hard rules** — don't write own scanpy code, don't silently drop
  cells, don't override profile defaults silently, don't assume raw
  counts, surface overrides + decision-log rationale, surface
  self-check before claiming findings.
- **Failure modes** — concrete recipes for the QC / cluster /
  annotation triggers.
- **End-to-end example dialogue** — the closest thing to a worked
  example; mirrors the structure of my expected reply.
- **Anti-patterns** — explicit list of "looks right but isn't" traps.
- **Source of truth & provenance** — repo links and where the defaults
  come from (Liu lab AIO + clinical practice).

## What I did

(chronological, hospital-server commands stripped of ssh wrapper for readability)

1. `ls ~/Gtest/R_data/南方医科大软骨下骨/seurat/BML_1/` to see what's
   actually in the cellranger directory. Found standard 10x trio:
   `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` (~72MB matrix).
   Note: the parent directory name contains 软骨下骨 = subchondral bone.
   The user said "OA cartilage" but the path says subchondral bone —
   I'll pass `--tissue "OA cartilage"` per the user's words but flag
   to the user that the panel auto-pick may swap to celltype_broad if
   the data is immune-rich.
2. `conda env list` — confirmed there is already a `scellrun-v1demo`
   env (which I was told to NOT touch) and various others. No
   `scellrun-coldval` env yet.
3. `conda create -n scellrun-coldval python=3.11 -y`. Clean env per
   SKILL.md "Environment hygiene". Per-task name (`coldval`).
4. `pip install scellrun==1.0.1` — **failed**. PyPI only has up to
   `1.0.0`. The validation brief says install `1.0.1` but that
   version isn't on PyPI yet (and the brief explicitly says don't
   tag/push). SKILL.md itself is "last validated against scellrun
   v1.0.0", so installing 1.0.0 is consistent with SKILL.md's pin.
5. `pip install scellrun==1.0.0`. Slow (PyPI international link).
6. While install ran, listed
   `~/Gtest/R_data/南方医科大软骨下骨/seurat/` and noticed 9 sibling
   sample dirs: `BML_1..5`, `NC_1..3`, `NC_SCB_1`. The user only asked
   about `BML_1`, so I'll do BML_1 first and surface the multi-sample
   question in my reply (SKILL.md "Finding the data" → "scellrun does
   **not** auto-merge multiple cellranger outputs").

## What I did (continued from step 7)

(Resumed in a new session after the Anthropic limit reset. Confirmed
the prior agent's `scellrun-coldval` env still exists with `scellrun
1.0.0` installed — no need to recreate. v1.0.2 is now on PyPI but
I stayed on 1.0.0 for resumption fidelity, as instructed.)

7. **Decision: BML_1 alone, surface siblings.** SKILL.md "Finding
   the data" is explicit — scellrun does not auto-merge cellranger
   outputs and `analyze` takes one h5ad. The user named BML_1
   specifically; if they want a multi-sample analysis I'll need to
   ask whether to merge BML_* (5 OA), NC_* (3 normal cartilage), and
   NC_SCB_1 (1 normal subchondral bone) and which sample column
   names they want. I'll mention this in the reply but not auto-do
   it.

8. **Convert.** `scellrun scrna convert . -o ~/scellrun_coldval/BML_1.h5ad`
   from the BML_1 dir. Auto-detected the 10x trio. Output:

   ```
   scellrun scrna convert  •  format=auto
   reading .
   loaded: 12,451 cells × 33,538 genes
   wrote: /home/data/liuxiyu/scellrun_coldval/BML_1.h5ad
   ```

   Clean. No warnings. SKILL.md's one-liner pattern worked verbatim.

9. **Analyze.**

   ```bash
   scellrun analyze BML_1.h5ad \
       --profile joint-disease \
       --tissue "OA cartilage" \
       --auto-fix \
       --no-ai
   ```

   `--no-ai` because the hospital server has no `ANTHROPIC_API_KEY`
   exported. SKILL.md's anti-pattern list flags exactly this case
   ("--ai is on by default, just leave it" → no, only when key is
   set). Stdout (5/5 stages):

   ```
   scellrun analyze  •  profile=joint-disease  species=human
                       tissue=OA cartilage  ai=False
   [1/5] qc: 10,835 / 12,451 cells passed (87.0%)
   [2/5] integrate: 10,835 cells, method=none, clusters per
         res={0.1: 10, 0.3: 13, 0.5: 19, 0.8: 23, 1.0: 25}
   [3/5] markers: 48,790 markers across 5 resolutions
   [4/5] annotate: 13 clusters labeled at res=0.3 via panel
         chondrocyte_markers
   [5/5] report: scellrun_out/run-20260430-183440/05_report/index.html
   ```

   ~6 minutes wall. Deluge of `PerformanceWarning: DataFrame is
   highly fragmented` from scanpy's `rank_genes_groups` during the
   markers stage — cosmetic, but a dozen+ identical warnings fly past
   and could panic a clinician reading the terminal. SKILL.md does
   not mention this output noise.

10. **Decision log.** `00_decisions.jsonl` has 27 rows, schema_version=1.
    Highlights I would surface:
    - `analyze.method_downgrade` — `value=none`, `default=harmony`,
      `source=auto`, rationale: "no sample/batch column
      (orig.ident/sample/batch/donor) in obs — single-sample input;
      auto-downgraded --method from harmony to none." Matches
      SKILL.md "single-sample auto-degrade" verbatim.
    - `qc.profile` — `value=joint-disease`, `source=user`. Override
      flagged amber in the HTML.
    - `qc.max_pct_mt` — `value=20.0`, rationale verbatim: "mt%
      ceiling 20.0% — joint tissue is stress-prone, the textbook
      10% silently drops real chondrocytes (PI cohort 2024-2026,
      AIO PM=20)". This is exactly the rationale SKILL.md's
      anti-patterns section calls out.
    - `analyze.chosen_resolution_for_annotate` — `value=0.3`,
      rationale: "fewest singletons → most balanced (every
      resolution fragmented) — picked res=0.3: n_clusters=13,
      largest=31.5%, smallest=0.2%, singletons=2"
    - `analyze.annotate.auto_panel` — `value=chondrocyte_markers`,
      rationale: "fine-subtype panel preferred (chondrocyte hits
      dominate or are tied)"
    - **No `self_check.*` rows. No `auto_fix.*` rows.** Despite
      what's clearly a panel-tissue mismatch (see step 11), no
      trigger fired.

11. **Annotation table** (`04_annotate/annotations.csv`, my read):

    | clus | size | panel_label | top markers | what they actually are |
    | --- | --- | --- | --- | --- |
    | 0 | 31.5% | preInfC | IFI27 GNG11 RAMP2 A2M SPARCL1 ENG TM4SF1 HSPG2 NPDC1 IGFBP7 | endothelial (RAMP2, ENG, A2M) |
    | 1 | ~17% | InfC | LYZ AIF1 HLA-DRA HLA-DRB1 HLA-DPA1 TYROBP FCER1G CD74 | macrophage / monocyte |
    | 2 | — | Unassigned | RPS27 IL32 RPS29 CXCR4 B2M CCL5 BTG1 RPS3 RPS15A ZFP36L2 | T cell / lymphocyte |
    | 3 | — | Unassigned | NKG7 KLRD1 CTSW GZMA CST7 KLRB1 HCST CD7 CCL5 CD247 | NK / CD8 T |
    | 4 | — | FC (1.0/1.0) | COL1A2 DCN LUM COL1A1 COL3A1 COL6A2 COL6A1 PCOLCE NNMT MGP | fibroblast (panel called FC = fibrocartilage chondrocyte? overloaded) |
    | 5 | — | ProC | IGFBP7 CALD1 NDUFA4L2 ADIRF RGS5 TPM2 MYL9 TAGLN NOTCH3 MCAM | pericyte / mural cell |
    | 6 | — | Unassigned | CD37 RPS8 CD79A RPL13A RPS23 RPLP2 RPS27 RPL11 | B cell |
    | 7 | — | InfC (0.33) | GZMB PLD4 JCHAIN PLAC8 SEC61B IRF8 IL3RA C12orf75 GPR183 IRF7 | plasmacytoid DC |
    | 8 | — | InfC (0.33) | ACTB TUBA1B H2AFZ GSTP1 ACTG1 HMGN2 TMSB10 JPT1 HLA-DPA1 CD74 | proliferating immune (housekeeping + HLA) |
    | 9 | — | HTC (0.33) | ACP5 MMP9 CTSK RGS10 TCIRG1 SPP1 CST3 CD68 ATP6V0E1 | osteoclast (ACP5+CTSK+MMP9 trio) |
    | 10 | — | Unassigned | KLRB1 IL7R LTB BTG1 IL32 RPLP1 CXCR4 SPOCK2 TNFAIP3 CD69 | memory T / ILC |
    | 11 | — | Unassigned | MZB1 SSR4 DERL3 JCHAIN SEC11C XBP1 FKBP11 IGHG1 CD79A PRDX4 | plasma cell |
    | 12 | — | Unassigned | CPA3 TPSAB1 TPSB2 HPGD HPGDS MS4A2 ANXA1 CTSG FTH1 UBB | mast cell |

    The data is overwhelmingly **bone marrow / immune-rich**, not
    cartilage. This matches the directory name (软骨下骨 =
    subchondral bone), and BML almost certainly stands for "bone
    marrow lesion" — a known OA imaging finding. The user said "OA
    cartilage" but the sample is sBM/BML aspirate. Six clusters are
    Unassigned (panel doesn't know NK/B/T/mast/plasma); the InfC
    calls on clusters 1, 7, 8 are wrong (cluster 1 is myeloid, 7 is
    pDC, 8 is generic activated immune); cluster 9 (osteoclast)
    being called HTC is wrong; clusters 0 and 5 (endothelial,
    pericyte) get chondrocyte-panel labels they shouldn't.

    **Per SKILL.md's auto-panel rule**, when "more than 50% of
    clusters fall in the celltype_broad-only bucket, scellrun swaps
    the panel to celltype_broad". Here >50% of clusters are clearly
    immune-only (T/NK/B/mast/plasma/pDC/macrophage = 7+/13). The
    swap did NOT happen. The auto-pick rationale "fine-subtype
    panel preferred (chondrocyte hits dominate or are tied)" is
    factually wrong on this dataset — chondrocyte hits do not
    dominate. The panel kept matching anything (endothelial got
    "preInfC", pericyte got "ProC") rather than reporting "no
    chondrocyte hits".

12. **Self-check missed it.** SKILL.md trigger table promises
    `annotate_panel_tissue_mismatch` for "chondrocyte panel chosen
    but most clusters have immune-marker top genes". That is
    exactly this run. No trigger row, no suggest row, no auto_fix
    row. The pipeline finished "clean" with `--auto-fix` on, never
    second-guessing the panel.

## Where SKILL.md got me through cleanly

- The decision tree said "10x mtx → `scrna convert` then `analyze`".
  Worked verbatim, no improvising.
- The end-to-end one-liner template was correct down to flag
  spelling. I needed to swap `--ai` for `--no-ai` on the hospital
  server, which SKILL.md's flag description anticipated ("Pass
  `--no-ai` explicitly if the user does not want LLM calls").
- "Single-sample auto-degrade" section explained the
  `method_downgrade` row before I saw it. No surprise, no scramble.
- "Reading the decision log" with the `source` translation table
  (auto / user / ai) made it trivial to walk the 27 rows.
- The hard-rule list ("don't write your own scanpy", "don't silently
  filter", "quote rationale verbatim") matches what I'd want to do
  anyway, and the rationale strings in `00_decisions.jsonl` ARE
  quotable verbatim — they're written for the senior reviewer
  audience SKILL.md promises.
- "Finding the data" → "does not auto-merge" was the right
  guidance for the multi-sample sibling situation.

## Where SKILL.md was INSUFFICIENT

- **scanpy `PerformanceWarning: DataFrame is highly fragmented`
  spam** is unmentioned. A dozen+ identical warnings fire during
  the markers stage on a real dataset. SKILL.md should either say
  "ignore these, cosmetic" or push for a `warnings.filterwarnings`
  inside `scrna markers`.
- **The `FC` and `HTC` panel labels are not glossed.** A clinician
  asking "what does HTC mean?" — SKILL.md only names "HomC", "RegC",
  "preFC", "InfC" inline in the example dialogue, and "Fan 2024
  11-subtype panel" without listing all 11. I had to infer "HTC =
  hypertrophic chondrocyte", "FC = fibrocartilage chondrocyte",
  "ProC = proliferating chondrocyte" from context, but the panel
  also called osteoclasts as HTC and pericytes as ProC. A glossary
  table or a pointer to `profiles show joint-disease` (which I did
  not run, sticking to "read SKILL.md once") would help.
- **Cluster sizes are not in the annotations.csv.** The CSV gives
  cluster ID but no cell count. I only saw cluster 0's size (31.5%)
  via the decision log's res-pick rationale. SKILL.md's example
  dialogue cites "% of cells" per cluster as if it were a column —
  it isn't. The agent has to cross-reference `02_integrate/`
  artifacts.
- **No `--no-ai` worked-example stdout in SKILL.md.** Hospital
  servers without API keys are common. Sample output for the
  `--no-ai` path (where `use_ai=false` and AI second-opinion is
  skipped) would help me know what's expected vs missing.
- **What to do when self-check should have fired but didn't.**
  SKILL.md assumes self-check is the safety net, but here it
  silently passed a clearly-mis-paneled run. There's no guidance
  for "the auto-pick made an obviously wrong choice and the
  self-check didn't catch it" — am I supposed to override with
  `--panel celltype_broad` and re-run? Manually inspect markers and
  flag to user? I went with "flag to user, recommend re-run with
  `--panel celltype_broad`", but that's a guess.
- **"BML" sample naming.** SKILL.md has no rules for interpreting
  Chinese path names or sample IDs. The 软骨下骨/BML clue was the
  key tell that this is bone marrow lesion aspirate (subchondral
  bone marrow), not cartilage. The user said "cartilage". A
  clinician-facing tool should probably warn when the panel-data
  mismatch is THIS extreme.

## Where SKILL.md was AMBIGUOUS

- **"More than 50% of clusters in celltype_broad-only bucket → swap
  to celltype_broad"** vs the actual rationale "chondrocyte hits
  dominate or are tied". My data has near-zero true chondrocyte
  signal — but the auto-pick said hits "dominate or are tied".
  How is "tied" computed? If 0 chondrocyte hits = 0 broad hits =
  tie → defaults to chondrocyte_markers, that's a footgun. SKILL.md
  doesn't define the tie-break. (I suspect the v1.0.0 implementation
  has a tie-default-to-chondrocyte that the doc doesn't reflect.)
- **"`source=user` for `panel`"** in the decision log. The agent
  (me) did not pass `--panel`. The auto-pick step records itself
  as `source=user, rationale="--panel 'chondrocyte_markers' forced"`.
  That's misleading — no human or agent forced anything; the
  orchestrator selected it auto-pick-style. SKILL.md says
  `source=user` means "the caller passed a flag different from
  default". This row breaks that contract. (See "Where SKILL.md was
  WRONG" — the orchestrator-passes-flag-to-stage-CLI artifact.)
- **`--auto-fix` did nothing.** SKILL.md says auto-fix re-runs a
  stage if a self-check fires. No self-check fired. So `--auto-fix`
  was a no-op. SKILL.md should be explicit: `--auto-fix` is only
  useful if a trigger codes table entry actually fires. Passing it
  blindly is fine but does nothing on a "clean" run.
- **What does `panel_score=0.5, panel_margin=0.5` mean for cluster
  0?** When panel_score == panel_margin, every other panel got 0.
  But cluster 4 has panel_score=1.0/margin=1.0 (FC). SKILL.md
  describes "margin to runner-up" but doesn't pin down what it
  means when only one panel entry hits at all (margin = score).

## Where SKILL.md was WRONG

- **The auto-panel-swap promise does not hold for v1.0.0 on
  immune-rich subchondral-bone data.** SKILL.md says "scellrun will
  swap panels rather than blindly run the chondrocyte panel" and
  "subchondral-bone with mostly immune clusters → celltype_broad".
  This run is exactly that case (BML = bone marrow lesion aspirate)
  and the swap did NOT happen. Either the auto-pick threshold is
  miscalibrated (50% of clusters bucketed wrong), or the
  "tie/dominate" comparison is off. **This is the most important
  gap.** SKILL.md is selling a safety property the binary does not
  deliver here.
- **`annotate_panel_tissue_mismatch` did not fire.** SKILL.md's
  trigger codes table lists this exact code with structured fix
  `{"panel_name": "celltype_broad"}`. Conditions clearly met (most
  clusters have immune-marker top genes; chondrocyte panel chosen).
  Either it's not actually wired in v1.0.0, or its threshold is
  stricter than SKILL.md describes.
- **`source=user` semantics on the annotate panel row.** SKILL.md:
  "`user` = override (caller passed a flag different from
  default)". v1.0.0 actual: when the orchestrator runs `analyze`,
  the auto-pick chooses a panel name and passes it as `--panel
  chondrocyte_markers` to the per-stage CLI, which logs it as
  `source=user`. This is an internal-vs-external distinction
  SKILL.md doesn't acknowledge. Agents quoting "this was an
  override" to the user (per SKILL.md hard rule 5) will be wrong —
  the user did NOT override it; the orchestrator did its
  auto-pick.
- **v1.0.2 README rewrite (not yet verified).** Brief mentions
  v1.0.2 README rewrite. I did not read README.md per cold-read
  rules and did not upgrade. SKILL.md is pinned to v1.0.0; if
  v1.0.2 fixed the auto-panel-swap or the self-check trigger, the
  pin needs bumping. (Out of scope for this validation but flagged.)

## My final reply to the user (drafted)

> I ssh'd to the lab server, set up a per-user conda env
> (`scellrun-coldval`, scellrun 1.0.0), converted the BML_1
> cellranger output to h5ad, and ran the joint-disease pipeline
> end-to-end. **The numbers are clean but the cell-type labels are
> off** — I need to flag this before walking the report.
>
> **Run summary** (`~/scellrun_coldval/scellrun_out/run-20260430-183440/05_report/index.html`):
>
> - **QC:** 10,835 / 12,451 cells passed (87.0%). pct_mt under the
>   joint-tissue 20% cap. Doublets flagged (not dropped) per scellrun
>   policy. No QC self-check triggers.
> - **Integrate:** Single-sample, so `--method` auto-downgraded from
>   harmony to none — quoting the decision log: *"no sample/batch
>   column (orig.ident/sample/batch/donor) in obs — single-sample
>   input; auto-downgraded --method from harmony to none."* If you
>   want harmony across samples, see the multi-sample question
>   below.
> - **Resolution:** scellrun picked res=0.3 (13 clusters), rationale:
>   *"fewest singletons → most balanced (every resolution
>   fragmented) — picked res=0.3: n_clusters=13, largest=31.5%,
>   smallest=0.2%, singletons=2."*
> - **Annotation:** ⚠ panel mismatch. The auto-pick chose the Fan
>   2024 chondrocyte panel ("fine-subtype panel preferred —
>   chondrocyte hits dominate or are tied"), but only 1 of 13
>   clusters has a real chondrocyte signature (cluster 4: COL1A2,
>   DCN, LUM, COL1A1 — fibroblast/fibrocartilage). The other
>   12 clusters are bone-marrow / immune cells: macrophages
>   (cluster 1, LYZ+CD74+TYROBP), NK/T (cluster 3, NKG7+GZMA+KLRD1),
>   B cells (cluster 6, CD79A), plasma cells (cluster 11, MZB1+
>   IGHG1+JCHAIN), mast cells (cluster 12, TPSAB1+CPA3+MS4A2),
>   pDC (cluster 7, IRF7+IRF8+IL3RA), osteoclasts (cluster 9,
>   ACP5+CTSK+MMP9+CD68), endothelial (cluster 0, RAMP2+ENG+A2M),
>   pericytes (cluster 5, NOTCH3+MCAM+RGS5), and a few unresolved.
>   Six are correctly flagged "Unassigned"; the others got
>   chondrocyte labels they don't deserve (cluster 9 osteoclast →
>   HTC, cluster 1 macrophage → InfC, cluster 0 endothelial →
>   preInfC, cluster 5 pericyte → ProC).
>
> **What's going on:** the path is `软骨下骨/seurat/BML_1`. 软骨下骨
> = subchondral bone, BML almost certainly = bone marrow lesion. The
> sample is bone-marrow aspirate from an OA subchondral lesion, not
> cartilage. That's why it's immune-rich.
>
> **What to do:**
>
> 1. **Re-run annotate with the broad panel.** scellrun ships a
>    `celltype_broad` panel on this profile (15-group: chondrocytes,
>    fibroblasts, endothelial, pericytes, macrophages, T/B/NK/plasma,
>    etc.) which will label these clusters correctly. One command:
>
>    ```
>    scellrun scrna annotate \
>        scellrun_out/run-20260430-183440/02_integrate/integrated.h5ad \
>        --profile joint-disease --panel celltype_broad \
>        --resolution 0.3 --tissue "subchondral bone marrow" \
>        --force \
>        -o scellrun_out/run-20260430-183440/04_annotate_broad
>    ```
>
> 2. **Multi-sample question.** I noticed 9 sibling sample dirs in
>    the same parent: BML_1..5 (5 OA samples), NC_1..3 (3 controls),
>    NC_SCB_1 (1 normal SCB). scellrun does not auto-merge cellranger
>    outputs — if you want a unified analysis with Harmony across
>    samples, I'd merge them upstream into one h5ad with an
>    `orig.ident` / `sample` column, then re-run `analyze`. Want me
>    to do that next? If so, which samples in / out and what's the
>    OA-vs-NC contrast you care about?
> 3. **Tissue-label.** Going forward, calling this dataset
>    "subchondral bone" or "BML aspirate" rather than "OA cartilage"
>    will let the auto-pick land on `celltype_broad` from the start
>    (and let the LLM resolution recommender, when --ai is on, weigh
>    immune-rich subtypes correctly).
>
> Decision log full table is in the report. Happy to drill into any
> cluster, or to re-run with `--panel celltype_broad` + multi-sample
> merge if you confirm.

## Verdict

Validation does NOT fully pass. SKILL.md got me through every
mechanical step (env, convert, analyze, log read, report read) and
the rationale strings in the decision log were quotable verbatim,
which is the doc's main promise. But the **auto-panel-swap safety
net SKILL.md sells did not engage on a textbook case for it**
(immune-rich subchondral bone marrow data + joint-disease profile),
and the **`annotate_panel_tissue_mismatch` self-check that's listed
in the trigger table did not fire**. A fresh agent following only
SKILL.md would either (a) trust the panel labels and tell the
clinician their bone marrow sample is full of chondrocytes — a
real clinical mis-call — or (b) notice the markers don't match and
go off-script. I caught it in (b) by reading the markers
myself, but SKILL.md never tells me to do that as a sanity check;
it tells me to trust the panel + self-check pipeline.

**Recommendation: hold for SKILL.md improvements.** Three concrete
gaps to close before tagging:

1. Either fix the v1.0.x auto-pick threshold so it actually swaps
   to celltype_broad on immune-rich data (and document the swap
   triggered in the rationale), OR add a SKILL.md section "the
   auto-pick is fallible; always sanity-check the top markers
   against the panel labels before quoting them to the user".
2. Document the `source=user` artifact: when `analyze` orchestrates
   per-stage CLIs and passes auto-picked values via flags, the row
   logs `source=user` even though no human user picked it. SKILL.md
   hard rule 5 ("when value differs from default, tell the user
   this was an override") will produce false claims otherwise.
3. Add a glossary of the chondrocyte_markers labels (HomC, RegC,
   HTC, FC, ProC, InfC, preInfC, preFC, etc.) — agents currently
   have to infer from context, and the labels can be wrong-tissue
   when the panel auto-pick fails.

