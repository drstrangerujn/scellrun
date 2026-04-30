# Differentiation research — what's the case for scellrun vs vanilla LLM agent + scanpy?

PI's challenge: a modern LLM agent told to "analyze this scRNA-seq" can already write scanpy boilerplate and produce a report. Where's our edge?

## Quote-worthy field statements (use in README)

### From "Perspectives on rigor and reproducibility in single cell genomics" (PMC9122178)

- **"in my group's experience, it is not unusual for reanalysis to find 20% fewer or more clusters in datasets"** — same raw data, different analyst, different answer. This is the failure scellrun fixes by encoding the choices in code instead of leaving them implicit.
- **Of ~50 single-cell papers in high-impact journals, "just a handful" reported external validation.** This is the "agent improvises a notebook → no audit trail" problem at scale.
- Recommended remedy (the field's own consensus): "transparent pipeline documentation, reporting reproducibility metrics, restricting downstream analysis to 'core cells' that consistently cluster across permutations."

### From the nf-core scrnaseq pipeline (the standard others reach for)

- **Differentiator: containers, parameter files, modular processes, automated CI on full-sized AWS datasets.**
- The case they make: "Custom workflows often rely on ad hoc approaches, hindering reproducibility, limiting configurability and reporting, and complicating integration with broader bioinformatics workflows. Adapting to new data, parameters, or references can be slow, and inconsistent reporting and poor modularity can limit scalability."
- This is the same argument we make, but at a different layer. nf-core fixes "different machines" reproducibility; scellrun fixes "different analysts (or different agents)" reproducibility.

### CellAtria (npj Artificial Intelligence, 2025) — direct comparable

- An "agentic AI framework for ingestion and standardization of single-cell RNA-seq data analysis."
- Markets itself on "fully transparent, auditable, reliable analytical steps" (basically our pitch).
- Confirms the audit/transparency framing isn't ours alone — the field is converging here.

## What scellrun actually does that vanilla agent + scanpy doesn't

1. **Decision log (`00_decisions.jsonl`).** Every threshold, every panel, every resolution recorded with `value`, `default`, `source` (auto/user/ai), and a `rationale` string. Reviewer asks "why mt% 20?" — answer is in the file, defensible. Vanilla agent's notebook has no such record.

2. **Same data → same answer.** scellrun is deterministic: run twice, get the same clusters and same labels. Two LLM agents handed the same h5ad will produce two different mt% thresholds, two different HVG counts, two different cluster counts (cf. the "20% more or fewer clusters on reanalysis" finding above).

3. **Profiles encode community practice.** `joint-disease` profile bakes in Fan 2024 chondrocyte panel + tighter hb% for avascular cartilage. Next user picking that profile inherits that consensus. No agent can derive Fan 2024 by introspection — it has to be told.

4. **Defaults are tested on real cohorts, not improvised.** `max_pct_mt=20` for joint tissue isn't textbook 10% — it was set after losing real chondrocytes at 10% (PI cohort 2024–2026, OARSI working group). A vanilla agent quotes textbook unless told otherwise.

5. **Self-check + auto-fix.** scellrun has telemetry on its own work (QC pass < 60% → suggest threshold relaxation, ≤2 clusters → wider sweep, panel margin all <0.05 → switch panel). A vanilla agent doesn't know what it doesn't know.

6. **Provenance trail in three tiers** — data / inference / literature — visible in every report. Vanilla agent's report has plots without their evidence chain.

7. **Failed retries preserved at `NN_stage.failed-N/`.** When auto-fix kicks in, the original failed state stays on disk for comparison. Agent improvising in a Jupyter notebook overwrites cells.

## What scellrun is NOT

- Not a workflow manager (use nf-core if you need orchestration across stages on a cluster).
- Not a replacement for scanpy (it calls scanpy under the hood with opinionated parameters).
- Not a model (no fine-tuned LLM; the `--ai` calls are deterministic prompts to Claude Haiku).
- Not for the user (that's the LLM agent's job; scellrun is the agent's tool).

## README rewrite — what the new opening should say

Current opening: "scellrun lowers the bar to running a defensible single-cell analysis." (User-facing.)

New opening should reframe to:

> **scellrun is what stops two analysts from getting two different answers on the same data.**
>
> Single-cell analysis has a known reproducibility problem: independent reanalyses commonly find 20% more or fewer clusters from identical raw data, because mt% thresholds, HVG counts, integration methods, and clustering resolutions are usually picked ad-hoc and not recorded. Modern LLM agents can do the work — but two agents handed the same h5ad will pick two different mt% values and produce two different reports, with no audit trail.
>
> scellrun fixes this by encoding the choices: every threshold has a tested default with a one-sentence rationale, every decision is recorded to a `00_decisions.jsonl` file, every report carries a three-tier provenance trail (data / inference / literature), and tissue-specific working practice ships as `profiles/` (cartilage, synovium, brain — contribute yours). When an LLM agent uses scellrun, the resulting report is the same regardless of which agent ran it.

(This is the working draft. Subagent will refine + wire into the full README.)

## Sources

- [Perspectives on rigor and reproducibility in single cell genomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC9122178/)
- [nf-core scrnaseq pipeline](https://nf-co.re/scrnaseq/usage)
- [CellAtria — agentic AI framework for scRNA standardization](https://www.nature.com/articles/s44387-025-00064-0)
- [scnanoseq nf-core pipeline (related)](https://academic.oup.com/bioinformatics/article/41/9/btaf487/8247965)
- [scprocess — recent atlas-scale pipeline](https://www.biorxiv.org/content/10.64898/2026.03.09.710141v1)
