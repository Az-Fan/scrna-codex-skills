---
name: 17-scrna-infer-grn
description: Infer a core scRNA-seq gene regulatory network with GRNBoost2, cisTarget motif pruning and single-cell AUCell regulon activity from Seurat RDS/QS or raw matrix bundles. Support reviewed metacell memberships or explicit single-cell inference, resource coverage audits and staged input review. Use for constructing TF-target candidates and regulons; downstream CSI/RSS, differential regulons, pathway-driver TFs and target-gene screening are separate analyses.
---

# Infer GRN and score regulons

Produce a traceable TF-target network, motif-supported regulons and aligned single-cell activity. GRNBoost2 importance is a predictive association weight; motif support and AUC do not establish causal regulation or TF protein activity.

## Workflow

1. Confirm the QC-filtered input, population scope, sample column, species, genome and identifiers. Require uncorrected raw RNA counts. Do not infer scope from sc06's mouse/EC example or silently change populations.
2. Read [references/grn-contract.md](references/grn-contract.md), then adapt [references/config.example.json](references/config.example.json). Locate existing TF lists, motif annotations and cisTarget ranking databases; explicitly match their species/genome/identifier declarations. Do not download or install resources during analysis.
3. Select `inference.mode=metacell` with a reviewed `metacell.column`, or explicitly select `single_cell`. Metacell memberships are split by biological sample and configured cell type, then raw counts are summed. Use skill 06 to generate candidate high-resolution clustering if memberships do not yet exist; resolution review and clustering are not duplicated here.
4. Run `scripts/check_dependencies.py --config <config>` and `scripts/run.py --config <config>` to review the resolved runtimes. Default `workflow.action=prepare` with `--execute` exports input, memberships, size/purity audits and resource coverage, then waits for input review.
5. After review, set `action=infer` and document `workflow.review_reason` in the same output/config. Prepared input, resources, parameters and package versions must match. If the user has already chosen the full input/membership/parameters, `action=run` performs all core stages directly. Read [references/long-running-execution.md](references/long-running-execution.md) for the bundled tmux runner.
6. Run GRNBoost2 → cisTarget motif pruning → regulon export → AUCell on the original single cells. Inspect complete `task_status.tsv`, resource/regulon coverage and the sample summary. Preserve all regulon names, target weights and motif context. One declared minimum regulon size applies to both exported and scored regulons.
7. Deliver `RESULTS.md` and `handoff.json`. Stop after the core network/activity handoff. Ask which downstream question to address when the user is ready; do not automatically run CSI, network community detection, RSS, differential activity, pathway/TF correlations, variance decomposition, hub TFs or target screening.

## Constraints

- Preserve source objects and cell IDs. Matrix orientation is declared; duplicate gene IDs, mismatched metadata, negative/non-integer counts and zero-count cells stop execution.
- Keep inference-unit aggregation separate from activity: metacells infer the network; AUCell scores every original cell with the full original gene universe.
- Record low coverage, constant genes and regulons excluded by size. Empty motif-supported results are a reported failure, not permission to relax thresholds automatically.
- Sample summaries are descriptive. Cells are not biological replicates; no condition-level significance or causal TF claims are produced.
- Prepared/native artifacts stay under `_provenance/`; source data are unchanged. Existing inference results are preserved and parameter experiments use new output directories.

Figure output: select one format per run, default PDF; follow [figure-output.md](references/figure-output.md). Never automatically export both PDF and PNG.
