---
name: 05-scrna-benchmark-integration
description: Run a configurable scRNA-seq batch-correction benchmark from a Seurat object, comparing selected methods and parameter grids against a required uncorrected baseline across user-selected batch variables, biological-conservation metrics, and diagnostic plots. Use when choosing among Harmony, Seurat RPCA, scVI, scANVI, BBKNN, or supplied embeddings; diagnosing batch effects; testing parameter sensitivity; or auditing overcorrection without treating UMAP mixing as sufficient evidence.
---

# Benchmark scRNA Integration

Compare only explicitly selected methods and parameters. Always retain `none` as the uncorrected baseline.

## Workflow

1. Inspect assays, reductions, metadata, cell counts, batch sizes, biological labels, and condition-batch confounding.
2. Ask the user to select methods, method parameters, batch variables, biological labels, metrics, plots, and optional gene programs. Never infer batch fields from sample names.
3. Create a project-local JSON config from [references/config.example.json](references/config.example.json).
4. Run `scripts/run.py --config <config>` to validate and write the execution manifest. Resolve every error before execution.
5. Run `scripts/run.py --config <config> --execute`. Keep cells, feature space, dimensions, neighbor settings, and seeds comparable across scenarios.
6. Review `recommendation_status.json`, method status, and skipped-metric tables before interpreting scores. An unresolved decision is a valid benchmark outcome, not permission to rank by an ad hoc substitute metric.
7. Report batch removal and biological conservation separately. Recommend a method only with explicit tradeoffs and intended downstream use.

If execution may exceed 10 minutes, has uncertain duration, or could outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed `--execute` command with `scripts/run_in_tmux.py`. Keep validation and dry runs in the foreground.

Read [references/benchmark-criteria.md](references/benchmark-criteria.md) for method semantics, metric requirements, scoring, and output interpretation.

## Selection rules

- Treat every method-parameter combination as a named scenario.
- Run Harmony, RPCA, scVI, scANVI, or BBKNN only when selected. Allow precomputed embeddings when a method was run elsewhere.
- Evaluate each configured `metadata.batch_variables` separately; do not collapse distinct technical effects into one synthetic batch.
- Require a trustworthy biological label for label-dependent conservation metrics. Skip unavailable metrics with a recorded reason.
- Generate only configured plots. Keep UMAP, marker, and program-retention plots diagnostic rather than ranking evidence by themselves.
- Accept optional user-supplied gene programs; do not embed tissue-specific marker sets in the generic workflow.

## Guardrails

- Detect perfect or near-perfect batch-condition confounding before correction. Do not reward removal of confounded biology.
- Leave the recommendation unresolved when all requested ranking metrics are unavailable or batch is perfectly confounded with condition. Diagnostic UMAPs and user-added surrogate metrics do not silently replace the configured benchmark contract.
- Never use integrated values as raw counts for differential expression or pseudobulk testing.
- Do not compare a graph-only representation as if it were an ordinary latent embedding; flag BBKNN and similar methods accordingly.
- Do not declare a winner from one aggregate score. Show component metrics and a batch-versus-biology tradeoff view.
- Preserve raw counts, the input object, the uncorrected representation, configuration, seeds, package versions, and method failures.
