---
name: run-scrna-differential-analysis
description: Audit scRNA-seq differential designs and run sample-level pseudobulk DESeq2 comparisons from raw counts when replication is adequate. Use for condition comparisons, biological-replicate checks, sample-by-condition cell counts, contrast validation, or determining whether a dataset supports formal differential claims.
---

# Run scRNA Differential Analysis

Choose the statistical unit before choosing the test.

## Workflow

1. Define numerator, denominator, cell population, sample column, covariates, and biological question.
2. Audit independent replicate counts and condition-batch or patient-condition confounding.
3. Tabulate cells for every sample-by-population combination and define transparent minimums.
4. Use raw counts for pseudobulk aggregation. Build and validate the sample-level design matrix.
5. Run pseudobulk DESeq2 for formal inference when replication supports it.
6. Apply independent filtering and multiple-testing correction; retain complete result tables.
7. Produce effect-size plots, sample-level expression views, and a limitations summary.

## Guardrails

- Never treat cells as independent biological replicates.
- Stop formal inference when the design is rank-deficient or a condition lacks replication; offer descriptive exploration instead.
- Do not use integrated assay values for pseudobulk counts.
- Keep gene-specific project hypotheses separate from the reusable comparison engine.

Read [references/statistical-design.md](references/statistical-design.md) before executing a comparison.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>`, then rerun with `--execute`. The default is audit-only. Set `analysis.mode` to `pseudobulk` only after reviewing replication and design; the executor refuses formal pseudobulk when either comparison group has fewer than two independent samples.
