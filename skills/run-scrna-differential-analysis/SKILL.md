---
name: run-scrna-differential-analysis
description: Design and run scRNA-seq differential expression, cluster marker, pseudobulk, enrichment, and pathway-focused analyses. Use for condition comparisons within cell types, marker discovery, sample-level pseudobulk inference, design-matrix review, enrichment of DEG results, or determining whether a public dataset supports formal differential claims.
---

# Run scRNA Differential Analysis

Choose the statistical unit before choosing the test.

## Workflow

1. Define numerator, denominator, cell population, sample column, covariates, and biological question.
2. Audit independent replicate counts and condition-batch or patient-condition confounding.
3. Tabulate cells for every sample-by-population combination and define transparent minimums.
4. Use raw counts for pseudobulk aggregation. Build and validate the sample-level design matrix.
5. Run pseudobulk DESeq2 or edgeR for formal inference when replication supports it.
6. Use cell-level FindMarkers only as exploratory or marker-oriented evidence and label it accordingly.
7. Apply independent filtering and multiple-testing correction; retain complete result tables.
8. Run enrichment from a declared universe and ranked or thresholded gene list.
9. Produce effect-size plots, sample-level expression views, and a limitations summary.

## Guardrails

- Never treat cells as independent biological replicates.
- Stop formal inference when the design is rank-deficient or a condition lacks replication; offer descriptive exploration instead.
- Do not use integrated assay values for pseudobulk counts.
- Keep gene-specific project hypotheses separate from the reusable comparison engine.

Read [references/statistical-design.md](references/statistical-design.md) before executing a comparison.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>`. Audit replication and design rank before adding `executor.argv`. Use an argv array for a versioned differential-analysis entrypoint; refuse formal pseudobulk execution when the design audit fails.
