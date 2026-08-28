# Interpretation and inference

## Cell-level interpretation

Use cell scores to visualize heterogeneity, co-activity, and localization on an existing embedding. Check whether high scores are driven by a few genes, low coverage, sequencing depth, stress, cell cycle, or cell-type composition. A program score is evidence of coordinated expression, not direct proof of pathway flux, protein activation, or causal signaling.

## Group summaries

Summaries created by `summarize_by` are descriptive means over cells. Include `sample` together with cell type and condition when the biological design contains replicates. Inspect per-sample consistency before interpreting condition patterns.

## Condition inference

Do not use pooled cells as independent replicates. For formal comparisons, use independent samples as the inferential unit, stratify by the intended cell population, and use a paired or blocked design when applicable. Report effect sizes, uncertainty, multiple-testing correction, numbers of samples, and numbers of cells contributing to each sample summary.

Do not call a pathway condition-specific when the result is driven by one sample, a shifting cell-type mixture, or a method-specific scale artifact.

## Cross-method comparisons

Compare methods only when useful to answer an explicit robustness question. Hold the object, cell set, signatures, and identifiers fixed. Report within-signature Spearman/Pearson correlations, group-level concordance, rank stability, missing genes, and method parameters. Do not require raw score equality between algorithms.
