---
name: review-scrna-qc
description: Diagnose and review sample-aware quality control for scRNA-seq count matrices or Seurat/AnnData objects. Use for QC reporting, threshold proposals, filtering sensitivity analysis, doublet or ambient-RNA review, and before annotation, integration, or focused subpopulation analysis.
---

# Review scRNA QC

Separate deterministic diagnostics from project-specific filtering decisions.

## Workflow

1. Verify raw counts and sample identity with `$standardize-scrna-input` when needed.
2. Calculate library size, detected genes, mitochondrial fraction, and available organism/tissue-specific metrics.
3. Summarize distributions per biological sample, condition, and batch before pooling.
4. Flag outliers using both distribution-aware rules and interpretable candidate thresholds.
5. Review doublet and ambient-RNA evidence when tools or results are available.
6. Compare candidate filtering scenarios by sample retention, cell-type composition, and embedding structure.
7. Generate a decision table that distinguishes retain, exclude, and manual-review cells or clusters.
8. Apply filtering only after thresholds or exclusions are confirmed; retain original object and cell IDs.

## Guardrails

- Do not reuse fixed thresholds merely because they worked for another dataset.
- Do not remove a cluster solely for high mitochondrial fraction; inspect markers, sample concentration, counts, and doublet evidence together.
- Report whether exclusions disproportionately affect a condition or batch.
- Treat filtering as a sensitivity analysis when borderline cells could affect conclusions.

Read [references/qc-review.md](references/qc-review.md) for required tables and plots.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json), run `scripts/run.py --config <config>` in the selected compute context, inspect the manifest, then rerun with `--execute`. Use the bundled default R executor unless the environment requires an explicit argv override. Keep threshold approval outside the executor; the default executor never filters cells.
