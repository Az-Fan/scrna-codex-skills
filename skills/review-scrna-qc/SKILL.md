---
name: review-scrna-qc
description: Diagnose and review sample-aware quality control for Seurat RDS/QS objects. Use for per-sample QC reporting, threshold proposals, and before annotation, integration, or focused subpopulation analysis; the default executor reports diagnostics and never filters cells.
---

# Review scRNA QC

Separate deterministic diagnostics from project-specific filtering decisions.

## Workflow

1. Verify raw counts and sample identity with `$standardize-scrna-input` when needed.
2. Calculate library size, detected genes, mitochondrial fraction, and available organism/tissue-specific metrics.
3. Summarize distributions per biological sample, condition, and batch before pooling.
4. Flag outliers using both distribution-aware rules and interpretable candidate thresholds.
5. Compare candidate filtering scenarios by sample retention and condition or batch balance.
6. Generate a decision table that distinguishes retain, exclude, and manual-review cells or clusters.
7. Apply filtering only after thresholds or exclusions are confirmed; retain original object and cell IDs.

## Guardrails

- Do not reuse fixed thresholds merely because they worked for another dataset.
- Do not remove a cluster solely for high mitochondrial fraction; inspect markers, sample concentration, counts, and doublet evidence together.
- Report whether exclusions disproportionately affect a condition or batch.
- Treat filtering as a sensitivity analysis when borderline cells could affect conclusions.

Read [references/qc-review.md](references/qc-review.md) for required tables and plots.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json), run `scripts/run.py --config <config>` in the selected compute context, inspect the manifest, then rerun with `--execute`. Use the bundled default R executor unless the environment requires an explicit argv override. Keep threshold approval outside the executor; the default executor never filters cells.
