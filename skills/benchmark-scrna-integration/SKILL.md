---
name: benchmark-scrna-integration
description: Build an uncorrected scRNA-seq integration baseline and optionally compare it with Harmony while preserving biological structure and sample provenance. Use when diagnosing batch effects, establishing a required no-correction baseline, or evaluating whether Harmony improves batch mixing without obvious biological loss.
---

# Benchmark scRNA Integration

Treat no correction as a required baseline and method selection as evidence-based.

## Workflow

1. Inspect the experimental design and identify possible condition-batch confounding.
2. Build an uncorrected baseline with fixed preprocessing choices.
3. Run only methods supported by the data size, environment, and annotation availability.
4. Keep feature selection, neighbor count, dimensions, seeds, and evaluation labels comparable.
5. Evaluate sample mixing, batch predictability, cell-type separation, local structure, cluster stability, and condition preservation.
6. Inspect embeddings by sample, batch, condition, cell type, and QC metrics.
7. Compare sensitive Harmony parameters when the conclusion depends on them.
8. Produce a benchmark table and a recommendation with explicit tradeoffs.

## Guardrails

- Do not select a method from UMAP appearance alone.
- Do not reward mixing when batch is confounded with real biology.
- Preserve the uncorrected representation and raw counts.

Read [references/benchmark-criteria.md](references/benchmark-criteria.md) for the evaluation dimensions.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>` first. Require a valid batch field and an uncorrected baseline, then rerun with `--execute`. The default executor always produces the baseline and optionally runs configured supported methods; preserve identical inputs across methods.
