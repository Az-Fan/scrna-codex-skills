---
name: benchmark-scrna-integration
description: Compare uncorrected and batch-integrated scRNA-seq representations, including Harmony, Seurat RPCA, scVI, scANVI, and BBKNN when available. Use when diagnosing batch effects, selecting an integration method or parameters, checking overcorrection, or validating that biological structure and sample mixing are both preserved.
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
7. Compare sensitive parameters such as Harmony theta or RPCA anchors when the conclusion depends on them.
8. Produce a benchmark table and a recommendation with explicit tradeoffs.

## Guardrails

- Do not select a method from UMAP appearance alone.
- Do not reward mixing when batch is confounded with real biology.
- Do not use scANVI labels derived from the evaluation target without documenting circularity.
- Preserve the uncorrected representation and raw counts.

Read [references/benchmark-criteria.md](references/benchmark-criteria.md) for the evaluation dimensions.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>` first. Require a valid batch field and an uncorrected baseline. Use `executor.argv` only for a versioned benchmark entrypoint in the selected environment; preserve identical inputs across methods.
