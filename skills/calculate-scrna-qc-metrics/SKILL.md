---
name: calculate-scrna-qc-metrics
description: Calculate comprehensive per-cell QC metrics for mouse scRNA-seq samples from STARsolo Gene/filtered and Velocyto/filtered matrices, following the sc06 implementation. Use when computing counts, detected genes, mitochondrial, chromosome-Y, nuclear, ambient-RNA, doublet, cell-cycle, and hemoglobin metrics before filtering.
---

# Calculate scRNA QC Metrics

Use the project's existing pixi environment. Never create an environment or install packages.

## Required conversation

Before execution, inspect paths read-only and obtain:

- the existing pixi project directory or `pixi.toml`;
- the STARsolo root containing one directory per sample;
- a GTF matching the STARsolo reference;
- sample IDs and batch IDs;
- the result output directory.

Always ask where results should be saved when the user has not provided an output directory. Do not invent or silently default `output_dir`.

## Execution

1. Create a project-local JSON config from [references/config.example.json](references/config.example.json).
2. Run `scripts/run.py --config <config>` for a dry run.
3. Confirm that the manifest uses the intended existing pixi environment, inputs, samples, and output directory.
4. Run `scripts/run.py --config <config> --execute` only after confirmation.

The executor computes metrics and writes new files; it never filters cells or changes the pixi environment.

## Boundaries

- Support the mouse STARsolo layout implemented and tested in `sc06`.
- Require `Gene/filtered` and `Velocyto/filtered` matrices for every sample.
- Treat the `n_genes` threshold as a diagnostic flag only.
- Stop and report missing packages or files; never install or repair dependencies automatically.
- Do not infer sample batches from sample names.

Read [references/input-output.md](references/input-output.md) when preparing or reviewing a run.
