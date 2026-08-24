---
name: calculate-scrna-qc-metrics
description: Calculate available per-cell QC metrics from STARsolo matrices or a Seurat RDS/QS object, following the sc06 implementation. Use before filtering to add counts, detected genes, mitochondrial, chromosome-Y, nuclear, ambient-RNA, doublet, cell-cycle, and hemoglobin metrics while recording metrics that cannot be computed.
---

# Calculate scRNA QC Metrics

Use the project's existing pixi environment. Never create an environment or install packages.

## Required conversation

Before execution, inspect paths read-only and obtain:

- the existing pixi project directory or `pixi.toml`;
- either a STARsolo root or a Seurat RDS/QS object;
- species and relevant metadata columns;
- optional GTF, Velocyto matrices, sample IDs, and batch IDs when available;
- the result output directory.

Always ask where results should be saved when the user has not provided an output directory. Do not invent or silently default `output_dir`.

## Execution

1. Create a project-local JSON config from [references/config.example.json](references/config.example.json).
2. Run `scripts/run.py --config <config>` for a dry run.
3. Confirm that the manifest uses the intended existing pixi environment, inputs, samples, and output directory.
4. Run `scripts/run.py --config <config> --execute` only after confirmation.

The executor computes metrics and writes new files; it never filters cells or changes the pixi environment.

## Boundaries

- Support STARsolo matrices and Seurat RDS/QS objects with raw counts.
- Use optional GTF and Velocyto inputs when available; record skipped metrics when they are absent.
- Treat the `n_genes` threshold as a diagnostic flag only.
- Stop only when the primary object or count matrix is unreadable. Record optional missing inputs or packages as skipped metrics.
- Never install or repair dependencies automatically.
- Do not infer sample batches from sample names.

Read [references/input-output.md](references/input-output.md) when preparing or reviewing a run.
