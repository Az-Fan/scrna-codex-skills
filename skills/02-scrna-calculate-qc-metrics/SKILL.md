---
name: 02-scrna-calculate-qc-metrics
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
- whether to run DecontX now, and an optional broad-cell-population metadata column when reliable labels already exist;
- the number of independent samples to process concurrently;
- the result output directory.

Always ask where results should be saved when the user has not provided an output directory. Do not invent or silently default `output_dir`.

## Execution

1. Create a project-local JSON config from [references/config.example.json](references/config.example.json).
2. Run `scripts/run.py --config <config>` for a dry run.
3. Confirm that the manifest uses the intended existing pixi environment, inputs, samples, and output directory.
4. Run `scripts/run.py --config <config> --execute` only after confirmation.

If execution may exceed 10 minutes, has uncertain duration, or could outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed `--execute` command with `scripts/run_in_tmux.py`. Keep validation and dry runs in the foreground.

The executor computes metrics and writes new files; it never filters cells or changes the pixi environment.
Keep user-facing data and figures at the result root. Put executor-only status and provenance records under `_provenance/`; inspect them when auditing a run, but do not present them as primary results.

For Seurat v5 input, let the executor join multiple raw-count layers in memory and verify that the resulting matrix covers every object cell. Do not require callers to modify or resave their source object first.

Set `ambient_rna.method` to `skip` for a fast run that records DecontX as deliberately skipped. Set `ambient_rna.cluster_column` only for complete, broad population labels; supplying it bypasses DecontX's cell-population estimation, although celda still generates a UMAP for its result object. Use `parallel.workers` to process independent samples concurrently, bounded by available memory and the number of samples.

## Boundaries

- Support STARsolo matrices and Seurat RDS/QS objects with raw counts.
- Use optional GTF and Velocyto inputs when available; record skipped metrics when they are absent.
- Treat the `n_genes` threshold as a diagnostic flag only.
- Stop only when the primary object or count matrix is unreadable. Record optional missing inputs or packages as skipped metrics.
- Never install or repair dependencies automatically.
- Do not infer sample batches from sample names.

Read [references/input-output.md](references/input-output.md) when preparing or reviewing a run.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
