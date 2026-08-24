---
name: review-scrna-qc
description: Auto-detect available QC metrics in a Seurat RDS/QS object and generate a comprehensive sample-, cluster-, annotation-, and UMAP-aware QC figure atlas plus threshold and retention audit tables. Use before or after filtering to review QC without modifying or filtering the object.
---

# Review scRNA-seq QC

Use the existing project pixi environment. Never create, install, or update an environment.

## Required inputs

- A Seurat `.rds` or `.qs` object containing raw cell metadata and available QC columns.
- The existing pixi project path.
- The metadata column identifying samples.
- An explicit output directory. If it is absent, ask the user where results should be saved before execution.

Condition, batch, cluster, and annotation columns are optional. Cluster and annotation columns are auto-detected when not configured. Missing metrics, group columns, or UMAP coordinates are skipped and recorded rather than treated as errors.

## Run

1. Copy `references/config.example.json` and adapt paths/column names.
2. Dry-run first:

   `python scripts/run.py --config CONFIG.json`

3. Show the resolved input, pixi manifest, environment, output directory, and command.
4. After confirmation, execute:

   `python scripts/run.py --config CONFIG.json --execute`

Read `references/qc-review.md` only when interpreting the threshold table or changing metric aliases.

## Guarantees

- Produce figures and review tables only; never filter cells or write a filtered Seurat object.
- Treat thresholds as candidates requiring biological review.
- Keep approval, decision, and notes fields blank.
- Preserve the input object.
- Record unavailable metrics and skipped plots.

## Outputs

The output directory contains a multipage `qc_atlas.pdf`; sample, cluster, annotation, scatter, threshold, retention, and per-metric UMAP PNGs when supported; `plot_status.tsv`; metric availability; sample summaries and quantiles; candidate threshold and retention tables; and `run_manifest.json`.
