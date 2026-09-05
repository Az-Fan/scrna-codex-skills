---
name: 03-scrna-review-qc
description: Auto-detect available QC metrics in a Seurat RDS/QS object and generate a concise sample-, cluster-, annotation-, and UMAP-aware QC atlas with threshold and sample-summary tables, plus optional full diagnostics. Use before or after filtering to review QC without modifying or filtering the object.
---

# Review scRNA-seq QC

Use the existing project pixi environment. Never create, install, or update an environment.

## Required inputs

- A Seurat `.rds` or `.qs` object containing raw cell metadata and available QC columns.
- The existing pixi project path.
- The metadata column identifying samples.
- An explicit output directory. If it is absent, ask the user where results should be saved before execution.

Condition, batch, cluster, and annotation columns are optional. Cluster and annotation columns are auto-detected when not configured. Missing metrics, group columns, or UMAP coordinates are skipped rather than treated as errors.

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
- Keep the default result directory concise; generate supplemental audit tables and individual PNGs only when explicitly requested.

## Outputs

The result root contains only `qc_atlas.pdf`, `threshold_review.tsv`, and `qc_summary_by_sample.tsv`. Technical records are stored under `_provenance/`. Set `output.detail_level` to `full` only when machine-readable availability, quantile, retention, plot-status tables, and individual PNG files are needed; place those supplements under `details/`.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
