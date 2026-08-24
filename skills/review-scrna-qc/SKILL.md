---
name: review-scrna-qc
description: Generate a sample-aware scRNA-seq QC figure atlas, metric availability report, candidate threshold review table, and retention sensitivity tables from a Seurat RDS/QS object. Use after QC metrics have been calculated and before filtering; it never filters cells or modifies the input object.
---

# Review scRNA-seq QC

Use the existing project pixi environment. Never create, install, or update an environment.

## Required inputs

- A Seurat `.rds` or `.qs` object containing raw cell metadata and available QC columns.
- The existing pixi project path.
- The metadata column identifying samples.
- An explicit output directory. If it is absent, ask the user where results should be saved before execution.

Condition and batch columns are optional. Missing optional QC metrics are skipped and recorded rather than treated as errors.

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

The output directory contains a multipage `qc_atlas.pdf`, standalone PNG figures, metric availability, sample summaries and quantiles, candidate threshold review, candidate retention summaries, and `run_manifest.json`.
