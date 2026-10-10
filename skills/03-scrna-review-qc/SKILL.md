---
name: 03-scrna-review-qc
description: Auto-detect available QC metrics in a Seurat RDS/QS object and generate a concise sample-, cluster-, annotation-, and UMAP-aware QC atlas with threshold and sample-summary tables, plus optional full diagnostics. Use before or after filtering to review QC without modifying or filtering the object.
---

# Review scRNA-seq QC

Use the existing project pixi environment. Never create, install, or update an environment.
The runner requires an existing lock file and installed R interpreter, and uses `pixi run --frozen --no-install`. Read [references/compatibility.md](references/compatibility.md) when diagnosing dependencies.

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

3. Show the resolved input, pixi manifest, environment, output directory, and command. Existing authorization to run or revise the review is sufficient.
4. Execute within the authorized scope:

   `python scripts/run.py --config CONFIG.json --execute`

Read `references/qc-review.md` only when interpreting the threshold table or changing metric aliases.

## Guarantees

- Produce figures and review tables only; never filter cells or write a filtered Seurat object.
- Treat thresholds as candidates requiring biological review.
- Keep approval, decision, and notes fields blank.
- Preserve the input object.
- Keep sample, cluster and UMAP overview sections compact, with optional PNG previews. Preserve established annotation and secondary diagnostics in a supplemental PDF.
- Read [references/compact-figures.md](references/compact-figures.md) for display units, labels, rule overlays and pagination.

## Outputs

The result root contains `qc_atlas.pdf` (sample, cluster and shared-coordinate UMAP sections), `qc_supplement.pdf`, `threshold_review.tsv`, `qc_summary_by_sample.tsv`, and display/sample-composition tables. Large sample or cluster sets paginate without dropping groups. Missing clusters or UMAP produce availability panels. Set `output.figure_format=png` to use standalone PNG figures instead of PDF; `output.preview_png=true` is a compatibility alias for that choice. Technical records are stored under `_provenance/`. Set `output.detail_level` to `full` for availability, quantile, hypothetical retention and plot-status tables under `details/`; individual diagnostic PNGs additionally require `output.preview_png=true`.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.

Result handoff: start with `<output_dir>/RESULTS.md` from the execute runner. It links scientific deliverables and reports current run status and limitations; follow [output-layout.md](references/output-layout.md) for retained results and technical records.

Figure output: select one format per run, default PDF; follow [figure-output.md](references/figure-output.md). Never automatically export both PDF and PNG.
