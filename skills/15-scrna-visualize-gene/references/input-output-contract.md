# Input and output contract

## Required input

- `input.object`: Seurat RDS/QS object with an existing normalized expression layer.
- `genes`: non-empty array of `{symbol, label}` records. `label` is optional and defaults to `symbol`.
- `metadata.sample`: biological sample column.
- `expression.assay`: assay used for expression visualization.
- `output_dir` and `project.id`.

`metadata.condition`, `metadata.population`, and `metadata.reduction` are optional in principle, but their dependent figures are skipped when absent. `input.differential_table` must be a complete delimited DE table. `input.pseudobulk_data` may be an RDS list containing `normalized_counts` and `coldata`.

For a combined DE table, configure `differential_selection`, for example `{"population": "EC", "comparison_id": "stz_vs_control"}`. Keys name actual table columns and values select one exact level. The selected table must describe one task and contain one row per gene; ambiguous tasks or duplicated genes stop before result tables are written. Keep all tested genes within the selected task as the volcano background. Scope columns and numerator/denominator are retained in the target-gene summary when supplied.

Pseudobulk `coldata` row names must match the normalized-count sample IDs; metadata is reordered to matrix column order before plotting.

## Primary tables

- `gene_status.tsv`: requested symbols/labels and availability in object, DE table, and pseudobulk matrix.
- `target_gene_cell_expression_summary.tsv`: gene × sample × condition × population cell counts, means, medians, and detected fractions.
- `target_gene_sample_expression.tsv`: plotted sample-level values and their source.
- `target_gene_summary.tsv`: selected-gene DE values plus availability fields.
- `plot_status.tsv`: every fixed figure family, generated or skipped with reason.

## Provenance

`_provenance/` contains `run_manifest.json`, `session_info.txt`, `figure_colors.tsv`, and the shared figure status record. Figures are deterministic for a fixed object, config, and software environment.

`RESULTS.md` is the portable result entry point; `_provenance/result_delivery.json` records the current run and explicitly retained review-stage files. Comparison task directories have their own linked entry points when present.
