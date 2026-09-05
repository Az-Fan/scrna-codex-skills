# Input and output contract

## Input

Accept a Seurat `.rds` or `.qs` object. Require nonmissing `metadata.sample` and `metadata.condition`. Use `metadata.population` for a cell-type or cluster column; omit it to compare all cells. `metadata.cell_type` remains a backward-compatible alias.

`population.mode=all` runs every observed population after exclusions. `population.mode=selected` runs labels in `population.include`. Expand every population against every item in `comparisons`. A legacy singular `comparison` object remains accepted.

Supported methods are `pseudobulk_deseq2`, `seurat_wilcox`, `seurat_mast`, and `seurat_lr`. Formal pseudobulk requires `DESeq2`; plots require `ggplot2`. Enrichment additionally requires `clusterProfiler`, the matching organism annotation package, and valid gene identifiers.

For `analysis.stage=enrichment_only`, accept one `input.differential_table` or an array of `input.differential_tables`. Each array item requires `path` and may set `sheet`, `id`, `population`, `comparison_id`, `numerator`, and `denominator`. Read CSV by comma, TSV/TXT by tab, and XLS/XLSX with `readxl`.

Auto-detect common columns from DESeq2, Seurat, edgeR, and limma. Prefer explicit `enrichment.table_columns` mappings when names are ambiguous. Require a gene identifier plus either a signed statistic or log fold change. Derive ORA Up/Down only from adjusted P value and log-fold-change thresholds unless an explicit significance column is mapped. Never infer significance from a generic sign-only `direction` column.

If imported tables contain `population`/`celltype`/`cluster` and `comparison_id`/`contrast`, split them into independent enrichment tasks automatically. Write `standardized_input_table.tsv` and `input_column_mapping.tsv` for every task.

## Outputs

Write root-level `design_audit.tsv`, `task_status.tsv`, `all_comparisons.tsv`, `significant_all_comparisons.tsv`, an optional `enrichment_all_comparisons.tsv`, `DEG_count_summary.pdf`. Technical records are `_provenance/session_info.txt` and `_provenance/run_manifest.json`.

For every `population × comparison`, write a directory below `comparisons/` containing:

- `sample_cell_counts.tsv` and, for pseudobulk, `sample_design.tsv`, `effect_size_audit.tsv`, and `pseudobulk_data.rds`.
- `all_genes.tsv`, `significant_genes.tsv`, `upregulated_genes.tsv`, and `downregulated_genes.tsv`.
- `volcano.pdf`, optional `MA_plot.pdf`, and pseudobulk `pseudobulk_PCA.pdf` and `top_DE_heatmap.pdf`.
- Optional `enrichment/` identifier mapping, database-level status audit, full GO-BP/MF/CC, KEGG, Reactome, and Hallmark ORA/GSEA tables, and summary plots.
- `ERROR.txt` when that task cannot run.

DE `task_status.tsv` uses `completed`, `skipped_low_replicates`, `invalid_design`, `missing_dependency`, or `failed`; its `enrichment_status` records `completed`, `partial`, `empty`, `failed`, or `not_requested`. Enrichment-only task status also permits `partial` when at least one requested database succeeds and another fails. One failed population does not discard successful populations.

The complete table includes comparison direction, raw and adjusted P values, effect size, `direction`, `significance`, `tested`, `filter_reason`, method, inference level, and cell/sample counts. `Not_tested` usually reflects independent filtering or an unavailable P value.
