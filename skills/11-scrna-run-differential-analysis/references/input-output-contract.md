# Input and output contract

## Input

Accept a Seurat `.rds` or `.qs` object. Require nonmissing `metadata.sample` and `metadata.condition`. Use `metadata.population` for a cell-type or cluster column; omit it to compare all cells. `metadata.cell_type` remains a backward-compatible alias.

`population.mode=all` runs every observed population after exclusions. `population.mode=selected` runs labels in `population.include`. Expand every population against every item in `comparisons`. A legacy singular `comparison` object remains accepted.

Supported methods are `pseudobulk_deseq2`, `seurat_wilcox`, `seurat_mast`, and `seurat_lr`. Formal pseudobulk requires `analysis.counts_source.kind` (`raw_umi` or `raw_read`), an uncorrected assay with finite nonnegative integer counts, and `DESeq2`; plots require `ggplot2`. Enrichment additionally requires `clusterProfiler`, the matching organism annotation package, and valid gene identifiers.

For `analysis.stage=enrichment_only`, accept one `input.differential_table` or an array of `input.differential_tables`. Each array item requires `path` and may set `sheet`, `id`, `population`, `comparison_id`, `numerator`, and `denominator`. Read CSV by comma, TSV/TXT by tab, and XLS/XLSX with `readxl`.

Auto-detect common columns from DESeq2, Seurat, edgeR, and limma. Prefer explicit `enrichment.table_columns` mappings when names are ambiguous. Require a gene identifier plus either a signed statistic or log fold change. Derive ORA Up/Down only from adjusted P value and log-fold-change thresholds unless an explicit significance column is mapped. Never infer significance from a generic sign-only `direction` column.

If imported tables contain `population`/`celltype`/`cluster` and `comparison_id`/`contrast`, split them into independent enrichment tasks automatically. Write `standardized_input_table.tsv` and `input_column_mapping.tsv` for every task.

Optional `analysis.counts_source.reference_object` points to an original Seurat RDS/QS; `reference_assay` defaults to RNA. Every selected gene/cell count must match that reference. `counts_source_audit.tsv` records the declaration and reference SHA256. Neither a declaration nor a user-supplied reference independently authenticates original raw-count provenance. Seurat v4 `Assay` and v5 `Assay5` are read by their actual object class.

## Outputs

Write root-level `design_audit.tsv`, `task_status.tsv`, `all_comparisons.tsv`, `significant_all_comparisons.tsv`, an optional `enrichment_all_comparisons.tsv`, `DEG_count_summary.pdf`. Technical records are `_provenance/session_info.txt` and `_provenance/run_manifest.json`.

Figures default to PDF; `output.figure_format=pdf|png` selects exactly one format per run.

For every `population × comparison`, write a directory below `comparisons/` containing:

- `sample_cell_counts.tsv` and, for pseudobulk, `sample_design.tsv`, `effect_size_audit.tsv`, `replication_audit.tsv`, `pseudobulk_transform_audit.tsv`, `gene_filter_audit.tsv`, `deseq2_results_audit.json`, and `pseudobulk_data.rds`.
- `all_genes.tsv`, `significant_genes.tsv`, `upregulated_genes.tsv`, and `downregulated_genes.tsv`.
- `volcano.pdf`, optional `MA_plot.pdf`, and pseudobulk `pseudobulk_PCA.pdf` and `top_DE_heatmap.pdf`.
- Optional `enrichment/` identifier mapping, database-level status audit, full GO-BP/MF/CC, KEGG, Reactome, and Hallmark ORA/GSEA tables, and summary plots.
- `ERROR.txt` when that task cannot run.

DE `task_status.tsv` uses `completed`, `skipped_low_replicates`, `invalid_design`, `missing_dependency`, or `failed`; its `enrichment_status` records `completed`, `partial`, `empty`, `failed`, or `not_requested`. Enrichment-only task status also permits `partial` when at least one requested database succeeds and another fails. One failed population does not discard successful populations.

The complete table retains `tested` as the compatibility flag for an available P value, and adds `wald_tested` (NA for non-Wald methods), `multiple_testing_eligible` (finite adjusted P value), and `test_status`. A tested gene with unavailable adjusted P is `Filtered`, not `Not_tested`. `filter_reason` distinguishes DESeq2 `independent_filtering`, `zero_counts`, `cooks_outlier`, and `test_unavailable`; imported tables with unknown causes use `adjusted_p_unavailable`. `Effect_unavailable` records a tested/eligible gene without a usable effect size. DESeq2 independent filtering uses `alpha=analysis.padj_threshold`; `deseq2_results_audit.json` records its threshold and state counts. Cook filtering is identified by comparison with results from the same fitted model using `cooksCutoff=false`; inferential P values retain the primary DESeq2 settings. Genes excluded before model fitting remain in `gene_filter_audit.tsv`.

Table import preserves supplied `tested`, `wald_tested`, `multiple_testing_eligible`, and `filter_reason` fields. Boolean fields accept TRUE/FALSE or 1/0; contradictions such as eligible rows without finite padj are rejected. ORA defaults to `enrichment.universe_mode=multiple_testing_eligible`; `tested` is an explicit alternative background. GSEA retains every finite signed statistic regardless of adjusted-P eligibility. Rank-only tables can run GSEA; they have no default ORA background without adjusted P values.

PCA and the diagnostic heatmap default to actual DESeq2 VST (`varianceStabilizingTransformation(dds, blind=false)`), including small gene sets that cannot use the fast `vst()` default subsample. Set `analysis.pca_transform=log2_normalized` only to request the explicitly named `log2(normalized_counts+1)` alternative. Record the transform in the TSV audit, PCA subtitle, and `pseudobulk_data.rds`; neither transform replaces the raw counts used for DE. No silent transform fallback is performed.

Formal results from fewer than three retained samples in either group carry `low_replication_warning=true` and `inference_qualification=exploratory_low_confidence`. Two per group remains a technical minimum; `replication_audit.tsv` records actual retained replication and task status also carries the warning. This qualification describes evidence strength, while `inference_level` retains the statistical unit/method.

`RESULTS.md` is the portable result entry point; `_provenance/result_delivery.json` records the current run and explicitly retained review-stage files. Comparison task directories have their own linked entry points when present.
