---
name: scrna-run-differential-analysis
description: Run auditable scRNA-seq differential-expression workflows from a Seurat RDS/QS object for one or many cell types, clusters, or whole-object populations, using replicated sample-level pseudobulk DESeq2 for formal condition inference or explicitly exploratory Seurat cell-level tests. Produce complete direction-annotated tables, plots, design audits, and comprehensive GO, KEGG, Reactome, and Hallmark ORA/GSEA; also start directly from external CSV, TSV, TXT, or Excel differential tables without a Seurat object. Use for condition contrasts, batch DEG analysis, volcano or pseudobulk heatmap generation, pathway enrichment, or rerunning enrichment from existing differential results.
---

# Run scRNA Differential Analysis

Choose the statistical unit before choosing the test. Treat positive log2 fold change as numerator greater than denominator.

## Workflow

1. Inspect the Seurat metadata and confirm sample, condition, population, assay, species, gene identifier type, comparison direction, and covariates.
2. Copy [references/config.example.json](references/config.example.json) and define one comparison or a `comparisons` array.
3. Set `population.mode` to `all` for batch analysis or `selected` with explicit labels. Omit `metadata.population` only for a whole-object comparison.
4. Run `scripts/run.py --config <config>` without `--execute`. Review the manifest and validation messages.
5. Review the sample-by-population audit, replication, pairing, confounding, and model formula using [references/statistical-design.md](references/statistical-design.md).
6. Rerun with `--execute`. Let each population-by-comparison task complete, skip, or fail independently.
7. Inspect `task_status.tsv`, full result tables, sample-level plots, and limitations before interpreting genes.
8. Enable enrichment only after confirming species and gene identifiers. By default run the complete GO-BP/MF/CC, KEGG, Reactome, and Hallmark suite described in [references/enrichment-design.md](references/enrichment-design.md).

## Modes

- Use `pseudobulk_deseq2` by default for formal condition comparisons with biological replicates.
- Use `seurat_wilcox`, `seurat_mast`, or `seurat_lr` only for explicit cell-level exploration; report `inference_level=cell_level_exploratory`.
- Use `analysis.stage=differential` for DE only or `differential_and_enrichment` for both stages.
- Use `analysis.stage=enrichment_only` with `input.differential_table` or `input.differential_tables` to start directly from CSV/TSV/TXT/XLSX results; no Seurat object or sample metadata is required. Copy [references/config.enrichment-only.example.json](references/config.enrichment-only.example.json).

## Guardrails

- Never treat cells as independent biological replicates for formal condition claims.
- Never use integrated assay values as pseudobulk counts; select a raw-count RNA-like assay.
- Stop formal inference for rank-deficient designs, samples mapped to multiple conditions, nonconstant sample covariates, or insufficient independent samples.
- Preserve every returned gene. Use thresholds to annotate `Up`, `Down`, `NS`, or `Not_tested`, not to truncate the full table.
- Use tested genes as the ORA universe. Keep identifier mapping and unmapped genes auditable.
- Prefer sample-level normalized expression for formal heatmaps and PCA. Treat cell-level plots as descriptive.
- Do not overwrite the input Seurat object.

Read [references/input-output-contract.md](references/input-output-contract.md) for fields, outputs, task statuses, and compatibility notes.
