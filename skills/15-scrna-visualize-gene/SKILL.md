---
name: 15-scrna-visualize-gene
description: Create a fixed, auditable paper-style visualization set for one or more target genes in a Seurat scRNA-seq object, with optional differential-expression and pseudobulk context. Use for FeaturePlot, grouped expression, sample-aware expression, heatmaps, and highlighting specified genes in complete DE results; do not use it to run differential analysis or claim significance from cells as replicates.
---

# Visualize scRNA Target Genes

Use this skill when the question is “where and in which cells is this gene expressed?” or “how does a target gene relate to an existing differential result?” It fixes the visual language so repeated requests remain comparable.

## Fixed figure set

1. `target_gene_featureplots`: UMAP expression with a light-grey to deep-red scale, one independent scale per gene.
2. `target_gene_dotplot`: detected fraction and scaled mean expression across the configured population field.
3. `target_gene_violinplot`: cell-level expression by condition; descriptive only.
4. `target_gene_sample_expression`: one point per biological sample, using pseudobulk normalized counts when supplied and otherwise sample mean normalized expression.
5. `target_gene_sample_heatmap`: row-scaled sample-level expression when at least two genes and two samples are available.
6. `target_gene_de_effect`: target-gene log2 fold changes from an existing complete DE table.
7. `target_gene_volcano_highlight`: the complete DE background with requested genes highlighted.

The default is the repository's `paper_v1` white-background, 300-dpi style. The compact panel logic is informed by FigureYa UMAP, single-cell violin, and marker heatmap examples, but the statistical unit and provenance remain explicit.

## Workflow

1. Confirm the Seurat RDS/QS object, target gene symbols, normalized assay, reduction, biological sample field, and optional condition/population fields.
2. Optionally supply a complete DE table and pseudobulk data from `11-scrna-run-differential-analysis`.
3. Dry-run `scripts/run.py --config <config>`, review paths and fields, then execute.
4. Inspect `gene_status.tsv` and `plot_status.tsv`; missing genes and inapplicable figures are skipped with reasons.
5. Interpret cell-level panels descriptively. Use sample-level observations and the supplied formal DE result for condition claims.

## Guardrails

- Never normalize, impute, integrate, test differential expression, or modify the input object.
- Require an existing normalized expression layer for cell-level plots.
- Do not treat cells as biological replicates or add significance stars to violin plots.
- Preserve the complete supplied DE result as the volcano background; do not infer significance from the selected-gene subset.
- Record requested aliases separately from feature symbols. An alias changes labels, not feature lookup.
- Skip rather than fabricate a plot when its required reduction, grouping levels, samples, genes, or DE columns are unavailable.
- Keep per-sample values visible and state whether they came from pseudobulk normalized counts or descriptive sample means.

Read [references/input-output-contract.md](references/input-output-contract.md) for field and output schemas. Read [references/figure-design.md](references/figure-design.md) before changing the fixed figure family.

