---
name: 14-scrna-visualize-cell-composition
description: Create auditable scRNA-seq visualizations of cluster, annotation, cell-state, and hierarchical composition across samples, conditions, batches, and other metadata groups. Use for descriptive composition summaries and batch-effect diagnostics; do not use it as a formal cell-abundance test or to convert clusters into annotations.
---

# Visualize scRNA Cell Composition

Use this skill when the question is “how are clusters, annotations, or states distributed across samples, batches, or conditions?” It supports both unsupervised cluster columns and confirmed annotation columns without conflating their meanings.

## Modes

- `composition_summary`: sample-level counts and proportions, condition summaries, stacked bars, dot plots, and heatmaps.
- `batch_diagnostic`: sample/batch/condition composition, coverage and imbalance tables, and optional UMAP panels. A batch column with one observed level is reported as unavailable, not interpreted.
- `hierarchical_composition`: parent-to-child composition when a parent metadata column is supplied; proportions are recomputed within each parent.

## Workflow

1. Confirm the Seurat RDS/QS object or a sample-by-category count table, sample column, optional condition/batch columns, and one or more configured composition variables.
2. State whether each variable is a `cluster`, `annotation`, `state`, or generic metadata field. Cluster IDs remain cluster IDs in every output.
3. Declare the denominator (`all_input_cells` or a named parent population). Compute proportions within sample, never after pooling cells across samples.
4. Dry-run `scripts/run.py --config <config>`, review resolved fields and output paths, then execute.
5. Inspect the complete count/proportion tables, sample coverage audit, figure status, and batch diagnostic status before interpretation.

## Figures

The default paper-style set uses deterministic label colours and sample-level points:

- composition overview: sample 100% stacked bars plus condition summaries;
- sample dot plot: one point per sample, with group summaries where available;
- sample-by-category heatmap;
- absolute count companion plot;
- optional UMAP panels by sample, batch, condition, and selected composition variable.

The layout is inspired by FigureYa's compact modular biomedical summaries, while retaining the project's `paper_v1` white-background, 300-dpi, paginated contract. Do not use pie charts by default. Do not add significance stars; hand formal abundance results from `13-scrna-test-cell-abundance` only as explicitly labelled annotations.

## Guardrails

- Keep sample-level observations visible; pooled cells are descriptive only and are not biological replication.
- Always display or record the denominator, total cells per sample, and low-coverage category/sample combinations.
- Show counts beside proportions; relative composition is not absolute tissue abundance.
- Preserve zeros and missing combinations in the complete tables.
- Stop or mark unavailable when a requested metadata field is absent, has one level, or is confounded with the requested grouping.
- Do not filter, annotate, recluster, integrate, or modify the input object.
- Use `13-scrna-test-cell-abundance` for formal replicated inference and `05-scrna-benchmark-integration` for quantitative integration benchmarking.

Read [references/input-output-contract.md](references/input-output-contract.md) for schemas, modes, and output interpretation. Read [references/figure-design.md](references/figure-design.md) when changing the default figure family.

## Result organization

Keep primary figures and complete scientific tables at the result root. Store manifests, session information, field resolution, colour maps, and plotting status under `_provenance/`.
