---
name: annotate-scrna-cells
description: Cluster and annotate scRNA-seq cells using cluster markers, canonical markers, optional reference-based labels, and explicit human review. Use when assigning broad cell types or subtypes, comparing resolutions, preparing an annotation review table, or writing confirmed annotations back to Seurat or AnnData objects.
---

# Annotate scRNA Cells

Combine evidence; never turn one automated label source into final truth.

## Workflow

1. Confirm the assay, reduction, clustering basis, and metadata contract.
2. Evaluate several biologically plausible resolutions and record cluster stability and size.
3. Compute cluster markers with detection fraction, effect size, and adjusted significance.
4. Plot canonical positive and exclusion markers across clusters and samples.
5. Add reference-based annotation when an appropriate species/tissue reference exists.
6. Generate `annotation-review.tsv` with cluster, candidate label, evidence, conflicts, QC flags, and decision status.
7. Pause for human review when labels are ambiguous or imply deletion.
8. Write confirmed labels into a new metadata column and save annotation provenance.

## Guardrails

- Keep cluster IDs separate from biological labels.
- Detect sample-specific clusters and possible doublets before naming rare types.
- Do not hard-code tissue-, species-, or cell-type-specific markers into general execution logic.
- Preserve broad and fine labels as separate columns.

Read [references/annotation-review.md](references/annotation-review.md) before producing the review table.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>` before clustering. Use the manifest to record inputs and outputs, then rerun with `--execute` to generate clusters, markers, a UMAP, and the review table. Never place confirmed labels in the executor configuration before human review.
