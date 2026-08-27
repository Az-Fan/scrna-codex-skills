---
name: scrna-find-cluster-markers
description: Identify differential marker genes across existing Seurat clusters with FindAllMarkers and produce complete, ranked, annotation-ready marker tables, cluster summaries, and a marker dot plot. Use after clustering and before cell-type annotation when comparing clusters, selecting positive or bidirectional markers, or auditing whether every cluster has usable markers; do not use for condition-level differential expression across biological samples.
---

# Find scRNA Cluster Markers

Treat cluster markers as annotation evidence, not final cell-type labels.

## Workflow

1. Confirm that the input is a clustered Seurat object and select normalized expression; permit in-memory log-normalization only when the configured assay lacks a data layer.
2. Choose the assay and either a metadata cluster column or the active identities.
3. Check that at least two non-missing groups exist and record their cell counts.
4. Run `Seurat::FindAllMarkers` with explicit test, direction, detection, fold-change, and adjusted-significance settings.
5. Preserve the complete result table. Rank markers within each cluster and extract a configurable top-marker table.
6. Review marker yield per cluster and the dot plot before using the results for annotation.
7. Pass marker evidence to an annotation workflow without writing biological labels into the object.

## Guardrails

- Use cluster-level marker testing for discovery and annotation support only. Use a replicated pseudobulk workflow for condition-level biological claims.
- Do not run on integrated assay values unless the selected method explicitly supports the intended interpretation; prefer normalized RNA or SCT expression.
- Keep negative markers when they are requested; absence can be informative for distinguishing related cell types.
- Do not silently omit clusters with no passing markers. Retain them in the summary with a zero marker count.
- Treat tiny clusters, sample-specific clusters, and stress, cell-cycle, ambient-RNA, or doublet signatures as review flags rather than automatic cell types.

Read [references/marker-interpretation.md](references/marker-interpretation.md) before interpreting or handing markers to annotation.

## Execution

Copy [references/config.example.json](references/config.example.json), set the input and grouping fields, and run `scripts/run.py --config <config>` to validate and write a plan manifest. Review the manifest, then rerun with `--execute`. The executor reads but does not save changes to the input object and writes full markers, ranked top markers, a per-cluster summary, a dot plot, and a run manifest.
