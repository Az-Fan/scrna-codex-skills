---
name: 08-scrna-analyze-subset
description: Prepare a focused scRNA-seq population for reanalysis by selecting configured labels and exporting raw counts, metadata, cell IDs, sample summaries, provenance, and a derivative Seurat object. Use before reclustering, subtype annotation, integration comparison, or exclusion-sensitivity analysis.
---

# Analyze scRNA Subset

Apply the same auditable reanalysis pattern to any selected population.

## Workflow

1. Define the source object, formal annotation column, inclusion labels, and analysis question.
2. Export cell IDs, counts, metadata, and provenance; avoid copying unnecessary full objects at every stage.
3. Re-run subset-appropriate QC diagnostics without automatically applying whole-dataset thresholds.
4. Establish an uncorrected baseline and assess sample-specific structure.
5. Run targeted integration comparisons only when batch effects warrant them.
6. Compare resolutions, compute markers, and create a subtype review table.
7. Flag contamination, doublets, low-quality states, and sample-restricted clusters using multiple evidence types.
8. Compare exclusion scenarios before removing borderline clusters.
9. Write confirmed subtypes or states to a derivative object; update the parent only when explicitly requested.

## Guardrails

- Do not treat a low-quality state as contamination without marker and sample evidence.
- Do not save repeated full Seurat objects when counts, metadata, IDs, and embeddings suffice.
- Keep subtype identity, activation state, QC status, and exclusion decision in separate fields.

Read [references/subset-review.md](references/subset-review.md) for the expected checkpoints.

## Execution

Create a JSON config from [references/config.example.json](references/config.example.json) and run `scripts/run.py --config <config>`. Confirm the inclusion labels from the source object's metadata, then rerun with `--execute`. The default executor exports counts, metadata, a summary, cell IDs, and a derivative object; keep parent-object writeback as a separate explicit action.
