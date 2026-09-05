---
name: 09-scrna-export-subset
description: Select configured labels from a Seurat RDS/QS object and export an auditable derivative subset with raw counts, metadata, feature and barcode files, sample summaries, and provenance. Use to hand a focused population to QC review, reclustering, annotation, integration benchmarking, or differential analysis; do not use this exporter as if it performed those downstream analyses.
---

# Export scRNA Subset

Create one auditable handoff without modifying the parent object.

## Workflow

1. Confirm the source object, formal annotation column, exact inclusion labels, sample field, and output directory.
2. Create a config from [references/config.example.json](references/config.example.json).
3. Dry-run `scripts/run.py --config <config>` and verify the resolved labels and paths.
4. Execute with `--execute`.
5. Pass the derivative object to the appropriate downstream skill; subset-specific QC thresholds, integration, clustering, marker detection, and annotation remain separate decisions.

If export may exceed 10 minutes or outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed command with `scripts/run_in_tmux.py`.

## Outputs

Preserve the existing complete output contract: `subset_object.qs`, `subset_counts.mtx`, `subset_metadata.tsv`, `subset_summary.tsv`, `features.tsv`, `barcodes.tsv`, and `_provenance/run_manifest.json`.

## Guardrails

- Stop when no cell matches or required labels are missing.
- Preserve raw counts, complete selected-cell metadata, cell order, features, and barcodes.
- Do not write labels back to the parent object.
- Do not claim that export alone performed reanalysis.
- Avoid additional full-object copies outside the one derivative handoff object.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
