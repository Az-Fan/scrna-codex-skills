---
name: 10-scrna-score-programs
description: Score curated gene signatures and pathway activity in a Seurat scRNA-seq object with VISION, AUCell, UCell, AddModuleScore, or PROGENy; support inline, GMT, MSigDB, and scMetabolism gene-set resources; audit gene coverage, cache reproducible score matrices, attach namespaced assays, and summarize scores by sample, cell type, or condition. Use for Hallmark, metabolism, signaling, custom signatures, or comparable cell-level program scoring. Use the documented handoff rather than this executor for de novo cNMF discovery.
---

# Score scRNA Programs

Treat pathway analysis as explicit gene-level scoring with versioned resources and auditable coverage.

## Workflow

1. Confirm the Seurat object, species, assay, layer, signature source, scoring method, and requested grouping columns.
2. Read [references/scoring-methods.md](references/scoring-methods.md) when choosing a method or expression layer.
3. Copy [references/config.example.json](references/config.example.json), declare every requested task explicitly, and avoid adding unrequested methods for comparison.
4. Run `scripts/run.py --config <config>` to validate the contract and write a plan manifest. Review it, then rerun with `--execute`.
5. Inspect `signature_coverage.tsv` before interpreting scores. Resolve species or identifier problems when many signatures are skipped or have low coverage.
6. Inspect `_provenance/task_manifest.json` for methods, resources, cache keys, output assays, and cache hits.
7. Use cell-level scores for visualization. Read [references/interpretation-and-inference.md](references/interpretation-and-inference.md) before testing conditions.

If execution may exceed 10 minutes, has uncertain duration, or could outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed `--execute` command with `scripts/run_in_tmux.py`. Keep validation and dry runs in the foreground.

## Guardrails

- Keep the input object unchanged and write a derivative object.
- Use raw counts for VISION, AUCell, and usually UCell; use normalized expression for AddModuleScore and PROGENy.
- Never score integrated or batch-corrected values unless the user explicitly requests and justifies them.
- Require `species` and preserve resource file hashes and package versions.
- Skip signatures with too few matched genes by default; never hide them from the coverage table.
- Do not compare absolute score scales across different methods.
- Do not infer condition effects from pooled cells. Aggregate or model by independent biological sample.
- Treat cNMF as de novo program discovery with rank and stability review; do not present it as another curated signature scorer.

## Outputs

Write a derivative Seurat object, per-task matrices, coverage and grouped-summary tables, grouped heatmaps, task and run manifests, and session information. Attach each score matrix as a namespaced Seurat assay with programs as features and cells as columns. Use `assay_feature_mapping.tsv` because Seurat assay feature names replace underscores with hyphens while exported matrices preserve original signature names.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
