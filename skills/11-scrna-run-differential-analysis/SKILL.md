---
name: 11-scrna-run-differential-analysis
description: Run auditable scRNA-seq differential expression from a Seurat RDS/QS object for one or many populations, using replicated sample-level pseudobulk DESeq2 for formal inference or explicitly exploratory Seurat cell-level tests. Produce complete direction-annotated tables, design audits, sample-level diagnostics, and plots. Use for condition contrasts and batch DEG analysis; use 12-scrna-run-pathway-enrichment when starting from an existing differential table.
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
8. Prefer handing the complete result table to `12-scrna-run-pathway-enrichment`. The legacy `differential_and_enrichment` mode remains available so existing combined calculations and outputs are not lost.

If execution may exceed 10 minutes, has uncertain duration, or could outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed `--execute` command with `scripts/run_in_tmux.py`. Keep validation and dry runs in the foreground.

## Modes

- Use `pseudobulk_deseq2` by default for formal condition comparisons with biological replicates.
- Use `seurat_wilcox`, `seurat_mast`, or `seurat_lr` only for explicit cell-level exploration; report `inference_level=cell_level_exploratory`.
- Use `analysis.stage=differential` for the primary DE-only workflow.
- `analysis.stage=differential_and_enrichment` is retained as a compatibility mode and uses the same established enrichment implementation as skill 12.

## Guardrails

- Never treat cells as independent biological replicates for formal condition claims.
- Never use integrated assay values as pseudobulk counts; select a raw-count RNA-like assay.
- Stop formal inference for rank-deficient designs, samples mapped to multiple conditions, nonconstant sample covariates, or insufficient independent samples.
- Preserve every returned gene. Use thresholds to annotate `Up`, `Down`, `NS`, or `Not_tested`, not to truncate the full table.
- Use tested genes as the ORA universe. Keep identifier mapping and unmapped genes auditable.
- Prefer sample-level normalized expression for formal heatmaps and PCA. Treat cell-level plots as descriptive.
- Do not overwrite the input Seurat object.

Read [references/input-output-contract.md](references/input-output-contract.md) for fields, outputs, task statuses, and compatibility notes.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
