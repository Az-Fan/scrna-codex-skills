---
name: 12-scrna-run-pathway-enrichment
description: Run auditable GO, KEGG, Reactome, and Hallmark ORA/GSEA from one or more complete differential-result tables, preserving tested-gene universes, identifier mapping, complete tables, partial database failures, and compact readable plots. Use after differential analysis or with external CSV/TSV/TXT/XLSX results; do not use it to perform differential expression.
---

# Run scRNA Pathway Enrichment

This dedicated entry point reuses the enrichment implementation from differential analysis without rerunning DE.

## Workflow

1. Confirm species, gene identifier type, effect direction, test-statistic columns, whether the input contains all tested genes, and an explicit integer `random_seed`.
2. Create a config from [references/config.example.json](references/config.example.json). Keep `analysis.stage=enrichment_only`.
3. Dry-run `scripts/run.py --config <config>` and inspect table mappings and requested databases.
4. Execute with `--execute` and inspect `task_status.tsv`, mapping tables, complete results, and database-level statuses.

## Output and plotting contract

- Preserve every complete enrichment result table and the tested-gene ORA universe.
- Keep unmapped identifiers and resource versions auditable.
- Record individual databases as completed, empty, skipped, or failed; retain successful databases when another fails.
- Produce a compact overview plus separate database/direction plots using configurable top-N terms, wrapped labels, bounded dimensions, and pagination rather than one unbounded figure.
- Treat summary plots as views of the complete TSV results, never as substitutes for them.

Read [references/enrichment-design.md](references/enrichment-design.md) for resource and interpretation rules. Dry-run first and use `scripts/run_in_tmux.py` for long confirmed runs.
