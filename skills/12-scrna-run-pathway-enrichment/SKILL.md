---
name: 12-scrna-run-pathway-enrichment
description: Run auditable GO, KEGG, Reactome, and Hallmark ORA/GSEA from one or more complete differential-result tables, preserving test-eligible gene universes, identifier mapping, complete tables, partial database failures, and compact readable plots. Use after differential analysis, with external differential tables, or to replot completed enrichment tables as compact atlases, evidence/core-gene matrices and optional GO semantic modules without rerunning enrichment.
---

# Run scRNA Pathway Enrichment

This dedicated entry point reuses the enrichment implementation from differential analysis without rerunning DE.

For multiple `input.differential_tables`, give each table a distinct `id` when stable names matter; omitted IDs become `table_1`, `table_2`, etc. Multi-table task paths include this table ID before population/comparison, and sanitized collisions are rejected. Combined enrichment rows retain `task_id` and `input_table`. Existing numerator/denominator columns are preserved unless the table specification explicitly overrides them.

## Workflow

1. Confirm species, gene identifier type, effect direction, test-statistic columns, whether the input contains the complete test results and eligibility flags, and an explicit integer `random_seed`.
2. Create a config from [references/config.example.json](references/config.example.json). Keep `analysis.stage=enrichment_only`.
3. Dry-run `scripts/run.py --config <config>` and inspect table mappings and requested databases.
4. Execute with `--execute` and inspect `task_status.tsv`, mapping tables, complete results, and database-level statuses.

## Output and plotting contract

- Preserve every complete enrichment result table and the declared ORA universe. Default to `enrichment.universe_mode=multiple_testing_eligible`; use `tested` only when deliberately including independently filtered genes. GSEA uses all finite signed statistics, including those filtered from adjusted-P testing. Preserve supplied test-state fields when importing skill 11 results.
- Keep unmapped identifiers and resource versions auditable.
- Record individual databases as completed, empty, skipped, or failed; retain successful databases when another fails.
- Produce a compact overview plus separate database/direction plots using configurable top-N terms, wrapped labels, bounded dimensions, and pagination rather than one unbounded figure. Default to PDF; PNG is an explicitly selected alternative for standard/compact plots.
- Treat summary plots as views of the complete TSV results, never as substitutes for them.

Read [references/enrichment-design.md](references/enrichment-design.md) for resource and interpretation rules. Dry-run first and use `scripts/run_in_tmux.py` for long confirmed runs.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.

The standard enrichment run retains its existing lollipop/dot summaries. For richer presentation of completed results, use the table-only review modes below. Read [references/figure-style.md](references/figure-style.md) for the shared `paper_v1` contract.

Result handoff: start with `<output_dir>/RESULTS.md` from the execute runner. It links scientific deliverables and reports current run status and limitations; follow [output-layout.md](references/output-layout.md) for retained results and technical records.

Figure output: select one format per run, default PDF; follow [figure-output.md](references/figure-output.md). Never automatically export both PDF and PNG.


## Replot existing enrichment results

These presentation-only entrypoints use complete enrichment tables and a fresh output directory. They do not rerun DE, ORA or GSEA. Use an existing registered analysis R environment (normally 06-deg-analysis); do not install missing packages or substitute species resources. Copy [config.review.example.json](references/config.review.example.json), select one population/comparison, and explicitly select `task_id`/`source_table` if multiple input tasks would otherwise overlap.

- **Database atlas**, preferred for comprehensive review: `Rscript scripts/enrichment_review_atlas.R <config.json>`. Read [review-atlas.md](references/review-atlas.md). Produces one PDF book each for GSEA and ORA, content-sized pages, directional returned/finite/significant/selected/page counts, NES or fold enrichment, FDR, core/hit ratios and gene examples. Requested zero-significant and failed databases retain explicit pages. Supply the matching enrichment status, DE and SYMBOL/ENTREZID mapping tables.
- **Compact dots and exact-ID evidence matrix**: `Rscript scripts/enrichment_review_plot.R <config.json>`. Read [review-compact.md](references/review-compact.md). Optional lollipops use the same selection; default dot only. ORA Up/Down use separate positive-magnitude panels. Matrix colour is signed capped -log10 FDR, not effect size; grey nonsignificance and missing rows differ. PDF books default; PNG is an alternative.
- **GSEA information/core-gene matrix**, optional: `Rscript scripts/enrichment_review_dense.R <config.json>`. Read [review-dense.md](references/review-dense.md). Two paginated PDF books, complete memberships and unmapped-ID audit. Grey means not a leading-edge member, not unexpressed. Gene log2FC colours are explicitly clipped.
- **GO semantic modules**, optional: `Rscript scripts/enrichment_review_modules.R <config.json>`. Read [review-modules.md](references/review-modules.md). Requires GOSemSim, GO.db, the explicit species OrgDb, igraph and a declared Wang threshold/seed. Keeps every significant term and singleton, with independent GSEA/ORA and direction groups. Medoid labels are original GO terms; modules have no new P values or curated theme identities.

Default review display FDR is 0.05. Raising it is explicit exploratory display, with hollow points in compact plots. Redundancy filtering uses within-database/method/direction contribution-gene Jaccard only for display; it does not alter significance. Never pool symbol and Entrez identifier systems in overlap comparisons, infer cross-database equivalence, or interpret ORA input direction as pathway activation. These mapping-assisted modes currently require SYMBOL-to-ENTREZID mappings; other ID systems must be explicitly prepared, not guessed.

Inspect dense, sparse and empty/failure pages before delivery. Keep complete source tables unchanged and deliver RESULTS.md, selected values, selection/count audits and provenance. PDF book modes require Ghostscript (`gs`); atlas, dense and modules reject PNG requests explicitly. Run only the requested figure family, preserving one scientific format per run.
