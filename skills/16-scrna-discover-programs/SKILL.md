---
name: 16-scrna-discover-programs
description: Discover de novo gene expression programs with cNMF from uncorrected raw scRNA-seq counts in Seurat RDS/QS or a declared matrix bundle. Scan ranks, review stability and reconstruction error, reuse a completed discovery for selected consensus ranks, and export cell usage, gene spectra, top genes, sample summaries and optional local-GMT enrichment. Use for discovering coexpressed identity or state programs; use 10-scrna-score-programs for scoring predefined signatures.
---

# Discover scRNA expression programs

cNMF separates non-negative expression into gene programs and cell usage. A cell may use several programs; programs need not coincide with clusters or represent disease effects.

## Workflow

1. Confirm the object or matrix bundle, raw-count provenance, biological sample column, population scope, candidate ranks, HVG count, repetitions, seed and available CPU/RAM. Input should already have passed QC. Use the population requested; do not silently subset, balance or correct batches.
2. Read [references/discovery-contract.md](references/discovery-contract.md) for formats, rank review, density filtering and result interpretation. Adapt [references/config.example.json](references/config.example.json).
3. Run `scripts/check_dependencies.py --config <config>` and `scripts/run.py --config <config>` in the foreground. Review the resolved interpreters and parameters. For long execution read [references/long-running-execution.md](references/long-running-execution.md) and use the bundled tmux supervisor with `--execute`.
4. Default `workflow.action=discover` runs preparation, replicate factorizations, combination and k-selection diagnostics, then reports `awaiting_k_confirmation`. Show `k_selection.png` and complete `k_selection_stats.tsv`; discuss stability/error tradeoffs and obtain the selected rank(s).
5. Change `workflow.action` to `consensus`, document `workflow.selection_reason`, and list `cnmf.consensus_k`. Keep the discovery input, preparation parameters and output directory. Factorization is not repeated.
6. When the user already specified ranks to export, `workflow.action=run` performs discovery and those explicit consensuses in one execution. Never replace unspecified ranks with highest silhouette or lowest error.
7. Review `task_status.tsv`, top genes, density-filtering clustergrams, sample/population usage summaries and existing-embedding plots. Fill `program_review.tsv` with candidate labels and evidence before describing programs biologically. Multi-k exports include all-pair usage correlations; matching program numbers do not imply matching programs.

## Scientific constraints

- Require `input.counts_source=raw_umi`. Use uncorrected RNA counts, never normalized, scaled, integrated or SCT values.
- Preserve source objects and cell IDs. Resolve duplicate genes explicitly. Zero-count cells stop execution; zero-expression genes are omitted with a complete feature audit.
- Align metadata by cell ID rather than row position. Join Seurat v5 counts layers in memory.
- Export raw and row-normalized fractional usage. Gene score coefficients can be signed. Fractional usage is compositional and depends on k.
- Check stress, mitochondrial, cell-cycle, sample and batch patterns before naming biological states. Condition differences need biological replicates; this executor produces descriptive summaries and optional gene-set ORA, not formal condition tests or causal TF inference.
- Discovery refuses existing prepared artifacts. Consensus verifies source/parameter/version/artifact integrity and never overwrites exported k directories.
- Missing scientific runtimes stop execution; analysis does not install dependencies or modify environments.

## Outputs and handoff

Start with `<output_dir>/RESULTS.md`. Primary results include rank diagnostics, metadata, input/feature audits and per-k usage, spectra, top genes, program GMT, sample summaries, review tables and PNG figures. Optional local GMT produces complete program ORA, coverage and tested-universe files. Native cNMF replicates and prepared data stay under `_provenance/` with logs, version/binding records and workflow state; retain these for consensus reuse.

Use `programs.gmt` with 10 only when deliberately scoring the discovered signatures; those scores are not original NMF usage. Read [references/output-layout.md](references/output-layout.md) for result organization.
