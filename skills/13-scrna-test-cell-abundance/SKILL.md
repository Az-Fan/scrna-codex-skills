---
name: 13-scrna-test-cell-abundance
description: Test replicated scRNA-seq cell-type composition or local cell-state abundance changes with propeller, sccomp, scCODA, DCATS, or Milo, while auditing the denominator, sample-level design, compositional interpretation, and method-specific evidence. Use after QC and annotation when comparing cell populations across conditions; do not use for gene-expression differential analysis or absolute tissue cell counting.
---

# Test scRNA Cell Abundance

Choose the estimand before the method. Cells are observations; samples or donors are the inferential replicates.

## Workflow

1. Confirm the input object or sample-by-cell-type count table, sample, condition, cell-type labels, covariates, biological replicates, and contrast direction.
2. Explicitly define the denominator as all input cells or a named set of parent populations. State what one reported proportion means biologically.
3. Read [references/methods-and-interpretation.md](references/methods-and-interpretation.md) and choose methods that answer the actual question. Do not combine method-specific effect sizes as though they shared one scale.
4. Copy [references/config.example.json](references/config.example.json), then dry-run `scripts/run.py --config <config>`.
5. Review validation messages and the resolved `07-cell-abundance` executor. Execute only after the denominator, design, comparison, and Milo reduction or scCODA reference strategy are explicit.
6. Inspect `design_audit.tsv`, `cell_type_eligibility.tsv`, `task_status.tsv`, complete results, sample-level plots, model diagnostics, and reference-sensitivity results before interpreting findings.

Use `scripts/run_in_tmux.py` for sccomp, scCODA, Milo, or multi-method runs that may outlive the session.

## Method routing

- Use `sccomp` as the robust annotated-composition model when outliers, overdispersion, or differential variability matter.
- Use `propeller` as a fast and transparent replicated baseline or independent sensitivity analysis.
- Use `sccoda` when a joint Bayesian compositional model is desired, especially with few replicates. Its result is relative to a reference population; keep reference sensitivity enabled when the reference is not biologically predeclared.
- Use `milo` for condition-associated neighborhoods within or across discrete annotations. Require an explicit existing reduction and audit how that representation was constructed.
- Use `dcats` primarily when annotation uncertainty can be represented by a similarity matrix or when its beta-binomial formulation is explicitly requested.

## Guardrails

- Never use cell-level Fisher, chi-square, Wilcoxon, or logistic tests as formal replicated abundance inference.
- Never describe captured scRNA-seq proportions as absolute tissue abundance without an external absolute-counting assay.
- Retain zero counts. Mark untestable tasks rather than silently deleting samples or populations.
- Stop rank-deficient, condition-confounded, or insufficient-replicate designs.
- Do not silently use an integrated embedding for Milo. Report the named reduction, dimensions, neighbors, and correction history.
- Preserve all method results. Thresholds annotate significance or credibility; they do not truncate complete tables.
- Interpret scCODA against its recorded reference. Treat sccomp posterior FDR, scCODA inclusion probability, propeller/DCATS adjusted P values, and Milo spatial FDR as different evidence types.
- Do not infer proliferation, death, migration, or recruitment mechanisms from relative abundance alone.

Read [references/input-output-contract.md](references/input-output-contract.md) for schemas, outputs, statuses, and plotting limits.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
