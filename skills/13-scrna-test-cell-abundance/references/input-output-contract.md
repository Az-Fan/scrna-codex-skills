# Input and output contract

## Input

Use exactly one of `input.object` (annotated Seurat RDS/QS; required for Milo) or `input.counts_table` (UTF-8 TSV). A count table requires the configured sample, condition, cell-type, covariate and count columns. Repeated rows are summed and missing sample-by-cell-type combinations become explicit zeros.

Every sample must map to one condition and one value for each sample-level covariate. Positive effects mean numerator greater than denominator on the recorded method-specific scale.

`analysis.denominator.mode=all_input_cells` uses every input cell. `selected_cell_types` requires `analysis.denominator.include` and recomputes sample totals after selection. Always state `analysis.denominator.description`. An EC-enriched object estimates relative EC-subtype composition, not the EC fraction or absolute EC number in the original tissue.

## Methods

`analysis.methods` accepts `propeller`, `sccomp`, `sccoda`, `dcats`, and `milo`. Each comparison-method task runs independently.

- sccomp requires CmdStanR and CmdStan. The registered environment supplies both and caches the first model compilation inside the pixi environment.
- scCODA runs through the named `sccoda` pixi environment. Automatic reference sensitivity adds stable, widely observed candidate references.
- Milo requires `input.object` and an explicit `method_options.milo.reduction`; the skill does not silently recompute or integrate the manifold.
- DCATS optionally accepts a square, labelled cell-type similarity matrix.

## Outputs

Root outputs include complete `sample_cell_counts.tsv`, `sample_cell_proportions.tsv`, `design_audit.tsv`, `cell_type_eligibility.tsv`, `task_status.tsv`, `all_method_results.tsv`, `significant_method_results.tsv`, `method_concordance.tsv`. Technical records are `_provenance/session_info.txt`, `_provenance/run.log`, and `_provenance/run_manifest.json`.

Descriptive figures are bounded multipage PDFs: `sample_composition.pdf`, `cell_type_proportions_by_condition.pdf`, and `sample_proportion_heatmap.pdf`. Every sample-level point is labelled.

Under `comparisons/<comparison>/<method>/`, retain the official raw table, standardized `all_results.tsv`, `significant_results.tsv`, and `effect_summary.pdf`. Extra outputs include scCODA per-reference tables and posterior diagnostics, sccomp draws and optional fit, DCATS coefficients, and Milo neighborhood tables, object and DA plots.

Standardized scales remain explicit: propeller descriptive log2 mean-proportion ratio, sccomp logit composition effect, scCODA selected log2 fold change relative to its reference, DCATS beta-binomial log odds, and Milo neighborhood-count log2 fold change. Never average these scales.
