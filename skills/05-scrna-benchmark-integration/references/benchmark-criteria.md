# Benchmark criteria

## Method registry

Supported method names are `none`, `harmony`, `rpca`, `scvi`, `scanvi`, `bbknn`, and `precomputed`. Expand each array-valued parameter in `parameter_grid` into separate scenarios. `none` is injected when omitted. Treat scANVI as supervised and require an explicit label field. Treat BBKNN as graph correction and exclude it from embedding-only metrics unless the implementation explicitly supports its graph.

Precomputed scenarios include the reduction name in their ID. An optional method `id` can distinguish other repeated method specifications. Expanded IDs must remain unique after sanitization; duplicate parameter values or IDs are rejected before outputs are written.

Run R-native adapters in the R process. Run scVI, scANVI, BBKNN, and scIB-metrics through `benchmark.python_argv_prefix`. Use `['python3']` for an environment containing the required packages or an argv prefix such as `['pixi','run','-e','scvi','python']`; never encode it as a shell command string. The Python stage runs on CPU unless the execution environment and user configuration explicitly provide otherwise.

## Metadata roles

- `sample`: biological or library replicate identifier.
- `batch_variables`: one or more technical variables evaluated independently.
- `biological_labels`: one or more trusted annotations used to quantify biological conservation.
- `condition`: biological comparison variable used for confounding and preservation audits.

Reject a batch variable with fewer than two observed levels. Report cross-tabs and Cramer's V for batch-condition association. A perfect batch-condition mapping blocks automatic recommendation but does not block descriptive benchmarking.

## Metric groups

Batch-removal choices may include `ilisi`, `batch_asw`, `pcr_comparison`, `graph_connectivity`, and `kbet`. Biological-conservation choices may include `clisi`, `label_asw`, `isolated_labels`, `nmi`, and `ari`. Exact availability depends on representation type, labels, graph availability, sample size, and installed packages.

Write one row per scenario, batch variable, biological label, and metric. Record `completed`, `skipped_missing_label`, `skipped_incompatible_representation`, `missing_dependency`, or `failed`; never silently drop a requested metric.

The current scIB adapter requires an explicit biological label even for a batch-only request; absence is recorded as `skipped_missing_label`. Preserve `raw_value` (scIB's result before across-scenario min–max scaling) alongside `scaled_value`; intrinsic scIB metric transformations are retained. The compatibility column `value` equals the latter and `value_scale=scib_min_max_across_scenarios` describes its scale. Scaling is relative to the scenarios within that batch/label evaluation, so scaled values cannot establish comparability between separate runs.

## Ranking

Keep batch-removal and biological-conservation summaries separate. If `scoring.enabled` is true, calculate the configured weighted score only after scaling comparable completed metrics. Do not impute missing metrics. Also report Pareto-efficient scenarios so a weighted score is never the sole selection rule.

Automatic recommendation requires both metric groups, a completed uncorrected baseline, and complete finite evidence for every requested scenario/batch/label/metric combination. Missing, failed, duplicated, or skipped requested evidence makes the decision unresolved; partial tables remain available for review. The machine-readable status includes per-scenario expected and completed metric-row counts.

## Plots

Supported plot selections are `score_barplot`, `score_heatmap`, `ranking_plot`, `metric_tradeoff`, `umap_by_batch`, `umap_by_sample`, `umap_by_condition`, `umap_by_label`, `marker_dotplot`, and `program_retention`. Marker and program plots require user-provided gene sets and normalized uncorrected expression.

## Required outputs

Write `method_runs.tsv`, `metric_results_long.tsv`, `method_summary.tsv`, `method_ranking.tsv`, `skipped_metrics.tsv`, `design_confounding.tsv`, selected figures, a derivative benchmark object or embedding bundle, `recommendation.md`, machine-readable `recommendation_status.json`, and `_provenance/run_manifest.json`. The decision status must be `unresolved` when no requested metric completed or when batch and condition are perfectly confounded.
