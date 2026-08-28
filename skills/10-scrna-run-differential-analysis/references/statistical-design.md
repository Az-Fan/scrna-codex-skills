# Statistical design

## Checklist

Record the biological replicate, numerator, denominator, population column, selected populations, minimum cells per sample-population, excluded samples, covariates, pairing or blocking, count assay, model formula, multiple-testing method, and whether inference is formal or exploratory.

Use raw counts aggregated by `sample × population` for formal condition inference. Require both groups after minimum-cell filtering and at least `analysis.min_samples_per_group` independent samples per group. Filter genes by both total count and expression of at least `analysis.min_count_per_sample` counts in `analysis.min_samples_expressed` samples. Two samples per group is a technical minimum, not a guarantee of useful power; prefer three or more.

Keep one condition and one value for every modeled covariate per sample. Use a paired formula such as `~ patient + condition` only when patients contribute the required conditions. Reject a rank-deficient model rather than silently dropping terms. Do not include a batch term that is perfectly confounded with condition.

Use `seurat_*` methods only as exploratory sensitivity analyses or when no valid sample-level design exists. Their adjusted P values do not account for biological replication.

Positive `log2FoldChange` always means numerator exceeds denominator. Thresholds annotate results but do not determine which rows are retained in `all_genes.tsv`.

Use apeglm-shrunken log2 fold changes by default for stable effect-size reporting and visualization. Keep Wald P values for inference and use the signed Wald statistic, not shrunken fold change, as the preferred GSEA ranking metric.
