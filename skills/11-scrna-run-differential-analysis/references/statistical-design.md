# Statistical design

## Checklist

Record the biological replicate, numerator, denominator, population column, selected populations, minimum cells per sample-population, excluded samples, covariates, pairing or blocking, count assay, model formula, multiple-testing method, and whether inference is formal or exploratory.

Use raw counts aggregated by `sample × population` for formal condition inference. Require both groups after minimum-cell filtering and at least `analysis.min_samples_per_group` independent samples per group. Filter genes by both total count and expression of at least `analysis.min_count_per_sample` counts in `analysis.min_samples_expressed` samples. Two samples per group is a technical minimum, not a guarantee of useful power; prefer three or more.

Keep one condition and one value for every modeled covariate per sample. Use a paired formula such as `~ patient + condition` only when patients contribute the required conditions. Reject a rank-deficient model rather than silently dropping terms. Do not include a batch term that is perfectly confounded with condition.

Use `seurat_*` methods only as exploratory sensitivity analyses or when no valid sample-level design exists. Their adjusted P values do not account for biological replication.

Positive `log2FoldChange` always means numerator exceeds denominator. Thresholds annotate results but do not determine which rows are retained in `all_genes.tsv`.

Use apeglm-shrunken log2 fold changes by default for stable effect-size reporting and visualization. `analysis.lfc_shrink=true` is a strict dependency contract: fail if `apeglm` is unavailable or the coefficient cannot be identified. Only fall back when `analysis.allow_unshrunk_lfc=true` is explicitly set, and record the fallback in `effect_size_audit.tsv` and result columns. Setting `analysis.lfc_shrink=false` is preferable when unshrunk effects are intentionally reported. Replace only `log2FoldChange` and `lfcSE` with apeglm estimates; retain Wald P values, adjusted P values, and signed Wald `stat` from the unshrunk DESeq2 result. Use that Wald statistic, not shrunken fold change, as the preferred GSEA ranking metric.
