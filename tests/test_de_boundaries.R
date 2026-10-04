source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")
reject <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  stopifnot(is.character(message), grepl(pattern, message, fixed = TRUE))
}
minimums <- c(min_samples_per_group = 2, min_cells_per_sample_population = 1,
              min_total_count = 0, min_count_per_sample = 0, min_samples_expressed = 1)
for (field in names(minimums)) {
  for (bad in list(minimums[[field]] - 1, 1.5, TRUE, "2", NA_real_, Inf, c(2, 3), NULL)) {
    config <- list(analysis = setNames(list(bad), field))
    reject(validate_de_thresholds(config), paste0("analysis.", field))
  }
  validate_de_thresholds(list(analysis = setNames(list(minimums[[field]]), field)))
}
for (field in c("lfc_threshold", "padj_threshold")) {
  for (bad in list(-1, TRUE, "0.05", NA_real_, NaN, Inf, c(.1, .2), NULL))
    reject(validate_de_thresholds(list(analysis = setNames(list(bad), field))), paste0("analysis.", field))
}
for (bad in c(0, 1)) reject(validate_de_thresholds(list(analysis = list(padj_threshold = bad))), "padj_threshold")
stopifnot(validate_de_thresholds(list(analysis = list(lfc_threshold = 0)))$lfc == 0)
for (design in c("~ sex * condition", "~ sex + condition + sex:condition", "~ (sex + condition)^2",
                 "~ condition - 1", "~ 0 + condition", "~ log(age) + condition", "~ .", "~ sex"))
  reject(validate_de_design(list(analysis = list(design = design)), "condition", c("sex", "age")), "design")
stopifnot(identical(all.vars(validate_de_design(list(analysis = list(design = "~ patient + condition")),
  "condition", "patient")), c("patient", "condition")))
validate_de_design(list(analysis = list(design = "~ `batch*id` + `condition+type`")), "condition+type", "batch*id")

# Real fitted DESeq2 models: special labels and swapped directions must map to
# the exact coefficient used by results(), not another valid coefficient.
set.seed(124)
counts <- matrix(rnbinom(120 * 8, mu = 100, size = 12), 120, 8,
                 dimnames = list(paste0("g", 1:120), paste0("s", 1:8)))
for (labels in list(c("Ctrl", "STZ+Drug"), c("Low Dose", "High Dose"), c("WT/WT", "KO/KO"), c("A.1", "B(2)"))) {
  for (direction in list(labels, rev(labels))) {
    comparison <- list(denominator = direction[1], numerator = direction[2])
    md <- data.frame(condition = factor(rep(labels, each = 4), levels = direction),
                     sex = factor(rep(c("F", "M"), 4)), row.names = colnames(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, md, ~ sex + condition)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, fitType = "mean")
    coefficient <- resolve_apeglm_coefficient(dds, "condition", comparison)
    by_contrast <- DESeq2::results(dds, contrast = c("condition", direction[2], direction[1]))
    by_name <- DESeq2::results(dds, name = coefficient)
    stopifnot(isTRUE(all.equal(by_contrast$log2FoldChange, by_name$log2FoldChange)),
              isTRUE(all.equal(by_contrast$pvalue, by_name$pvalue)))
    shrunken <- DESeq2::lfcShrink(dds, coef = coefficient, type = "apeglm")
    stopifnot(identical(rownames(shrunken), rownames(by_contrast)), all(is.finite(shrunken$log2FoldChange)))
  }
}
cat("PASS strict R numeric boundaries, additive estimand designs, and real DESeq2/apeglm special-label coefficient equivalence\n")
