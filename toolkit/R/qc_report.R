args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript qc_report.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
config <- read_skill_config(args[[1]])
obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE), "auto")
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
assert_metadata(obj, sample_col)
meta <- obj[[]]
if (!"percent.mt" %in% colnames(meta)) {
  pattern <- cfg_get(config, "qc.mito_pattern")
  if (is.null(pattern)) {
    upper <- sum(grepl("^MT-", rownames(obj))); lower <- sum(grepl("^mt-", rownames(obj)))
    pattern <- if (upper >= lower && upper > 0) "^MT-" else if (lower > 0) "^mt-" else NULL
  }
  if (!is.null(pattern)) {
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = pattern)
    meta <- obj[[]]
  }
}
metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt"), colnames(meta))
if (!length(metrics)) stop("No standard QC metrics found")
summary <- do.call(rbind, lapply(split(meta, meta[[sample_col]]), function(x) {
  values <- unlist(lapply(metrics, function(m) c(median = median(x[[m]], na.rm = TRUE), q01 = quantile(x[[m]], .01, na.rm = TRUE), q99 = quantile(x[[m]], .99, na.rm = TRUE))))
  data.frame(cells = nrow(x), as.list(values), check.names = FALSE)
}))
summary[[sample_col]] <- rownames(summary); rownames(summary) <- NULL
out <- prepare_output(config)
summary_path <- file.path(out, "qc_summary_by_sample.tsv")
decision_path <- file.path(out, "qc_decision_template.tsv")
utils::write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(rule = character(), action = character(), reason = character(), approved = logical()), decision_path, sep = "\t", quote = FALSE, row.names = FALSE)
pdf(file.path(out, "qc_distributions.pdf"), width = 10, height = 4 * length(metrics)); par(mfrow = c(length(metrics), 1), mar = c(8, 4, 2, 1)); for (m in metrics) boxplot(meta[[m]] ~ meta[[sample_col]], las = 2, ylab = m, xlab = ""); dev.off()
write_run_manifest(config, "review-scrna-qc", out, c(summary_path, decision_path, file.path(out, "qc_distributions.pdf")), "No cells were filtered")
