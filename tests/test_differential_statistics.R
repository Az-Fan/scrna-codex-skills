source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")
root <- tempfile("scrna-statistics-"); dir.create(root)
comparison <- list(id = "B_vs_A", numerator = "B", denominator = "A")
thresholds <- list(padj = .05, lfc = .25, min_total_count = 0, min_count_per_sample = 0, min_samples_expressed = 1)

# Real fitted model with low-expression genes, a zero row and a Cook outlier.
set.seed(42)
n <- 1600L
mu <- c(rep(100, 400), rep(.3, 1200))
means <- matrix(mu, n, 8)
means[2:151, 5:8] <- means[2:151, 5:8] * 1.7
counts <- matrix(rnbinom(n * 8, mu = means, size = 20), n, 8,
  dimnames = list(paste0("gene", seq_len(n)), paste0("sample", 1:8)))
counts[1, ] <- c(1000000, 20, 20, 20, 20, 20, 20, 20)
counts[n, ] <- 0
obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
obj$sample <- colnames(counts); obj$condition <- rep(c("A", "B"), each = 4)
config <- list(analysis = list(lfc_shrink = FALSE))
# Figure rendering is covered separately; this test exercises the fitted statistics.
plot_pseudobulk <- function(...) invisible(NULL)
res <- run_pseudobulk(obj, obj[[]], "RNA", "sample", "condition", character(), comparison, thresholds, config, root)
res <- annotate_de_result(res, "all", comparison, "pseudobulk_deseq2", thresholds, obj[[]], "sample", "condition")
audit <- jsonlite::read_json(file.path(root, "deseq2_results_audit.json"), simplifyVector = TRUE)
stopifnot(audit$alpha == .05, any(res$filter_reason == "independent_filtering"),
  any(res$filter_reason == "cooks_outlier"), res$filter_reason[res$gene == paste0("gene", n)] == "zero_counts")
stopifnot(all(res$wald_tested == is.finite(res$pvalue)), all(res$multiple_testing_eligible == is.finite(res$padj)),
  all(res$significance[res$filter_reason == "independent_filtering"] == "Filtered"),
  all(res$significance[!res$tested] == "Not_tested"))
# An independently constructed DESeq2 result must agree at the requested FDR.
pb <- readRDS(file.path(root, "pseudobulk_data.rds"))
dds <- DESeq2::DESeqDataSetFromMatrix(as.matrix(pb$counts), pb$coldata, ~ condition)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
expected <- DESeq2::results(dds, contrast = c("condition", "B", "A"), alpha = .05)
stopifnot(isTRUE(all.equal(res$pvalue, expected$pvalue)), isTRUE(all.equal(res$padj, expected$padj)),
  isTRUE(all.equal(audit$filter_threshold, unname(S4Vectors::metadata(expected)$filterThreshold))))
stopifnot(pb$diagnostic_transform == "vst",
  isTRUE(all.equal(pb$diagnostic_expression, SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)))),
  all(pb$counts == counts), !any(res$low_replication_warning))
few <- annotate_de_result(res, "all", comparison, "pseudobulk_deseq2", thresholds, obj[[]][c(1, 2, 5, 6), ], "sample", "condition")
stopifnot(all(few$low_replication_warning), all(few$inference_qualification == "exploratory_low_confidence"))
config$analysis$pca_transform <- "log2_normalized"
alternative_dir <- file.path(root, "alternative"); dir.create(alternative_dir)
alternative <- run_pseudobulk(obj, obj[[]], "RNA", "sample", "condition", character(), comparison, thresholds, config, alternative_dir)
alternative_pb <- readRDS(file.path(alternative_dir, "pseudobulk_data.rds"))
stopifnot(alternative_pb$diagnostic_transform == "log2_normalized",
  isTRUE(all.equal(alternative_pb$diagnostic_expression, log2(alternative_pb$normalized_counts + 1))),
  identical(res$pvalue, alternative$pvalue), identical(res$padj, alternative$padj))

# TSV handoff must retain tested-but-filtered rows and their exact reason.
path <- file.path(root, "de.tsv"); write_tsv(res, path)
imported <- standardize_de_table(read_de_table(path), list())
stopifnot(identical(imported$tested, res$tested), identical(imported$multiple_testing_eligible, res$multiple_testing_eligible),
  identical(imported$filter_reason, res$filter_reason), identical(imported$significance, res$significance))
default <- ora_background(imported, list())
wider <- ora_background(imported, list(enrichment = list(universe_mode = "tested")))
stopifnot(setequal(default$genes, imported$gene[is.finite(imported$padj)]),
  length(wider$genes) > length(default$genes),
  all(imported$gene[imported$filter_reason == "independent_filtering"] %in% wider$genes))
external <- data.frame(gene = c("a", "b"), log2FoldChange = c(1, -1), pvalue = c(.02, NA), padj = c(NA, NA))
external <- standardize_de_table(external, list())
stopifnot(identical(external$significance, c("Filtered", "Not_tested")),
  external$filter_reason[1] == "adjusted_p_unavailable", !length(ora_background(external, list())$genes))
rank_only <- standardize_de_table(data.frame(gene = "ranked", stat = 3), list())
stopifnot(rank_only$tested, !rank_only$multiple_testing_eligible, rank_only$stat == 3)
bad <- res; bad$multiple_testing_eligible[!bad$tested] <- TRUE
failure <- tryCatch(standardize_de_table(bad, list()), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("Inconsistent", failure))
bad <- res; bad$tested <- "yes"
failure <- tryCatch(standardize_de_table(bad, list()), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("Invalid boolean", failure))
prefilter <- utils::read.delim(file.path(root, "gene_filter_audit.tsv"))
stopifnot(nrow(prefilter) == n, all(prefilter$entered_deseq2))
unlink(root, recursive = TRUE)
cat("PASS real DESeq2 alpha/filter states, VST and explicit log transform, low replication, TSV handoff and ORA universes\n")
