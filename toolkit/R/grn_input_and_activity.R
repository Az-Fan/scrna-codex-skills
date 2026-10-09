# GRN 输入聚合和单细胞 AUCell 评分；不修改源 Seurat 对象。
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: grn_input_and_activity.R config.json prepare|score")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "runtime.R"))
config <- read_skill_config(args[[1]])
out <- cfg_get(config, "output_dir", required = TRUE)
dest <- file.path(out, "_provenance", "grn_input")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
read_ids <- function(path, column = 1L) {
  tab <- utils::read.delim(path, header = FALSE, colClasses = "character", check.names = FALSE, na.strings = NULL)
  if (column > ncol(tab)) stop("Identifier column absent")
  ids <- tab[[column]]
  if (anyDuplicated(ids) || any(!nzchar(trimws(ids)))) stop("Empty or duplicated identifiers")
  ids
}
read_counts <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con))
  as(Matrix::readMM(con), "CsparseMatrix")
}
if (args[[2]] == "prepare") {
  if (cfg_get(config, "input.type", "seurat") == "seurat") {
    obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE))
    counts <- get_raw_counts(obj, cfg_get(config, "input.assay", "RNA"))
    if (!setequal(colnames(counts), colnames(obj))) stop("Counts and object cells disagree")
    counts <- counts[, colnames(obj), drop = FALSE]
    meta <- obj[[]][colnames(counts), , drop = FALSE]
  } else {
    genes <- read_ids(cfg_get(config, "input.features"), cfg_get(config, "input.feature_column", 1L))
    cells <- read_ids(cfg_get(config, "input.barcodes"))
    counts <- read_counts(cfg_get(config, "input.counts"))
    if (cfg_get(config, "input.orientation", "genes_by_cells") == "cells_by_genes") counts <- Matrix::t(counts)
    if (!identical(dim(counts), c(length(genes), length(cells)))) stop("Matrix dimensions/orientation disagree with IDs")
    dimnames(counts) <- list(genes, cells)
    tab <- utils::read.delim(cfg_get(config, "input.metadata"), colClasses = "character", check.names = FALSE, na.strings = NULL)
    if (anyDuplicated(tab[[1]]) || !setequal(tab[[1]], cells)) stop("Metadata IDs must exactly match barcodes")
    rownames(tab) <- tab[[1]]
    meta <- tab[cells, -1, drop = FALSE]
  }
  # Seurat counts layers can be valid dense matrices; standardize storage before
  # slot validation and Matrix Market export without changing count values.
  counts <- as(Matrix::Matrix(counts, sparse = TRUE), "CsparseMatrix")
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Counts IDs must be unique")
  if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) stop("Require finite non-negative integer raw UMI counts")
  if (any(Matrix::colSums(counts) <= 0)) stop("Zero-count cells; review QC first")
  roles <- c("sample", "cell_type", "condition", "batch")
  columns <- unique(unlist(lapply(roles, function(r) cfg_get(config, paste0("metadata.", r)))))
  if (!all(columns %in% colnames(meta))) stop("Configured metadata column absent")
  for (column in columns) if (anyNA(meta[[column]]) || any(!nzchar(trimws(as.character(meta[[column]]))))) stop("Missing metadata: ", column)
  sample <- cfg_get(config, "metadata.sample")
  for (role in c("condition", "batch")) {
    column <- cfg_get(config, paste0("metadata.", role))
    if (!is.null(column) && any(vapply(split(meta[[column]], meta[[sample]]), function(x) length(unique(x)) > 1L, logical(1)))) stop("A sample has multiple ", role, " values")
  }
  mode <- cfg_get(config, "inference.mode")
  if (mode == "metacell") {
    column <- cfg_get(config, "metacell.column")
    if (!column %in% colnames(meta) || anyNA(meta[[column]]) || any(!nzchar(trimws(as.character(meta[[column]]))))) stop("Missing metacell memberships")
    # 同一个聚类标签按样本及已有细胞类型拆分，避免混合不同生物学来源。
    strata <- unique(c(sample, cfg_get(config, "metadata.cell_type"), column))
    keys <- lapply(seq_len(nrow(meta)), function(i) vapply(strata, function(column) as.character(meta[[column]][i]), character(1)))
    key <- vapply(keys, function(x) jsonlite::toJSON(x, auto_unbox = FALSE), character(1))
    membership <- sprintf("metacell_%06d", match(key, unique(key)))
  } else membership <- colnames(counts)
  ids <- unique(membership)
  aggregate_matrix <- Matrix::sparseMatrix(i = seq_along(membership), j = match(membership, ids), x = 1,
                                          dims = c(length(membership), length(ids)))
  inference <- counts %*% aggregate_matrix
  dimnames(inference) <- list(rownames(counts), ids)
  sizes <- tabulate(match(membership, ids))
  audit <- data.frame(inference_id = ids, n_cells = sizes)
  for (column in columns) {
    audit[[paste0(column, "_purity")]] <- vapply(ids, function(id) {
      # Convert factors to labels so unused levels never contribute 0 * log(0).
      p <- prop.table(table(as.character(meta[[column]][membership == id]))); if (length(p) == 1L) return(1)
      # 用全输入类别数作为熵的统一分母。
      1 + sum(p * log(p)) / log(length(unique(meta[[column]])))
    }, numeric(1))
  }
  write_tsv(data.frame(cell_id = colnames(counts), inference_id = membership, meta, check.names = FALSE), "cell_membership.tsv")
  write_tsv(audit, "inference_unit_audit.tsv")
  if (any(sizes < cfg_get(config, "metacell.min_cells", 1L))) stop("Inference units smaller than metacell.min_cells; review memberships")
  if (ncol(inference) < 3L) stop("Need at least three inference units")
  sums <- Matrix::rowSums(inference)
  variable <- Matrix::rowSums(inference^2) - sums^2 / ncol(inference) > 1e-8
  write_tsv(data.frame(gene = rownames(counts), inference_status = ifelse(variable, "retained", "constant_or_zero")), "feature_status.tsv")
  if (sum(variable) < 3L) stop("Need at least three variable genes for GRN inference")
  Matrix::writeMM(counts, file.path(dest, "single_cell_counts.mtx"))
  Matrix::writeMM(inference[variable, , drop = FALSE], file.path(dest, "inference_counts.mtx"))
  writeLines(rownames(counts), file.path(dest, "genes.tsv"))
  writeLines(rownames(counts)[variable], file.path(dest, "inference_genes.tsv"))
  writeLines(colnames(counts), file.path(dest, "cells.tsv"))
  writeLines(ids, file.path(dest, "inference_ids.tsv"))
  write_tsv(data.frame(cell_id = rownames(meta), meta, check.names = FALSE), "cell_metadata.tsv")
  file.copy(file.path(out, "cell_metadata.tsv"), file.path(dest, "cell_metadata.tsv"))
  grDevices::png(file.path(out, "inference_unit_sizes.png"), width = 1600, height = 1000, res = 180)
  hist(sizes, main = "Cells per GRN inference unit", xlab = "Cell count", col = "grey60", breaks = "Sturges")
  grDevices::dev.off()
} else if (args[[2]] == "score") {
  counts <- read_counts(file.path(dest, "single_cell_counts.mtx"))
  dimnames(counts) <- list(readLines(file.path(dest, "genes.tsv")), readLines(file.path(dest, "cells.tsv")))
  edge <- utils::read.delim(file.path(out, "regulon_targets.tsv"), colClasses = "character", check.names = FALSE, na.strings = NULL)
  gene_sets <- lapply(split(edge$target, edge$regulon), unique)
  if (!length(gene_sets)) stop("No retained regulons to score")
  min_genes <- cfg_get(config, "regulons.min_genes", 20L)
  if (any(lengths(gene_sets) < min_genes) || any(!unlist(gene_sets) %in% rownames(counts))) stop("Regulon targets disagree with the scored gene universe")
  saveRDS(gene_sets, file.path(out, "regulon_list.rds"))
  seed <- cfg_get(config, "grn.seed", 777L)
  batch_size <- cfg_get(config, "aucell.batch_size", 500L)
  rank_max <- floor(nrow(counts) * cfg_get(config, "aucell.rank_fraction", 0.05))
  if (rank_max < 2L) stop("aucell.rank_fraction gives fewer than two ranked genes; choose explicitly for small inputs")
  batches <- split(seq_len(ncol(counts)), ceiling(seq_len(ncol(counts)) / batch_size))
  score <- lapply(seq_along(batches), function(i) {
    # 固定并记录随机种子；相同表达值的随机排名依赖批次设置。
    set.seed(seed + i - 1L)
    ranking <- AUCell::AUCell_buildRankings(counts[, batches[[i]], drop = FALSE], plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(gene_sets, ranking, aucMaxRank = rank_max, normAUC = TRUE, verbose = FALSE)
    AUCell::getAUC(auc)
  })
  score <- t(do.call(cbind, score))
  score <- score[colnames(counts), , drop = FALSE]
  if (any(!is.finite(score)) || any(score < 0)) stop("Non-finite or negative AUCell scores")
  write_tsv(data.frame(cell_id = rownames(score), score, check.names = FALSE), "regulon_activity.tsv")
  saveRDS(score, file.path(out, "regulon_activity.rds"))
  meta <- utils::read.delim(file.path(dest, "cell_metadata.tsv"), colClasses = "character", check.names = FALSE, na.strings = NULL)
  if (anyDuplicated(meta$cell_id) || !setequal(meta$cell_id, rownames(score))) stop("AUCell metadata cell IDs disagree")
  meta <- meta[match(rownames(score), meta$cell_id), , drop = FALSE]
  cols <- unique(unlist(lapply(c("sample", "cell_type", "condition", "batch"), function(r) cfg_get(config, paste0("metadata.", r)))))
  summaries <- stats::aggregate(as.data.frame(score), by = meta[, cols, drop = FALSE], FUN = mean)
  sizes <- stats::aggregate(rep(1L, nrow(meta)), by = meta[, cols, drop = FALSE], FUN = sum)
  summaries$n_cells <- sizes$x
  write_tsv(summaries, "regulon_activity_by_sample.tsv")
  jsonlite::write_json(list(rank_max = rank_max, n_genes = nrow(counts), seed = seed, batch_size = batch_size,
                            package_version = as.character(utils::packageVersion("AUCell"))), file.path(out, "_provenance", "aucell_record.json"), auto_unbox = TRUE, pretty = TRUE)
} else stop("Unknown stage")
writeLines(capture.output(sessionInfo()), file.path(out, "_provenance", "r_session_info.txt"))
