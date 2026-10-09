# 导出未校正 counts 和按 cell ID 对齐的元数据，不修改源对象。
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: export_cnmf_input.R config.json destination")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "runtime.R"))
config <- read_skill_config(args[[1]])
obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE))
assay <- cfg_get(config, "input.assay", "RNA")
if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
if (tolower(assay) %in% c("sct", "integrated")) stop("Use uncorrected RNA raw counts")
Seurat::DefaultAssay(obj) <- assay
counts <- get_raw_counts(obj, assay)
if (!setequal(colnames(counts), colnames(obj))) stop("Counts cell IDs do not match object")
counts <- counts[, colnames(obj), drop = FALSE]
meta <- obj[[]][colnames(counts), , drop = FALSE]
columns <- unlist(lapply(c("sample", "condition", "cell_type", "batch"), function(role) cfg_get(config, paste0("metadata.", role))))
assert_metadata(obj, columns)
dest <- args[[2]]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
Matrix::writeMM(counts, file.path(dest, "counts.mtx"))
writeLines(rownames(counts), file.path(dest, "features.tsv"))
writeLines(colnames(counts), file.path(dest, "barcodes.tsv"))
utils::write.table(cbind(cell_id = rownames(meta), meta), file.path(dest, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
reduction <- cfg_get(config, "metadata.reduction")
if (!is.null(reduction)) {
  if (!reduction %in% names(obj@reductions)) stop("Configured reduction is absent: ", reduction)
  embedding <- Seurat::Embeddings(obj, reduction)[colnames(counts), , drop = FALSE]
  if (ncol(embedding) < 2L) stop("Embedding requires at least two dimensions")
  utils::write.table(cbind(cell_id = rownames(embedding), embedding[, 1:2, drop = FALSE]), file.path(dest, "embedding.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("Exported", ncol(counts), "cells and", nrow(counts), "genes\n")
