args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript calculate_metrics.R config.json")

suppressPackageStartupMessages({
  library(Matrix)
  library(celda)
  library(Seurat)
  library(ggplot2)
  library(SingleCellExperiment)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("The existing pixi environment must provide jsonlite")
config <- jsonlite::read_json(args[[1]], simplifyVector = FALSE)
getv <- function(path, default = NULL) {
  value <- config
  for (key in strsplit(path, "\\.")[[1]]) {
    if (is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

starsolo_dir <- normalizePath(getv("input.starsolo_dir"), mustWork = TRUE)
gtf_file <- normalizePath(getv("input.gtf_file"), mustWork = TRUE)
out_root <- getv("output_dir")
if (is.null(out_root) || !nzchar(out_root)) stop("output_dir is required")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
min_genes_hq <- as.integer(getv("parameters.min_genes_hq", 500))
mito_pattern <- getv("parameters.mito_pattern", "^mt-")
artificial_fraction <- as.numeric(getv("parameters.doublet_artificial_fraction", 0.2))
set.seed(as.integer(getv("parameters.seed", 1)))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "doublet_finder.R"))

hbb_genes_mouse <- c(
  "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh1",
  "Hbb-y", "Hbb-bh2", "Hbb-b1", "Hbb-b2", "Hbm"
)
to_mouse_symbol <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

gtf <- readLines(gtf_file)
gtf <- gtf[!grepl("^#", gtf)]
fields <- strsplit(gtf, "\t", fixed = TRUE)
fields <- fields[vapply(fields, length, integer(1)) >= 9]
genes <- fields[vapply(fields, function(x) x[[3]] == "gene", logical(1))]
chrY <- genes[vapply(genes, function(x) x[[1]] %in% c("Y", "chrY"), logical(1))]
chrY_genes <- unique(vapply(chrY, function(x) sub('.*gene_name "([^"]+)".*', "\\1", x[[9]]), character(1)))

read_star_matrix <- function(directory, matrix_name = "matrix.mtx.gz") {
  counts <- as(Matrix::readMM(file.path(directory, matrix_name)), "CsparseMatrix")
  features <- read.delim(file.path(directory, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(file.path(directory, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- make.unique(features$V2)
  colnames(counts) <- barcodes$V1
  list(counts = counts, features = features)
}

samples <- getv("samples")
for (sample in samples) {
  sample_id <- sample$sample_id
  batch_id <- sample$batch_id
  message("Processing sample: ", sample_id)
  gene_dir <- file.path(starsolo_dir, sample_id, "Solo.out/Gene/filtered")
  velo_dir <- file.path(starsolo_dir, sample_id, "Solo.out/Velocyto/filtered")
  out_dir <- file.path(out_root, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gene_data <- read_star_matrix(gene_dir)
  counts <- gene_data$counts
  metadata <- data.frame(
    row.names = colnames(counts), sample_id = sample_id, batch_id = batch_id,
    n_genes = Matrix::colSums(counts > 0), n_UMIs = Matrix::colSums(counts)
  )
  safe_fraction <- function(selected) {
    if (!length(selected)) return(rep(0, ncol(counts)))
    Matrix::colSums(counts[selected, , drop = FALSE]) / pmax(metadata$n_UMIs, 1)
  }
  metadata$mito_frac <- safe_fraction(grep(mito_pattern, rownames(counts), value = TRUE))
  y_use <- intersect(chrY_genes, rownames(counts))
  metadata$chrY_frac <- if (length(y_use)) safe_fraction(y_use) else NA_real_

  velo_barcodes <- read.delim(file.path(velo_dir, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  spliced <- as(Matrix::readMM(file.path(velo_dir, "spliced.mtx.gz")), "CsparseMatrix")
  unspliced <- as(Matrix::readMM(file.path(velo_dir, "unspliced.mtx.gz")), "CsparseMatrix")
  colnames(spliced) <- velo_barcodes$V1
  colnames(unspliced) <- velo_barcodes$V1
  nuclear <- Matrix::colSums(unspliced) / (Matrix::colSums(spliced) + Matrix::colSums(unspliced))
  nuclear <- nuclear[is.finite(nuclear)]
  metadata$nuclear_frac <- NA_real_
  common <- intersect(rownames(metadata), names(nuclear))
  metadata[common, "nuclear_frac"] <- nuclear[common]

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = metadata)
  sce <- celda::decontX(sce)
  metadata$ambient_frac_decontx <- SummarizedExperiment::colData(sce)$decontX_contamination

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  doublet <- doubletFinder(counts, Seurat::VariableFeatures(seu), artificial_fraction)[[1]]
  metadata$doublet_score <- NA_real_
  metadata[names(doublet), "doublet_score"] <- doublet

  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  s_genes <- intersect(to_mouse_symbol(Seurat::cc.genes.updated.2019$s.genes), rownames(counts))
  g2m_genes <- intersect(to_mouse_symbol(Seurat::cc.genes.updated.2019$g2m.genes), rownames(counts))
  seu <- Seurat::CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes)
  metadata$phase <- seu@meta.data$Phase
  metadata$s_score <- seu@meta.data$S.Score
  metadata$g2m_score <- seu@meta.data$G2M.Score
  metadata$hbb_score <- safe_fraction(intersect(hbb_genes_mouse, rownames(counts)))
  metadata$is_HQ <- metadata$n_genes >= min_genes_hq

  Matrix::writeMM(counts, file.path(out_dir, "counts.mtx"))
  system2("gzip", c("-f", file.path(out_dir, "counts.mtx")))
  gz <- gzfile(file.path(out_dir, "features.tsv.gz"), "wt")
  write.table(gene_data$features, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(gz)
  gz <- gzfile(file.path(out_dir, "metadata.tsv.gz"), "wt")
  write.table(metadata, gz, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  close(gz)

  foreground <- metadata[metadata$is_HQ, , drop = FALSE]
  background <- metadata[!metadata$is_HQ, , drop = FALSE]
  plot <- ggplot2::ggplot(background, ggplot2::aes(x = nuclear_frac, y = n_UMIs)) +
    ggplot2::geom_point(size = .1, alpha = .2) +
    ggplot2::geom_point(data = foreground, color = "red", size = .1) +
    ggplot2::scale_y_log10() + ggplot2::labs(x = "Nuclear fraction", y = "Number of UMIs", title = sample_id) +
    ggplot2::theme_linedraw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))
  ggplot2::ggsave(file.path(out_dir, "qc_diagnosis.png"), plot, width = 4, height = 4, units = "in", dpi = 300)
}
