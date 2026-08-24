args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript calculate_metrics.R config.json")
for (pkg in c("jsonlite", "Matrix", "Seurat", "SeuratObject", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Core package missing from existing pixi environment: ", pkg)
}
config <- jsonlite::read_json(args[[1]], simplifyVector = FALSE)
getv <- function(path, default = NULL) {
  value <- config
  for (key in strsplit(path, "\\.")[[1]]) {
    if (is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

input_type <- getv("input.type", "auto")
if (input_type == "auto") input_type <- if (!is.null(getv("input.object"))) "seurat" else "starsolo"
out_root <- getv("output_dir")
if (is.null(out_root) || !nzchar(out_root)) stop("output_dir is required")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
species <- tolower(getv("parameters.species", "mouse"))
if (!species %in% c("mouse", "human")) stop("parameters.species must be mouse or human")
min_genes_hq <- as.integer(getv("parameters.min_genes_hq", 500))
mito_pattern <- getv("parameters.mito_pattern", if (species == "mouse") "^mt-" else "^MT-")
artificial_fraction <- as.numeric(getv("parameters.doublet_artificial_fraction", 0.2))
set.seed(as.integer(getv("parameters.seed", 1)))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "doublet_finder.R"))

hbb_genes <- if (species == "mouse") {
  c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh1", "Hbb-y", "Hbb-bh2", "Hbb-b1", "Hbb-b2", "Hbm")
} else c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
cycle_symbols <- function(x) {
  if (species == "human") return(x)
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}
status_row <- function(sample_id, metric, status, reason = "") {
  data.frame(sample_id = sample_id, metric = metric, status = status, reason = reason, stringsAsFactors = FALSE)
}
read_chrY <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines)]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1)) >= 9]
  genes <- fields[vapply(fields, function(x) x[[3]] == "gene" && x[[1]] %in% c("Y", "chrY"), logical(1))]
  unique(vapply(genes, function(x) sub('.*gene_name "([^"]+)".*', "\\1", x[[9]]), character(1)))
}
chrY_genes <- read_chrY(getv("input.gtf_file"))

get_counts <- function(obj) {
  assay <- Seurat::DefaultAssay(obj)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(SeuratObject::LayerData(obj, assay = assay, layer = "counts"))
  }
  Seurat::GetAssayData(obj, assay = assay, slot = "counts")
}
read_star <- function(directory) {
  counts <- as(Matrix::readMM(file.path(directory, "matrix.mtx.gz")), "CsparseMatrix")
  features <- read.delim(file.path(directory, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(file.path(directory, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- make.unique(features$V2)
  colnames(counts) <- barcodes$V1
  list(counts = counts, features = features)
}

calculate_one <- function(counts, metadata, sample_id, velocity_dir = NULL) {
  metadata$n_genes <- Matrix::colSums(counts > 0)
  metadata$n_UMIs <- Matrix::colSums(counts)
  statuses <- rbind(status_row(sample_id, "n_genes", "computed"), status_row(sample_id, "n_UMIs", "computed"))
  fraction <- function(genes) {
    if (!length(genes)) return(rep(NA_real_, ncol(counts)))
    Matrix::colSums(counts[genes, , drop = FALSE]) / pmax(metadata$n_UMIs, 1)
  }

  mito <- grep(mito_pattern, rownames(counts), value = TRUE)
  metadata$mito_frac <- fraction(mito)
  statuses <- rbind(statuses, status_row(sample_id, "mito_frac", if (length(mito)) "computed" else "skipped", if (length(mito)) "" else "no matching mitochondrial genes"))
  y_use <- intersect(chrY_genes %||% character(), rownames(counts))
  metadata$chrY_frac <- fraction(y_use)
  statuses <- rbind(statuses, status_row(sample_id, "chrY_frac", if (length(y_use)) "computed" else "skipped", if (is.null(chrY_genes)) "GTF not provided" else if (!length(y_use)) "no matching chrY genes" else ""))

  metadata$nuclear_frac <- NA_real_
  velocity_files <- if (!is.null(velocity_dir)) file.path(velocity_dir, c("spliced.mtx.gz", "unspliced.mtx.gz", "barcodes.tsv.gz")) else character()
  if (length(velocity_files) && all(file.exists(velocity_files))) {
    barcodes <- read.delim(velocity_files[[3]], header = FALSE, stringsAsFactors = FALSE)$V1
    spliced <- as(Matrix::readMM(velocity_files[[1]]), "CsparseMatrix")
    unspliced <- as(Matrix::readMM(velocity_files[[2]]), "CsparseMatrix")
    colnames(spliced) <- colnames(unspliced) <- barcodes
    nuclear <- Matrix::colSums(unspliced) / (Matrix::colSums(spliced) + Matrix::colSums(unspliced))
    nuclear <- nuclear[is.finite(nuclear)]
    common <- intersect(rownames(metadata), names(nuclear))
    metadata[common, "nuclear_frac"] <- nuclear[common]
    statuses <- rbind(statuses, status_row(sample_id, "nuclear_frac", "computed"))
  } else statuses <- rbind(statuses, status_row(sample_id, "nuclear_frac", "skipped", "Velocyto matrices not provided"))

  metadata$ambient_frac_decontx <- NA_real_
  if (requireNamespace("celda", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    ambient_error <- tryCatch({
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = metadata)
      sce <- celda::decontX(sce)
      metadata$ambient_frac_decontx <- SummarizedExperiment::colData(sce)$decontX_contamination
      NULL
    }, error = function(e) conditionMessage(e))
    statuses <- rbind(statuses, status_row(sample_id, "ambient_frac_decontx", if (is.null(ambient_error)) "computed" else "skipped", ambient_error %||% ""))
  } else statuses <- rbind(statuses, status_row(sample_id, "ambient_frac_decontx", "skipped", "celda or SingleCellExperiment unavailable"))

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  metadata$doublet_score <- NA_real_
  if (requireNamespace("RANN", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)) {
    doublet_error <- tryCatch({
      seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      score <- doubletFinder(counts, Seurat::VariableFeatures(seu), artificial_fraction)[[1]]
      metadata[names(score), "doublet_score"] <- score
      NULL
    }, error = function(e) conditionMessage(e))
    statuses <- rbind(statuses, status_row(sample_id, "doublet_score", if (is.null(doublet_error)) "computed" else "skipped", doublet_error %||% ""))
  } else statuses <- rbind(statuses, status_row(sample_id, "doublet_score", "skipped", "RANN or S4Vectors unavailable"))

  metadata$phase <- NA_character_; metadata$s_score <- NA_real_; metadata$g2m_score <- NA_real_
  s_genes <- intersect(cycle_symbols(Seurat::cc.genes.updated.2019$s.genes), rownames(counts))
  g2m_genes <- intersect(cycle_symbols(Seurat::cc.genes.updated.2019$g2m.genes), rownames(counts))
  if (length(s_genes) && length(g2m_genes)) {
    cycle_error <- tryCatch({
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      seu <- Seurat::CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes)
      metadata$phase <- seu@meta.data$Phase; metadata$s_score <- seu@meta.data$S.Score; metadata$g2m_score <- seu@meta.data$G2M.Score
      NULL
    }, error = function(e) conditionMessage(e))
    statuses <- rbind(statuses, status_row(sample_id, "cell_cycle", if (is.null(cycle_error)) "computed" else "skipped", cycle_error %||% ""))
  } else statuses <- rbind(statuses, status_row(sample_id, "cell_cycle", "skipped", "cell-cycle genes not found"))

  hbb_use <- intersect(hbb_genes, rownames(counts))
  metadata$hbb_score <- fraction(hbb_use)
  statuses <- rbind(statuses, status_row(sample_id, "hbb_score", if (length(hbb_use)) "computed" else "skipped", if (length(hbb_use)) "" else "hemoglobin genes not found"))
  metadata$is_HQ <- metadata$n_genes >= min_genes_hq
  list(metadata = metadata, status = statuses)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
plot_qc <- function(metadata, path, title = NULL) {
  x_metric <- if (any(is.finite(metadata$nuclear_frac))) "nuclear_frac" else "n_genes"
  plot <- ggplot2::ggplot(metadata, ggplot2::aes(x = .data[[x_metric]], y = n_UMIs, color = is_HQ)) +
    ggplot2::geom_point(size = .15, alpha = .5) + ggplot2::scale_y_log10() +
    ggplot2::labs(x = x_metric, y = "Number of UMIs", title = title) + ggplot2::theme_linedraw()
  ggplot2::ggsave(path, plot, width = 5, height = 4, units = "in", dpi = 300)
}

if (input_type == "starsolo") {
  root <- normalizePath(getv("input.starsolo_dir"), mustWork = TRUE)
  for (sample in getv("samples")) {
    sample_id <- sample$sample_id; batch_id <- sample$batch_id
    gene_dir <- file.path(root, sample_id, "Solo.out/Gene/filtered")
    velocity_dir <- file.path(root, sample_id, "Solo.out/Velocyto/filtered")
    raw <- read_star(gene_dir)
    metadata <- data.frame(row.names = colnames(raw$counts), sample_id = sample_id, batch_id = batch_id)
    result <- calculate_one(raw$counts, metadata, sample_id, velocity_dir)
    out <- file.path(out_root, sample_id); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    Matrix::writeMM(raw$counts, file.path(out, "counts.mtx")); system2("gzip", c("-f", file.path(out, "counts.mtx")))
    gz <- gzfile(file.path(out, "features.tsv.gz"), "wt"); write.table(raw$features, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(gz)
    gz <- gzfile(file.path(out, "metadata.tsv.gz"), "wt"); write.table(result$metadata, gz, sep = "\t", quote = FALSE); close(gz)
    write.table(result$status, file.path(out, "metric_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    plot_qc(result$metadata, file.path(out, "qc_diagnosis.png"), sample_id)
  }
} else {
  object_path <- normalizePath(getv("input.object"), mustWork = TRUE)
  ext <- tolower(tools::file_ext(object_path))
  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) stop("QS input requires qs in the existing pixi environment")
    obj <- qs::qread(object_path)
  } else obj <- readRDS(object_path)
  if (!inherits(obj, "Seurat")) stop("input.object is not a Seurat object")
  counts <- get_counts(obj); original_meta <- obj[[]]
  sample_col <- getv("metadata.sample"); batch_col <- getv("metadata.batch")
  sample_values <- if (!is.null(sample_col) && sample_col %in% colnames(original_meta)) as.character(original_meta[[sample_col]]) else rep("all_cells", nrow(original_meta))
  batch_values <- if (!is.null(batch_col) && batch_col %in% colnames(original_meta)) as.character(original_meta[[batch_col]]) else rep(NA_character_, nrow(original_meta))
  names(sample_values) <- names(batch_values) <- rownames(original_meta)
  combined_meta <- original_meta; all_status <- list()
  for (sample_id in unique(sample_values)) {
    cells <- names(sample_values)[sample_values == sample_id]
    meta <- data.frame(row.names = cells, sample_id = sample_id, batch_id = batch_values[cells])
    result <- calculate_one(counts[, cells, drop = FALSE], meta, sample_id)
    for (column in setdiff(colnames(result$metadata), c("sample_id", "batch_id"))) combined_meta[cells, column] <- result$metadata[cells, column]
    all_status[[sample_id]] <- result$status
  }
  for (column in setdiff(colnames(combined_meta), colnames(obj[[]]))) obj[[column]] <- combined_meta[colnames(obj), column]
  saveRDS(obj, file.path(out_root, "qc_metrics_object.rds"))
  gz <- gzfile(file.path(out_root, "metadata.tsv.gz"), "wt"); write.table(combined_meta, gz, sep = "\t", quote = FALSE); close(gz)
  write.table(do.call(rbind, all_status), file.path(out_root, "metric_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_qc(combined_meta, file.path(out_root, "qc_diagnosis.png"), basename(object_path))
}
