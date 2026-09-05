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
provenance_dir <- file.path(out_root, "_provenance")
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
species <- tolower(getv("parameters.species", "mouse"))
if (!species %in% c("mouse", "human")) stop("parameters.species must be mouse or human")
min_genes_hq <- as.integer(getv("parameters.min_genes_hq", 500))
mito_pattern <- getv("parameters.mito_pattern", if (species == "mouse") "^mt-" else "^MT-")
artificial_fraction <- as.numeric(getv("parameters.doublet_artificial_fraction", 0.2))
base_seed <- as.integer(getv("parameters.seed", 1))
if (is.na(base_seed)) stop("parameters.seed must be an integer")
workers <- as.integer(getv("parallel.workers", 1))
if (is.na(workers) || workers < 1L) stop("parallel.workers must be a positive integer")
ambient_method <- tolower(getv("ambient_rna.method", "decontx"))
if (!ambient_method %in% c("decontx", "skip")) stop("ambient_rna.method must be decontx or skip")
decontx_cluster_col <- getv("ambient_rna.cluster_column")
if (!is.null(decontx_cluster_col) && (!is.character(decontx_cluster_col) || length(decontx_cluster_col) != 1L || !nzchar(decontx_cluster_col))) {
  stop("ambient_rna.cluster_column must be null or a non-empty metadata column name")
}
set.seed(base_seed)
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
  attributes <- vapply(genes, `[[`, character(1), 9)
  has_name <- grepl('gene_name "[^"]+"', attributes)
  unique(sub('.*gene_name "([^"]+)".*', "\\1", attributes[has_name]))
}
chrY_genes <- read_chrY(getv("input.gtf_file"))

get_counts <- function(obj) {
  requested <- getv("parameters.assay")
  assay <- if (!is.null(requested)) requested else if ("RNA" %in% Seurat::Assays(obj)) "RNA" else Seurat::DefaultAssay(obj)
  if (!assay %in% Seurat::Assays(obj)) stop("Configured assay not found in Seurat object: ", assay)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    count_layers <- SeuratObject::Layers(obj[[assay]], search = "^counts($|\\.)")
    if (!length(count_layers)) stop("The selected assay has no counts layer: ", assay)
    if (length(count_layers) > 1L) {
      obj[[assay]] <- SeuratObject::JoinLayers(
        obj[[assay]], layers = "counts", new = "counts"
      )
    }
    counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  } else {
    counts <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  }
  if (!nrow(counts) || !ncol(counts)) stop("The selected assay has no usable counts matrix: ", assay)
  object_cells <- colnames(obj)
  missing_cells <- setdiff(object_cells, colnames(counts))
  extra_cells <- setdiff(colnames(counts), object_cells)
  if (length(missing_cells) || length(extra_cells)) {
    stop(
      "Counts cells do not match the Seurat object for assay ", assay,
      ": ", length(missing_cells), " missing and ", length(extra_cells), " extra"
    )
  }
  counts <- counts[, object_cells, drop = FALSE]
  stored_values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (length(stored_values) && (any(!is.finite(stored_values)) || any(stored_values < 0))) stop("Counts must be finite and non-negative in assay: ", assay)
  counts
}
read_star <- function(directory) {
  counts <- as(Matrix::readMM(file.path(directory, "matrix.mtx.gz")), "CsparseMatrix")
  features <- read.delim(file.path(directory, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(file.path(directory, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- make.unique(features$V2)
  colnames(counts) <- barcodes$V1
  list(counts = counts, features = features)
}

calculate_one <- function(counts, metadata, sample_id, velocity_dir = NULL, sample_seed = base_seed, decontx_z = NULL) {
  if (!identical(rownames(metadata), colnames(counts))) {
    if (!setequal(rownames(metadata), colnames(counts))) stop("Metadata cells do not match counts for sample: ", sample_id)
    metadata <- metadata[colnames(counts), , drop = FALSE]
  }
  set.seed(sample_seed)
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
  if (ambient_method == "skip") {
    statuses <- rbind(statuses, status_row(sample_id, "ambient_frac_decontx", "skipped", "disabled_by_config"))
  } else if (requireNamespace("celda", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    ambient_error <- tryCatch({
      if (!is.null(decontx_z)) {
        if (length(decontx_z) != ncol(counts)) stop("DecontX cluster labels do not match the sample cell count")
        if (anyNA(decontx_z) || any(!nzchar(as.character(decontx_z)))) stop("DecontX cluster labels contain missing or empty values")
        if (length(unique(decontx_z)) < 2L) stop("DecontX cluster labels must contain at least two broad cell populations")
      }
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = metadata)
      sce <- celda::decontX(sce, z = decontx_z, seed = sample_seed)
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
      score <- sc06SyntheticDoubletScore(counts, Seurat::VariableFeatures(seu), artificial_fraction)[[1]]
      metadata[names(score), "doublet_score"] <- score
      NULL
    }, error = function(e) conditionMessage(e))
    statuses <- rbind(statuses, status_row(sample_id, "doublet_score", if (is.null(doublet_error)) "computed" else "skipped",
                                          if (is.null(doublet_error)) "sc06 synthetic-doublet nearest-neighbour score; not canonical DoubletFinder classification" else doublet_error))
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

run_sample_jobs <- function(items, FUN) {
  if (!length(items)) return(list())
  n_workers <- min(workers, length(items))
  indices <- seq_along(items)
  if (n_workers == 1L || .Platform$OS.type == "windows") {
    if (workers > 1L && .Platform$OS.type == "windows") {
      warning("parallel.workers > 1 is not supported by this executor on Windows; using one worker")
    }
    return(lapply(indices, function(i) FUN(items[[i]], i)))
  }
  message("Running ", length(items), " samples with ", n_workers, " workers")
  results <- parallel::mclapply(
    indices,
    function(i) FUN(items[[i]], i),
    mc.cores = n_workers,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("Sample workers failed: ", paste(vapply(results[failed], as.character, character(1)), collapse = "; "))
  }
  results
}

if (input_type == "starsolo") {
  root <- normalizePath(getv("input.starsolo_dir"), mustWork = TRUE)
  process_star_sample <- function(sample, sample_index) {
    sample_id <- sample$sample_id; batch_id <- sample$batch_id
    gene_dir <- file.path(root, sample_id, "Solo.out/Gene/filtered")
    velocity_dir <- file.path(root, sample_id, "Solo.out/Velocyto/filtered")
    raw <- read_star(gene_dir)
    metadata <- data.frame(
      row.names = colnames(raw$counts),
      sample_id = rep(sample_id, ncol(raw$counts)),
      batch_id = rep(batch_id, ncol(raw$counts))
    )
    result <- calculate_one(raw$counts, metadata, sample_id, velocity_dir, base_seed + sample_index - 1L)
    out <- file.path(out_root, sample_id); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    Matrix::writeMM(raw$counts, file.path(out, "counts.mtx")); system2("gzip", c("-f", file.path(out, "counts.mtx")))
    gz <- gzfile(file.path(out, "features.tsv.gz"), "wt"); write.table(raw$features, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(gz)
    gz <- gzfile(file.path(out, "metadata.tsv.gz"), "wt"); write.table(result$metadata, gz, sep = "\t", quote = FALSE); close(gz)
    plot_qc(result$metadata, file.path(out, "qc_diagnosis.png"), sample_id)
    list(sample_id = sample_id, status = result$status)
  }
  sample_results <- run_sample_jobs(getv("samples"), process_star_sample)
  all_status <- lapply(sample_results, `[[`, "status")
  write.table(do.call(rbind, all_status), file.path(provenance_dir, "metric_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
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
  if (!is.null(decontx_cluster_col) && !decontx_cluster_col %in% colnames(original_meta)) {
    stop("Configured ambient_rna.cluster_column not found in Seurat metadata: ", decontx_cluster_col)
  }
  cluster_values <- if (!is.null(decontx_cluster_col)) as.character(original_meta[[decontx_cluster_col]]) else NULL
  names(sample_values) <- names(batch_values) <- rownames(original_meta)
  if (!is.null(cluster_values)) names(cluster_values) <- rownames(original_meta)
  combined_meta <- original_meta; all_status <- list()
  sample_ids <- unique(sample_values)
  process_seurat_sample <- function(sample_id, sample_index) {
    cells <- names(sample_values)[sample_values == sample_id]
    meta <- data.frame(row.names = cells, sample_id = sample_id, batch_id = batch_values[cells])
    result <- calculate_one(
      counts[, cells, drop = FALSE], meta, sample_id,
      sample_seed = base_seed + sample_index - 1L,
      decontx_z = if (!is.null(cluster_values)) cluster_values[cells] else NULL
    )
    list(sample_id = sample_id, cells = cells, result = result)
  }
  sample_results <- run_sample_jobs(sample_ids, process_seurat_sample)
  for (sample_result in sample_results) {
    sample_id <- sample_result$sample_id
    cells <- sample_result$cells
    result <- sample_result$result
    for (column in setdiff(colnames(result$metadata), c("sample_id", "batch_id"))) combined_meta[cells, column] <- result$metadata[cells, column]
    all_status[[sample_id]] <- result$status
  }
  qc_columns <- c("n_genes", "n_UMIs", "mito_frac", "chrY_frac", "nuclear_frac",
                  "ambient_frac_decontx", "doublet_score", "phase", "s_score", "g2m_score",
                  "hbb_score", "is_HQ")
  for (column in intersect(qc_columns, colnames(combined_meta))) {
    obj[[column]] <- combined_meta[colnames(obj), column]
  }
  saveRDS(obj, file.path(out_root, "qc_metrics_object.rds"))
  gz <- gzfile(file.path(out_root, "metadata.tsv.gz"), "wt"); write.table(combined_meta, gz, sep = "\t", quote = FALSE); close(gz)
  write.table(do.call(rbind, all_status), file.path(provenance_dir, "metric_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_qc(combined_meta, file.path(out_root, "qc_diagnosis.png"), basename(object_path))
}
