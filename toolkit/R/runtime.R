.scrna_runtime_started_at <- Sys.time()
.scrna_sha256 <- function(path) {
  if (is.null(path) || !file.exists(path) || !requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  if (dir.exists(path)) {
    members <- .scrna_directory_members(path)
    return(digest::digest(members, algo = "sha256"))
  }
  digest::digest(path, algo = "sha256", file = TRUE)
}

.scrna_directory_members <- function(path) {
  files <- sort(list.files(path, recursive = TRUE, full.names = FALSE, all.files = TRUE))
  files <- files[!file.info(file.path(path, files))$isdir]
  lapply(files, function(member) list(path = member,
    bytes = unname(file.info(file.path(path, member))$size),
    sha256 = .scrna_sha256(file.path(path, member))))
}

read_skill_config <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required")
  config <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  config_path <- normalizePath(path, mustWork = TRUE)
  input_config <- config[["input"]] %||% list()
  # Exact lookup prevents differential_table from partially matching the
  # differential_tables array and passing a list to file.exists().
  input <- input_config[["object"]] %||% input_config[["counts_table"]] %||% input_config[["path"]] %||%
    input_config[["differential_table"]]
  attr(config, "config_path") <- config_path
  attr(config, "input_record") <- list(
    path = input,
    bytes = if (!is.null(input) && file.exists(input)) unname(file.info(input)$size) else NA_real_,
    sha256 = .scrna_sha256(input)
  )
  if (!is.null(input) && dir.exists(input)) {
    attr(config, "input_record")$type <- "directory"
    attr(config, "input_record")$members <- .scrna_directory_members(input)
    attr(config, "input_record")$bytes <- sum(vapply(attr(config, "input_record")$members, function(x) x$bytes, numeric(1)))
  }
  if (is.null(input) && length(input_config[["differential_tables"]])) {
    tables <- lapply(input_config[["differential_tables"]], function(spec) {
      path <- spec[["path"]]
      list(path = path, bytes = unname(file.info(path)$size), sha256 = .scrna_sha256(path))
    })
    attr(config, "input_record") <- list(type = "tables", tables = tables,
      bytes = sum(vapply(tables, function(x) x$bytes, numeric(1))),
      sha256 = if (requireNamespace("digest", quietly = TRUE)) digest::digest(tables, algo = "sha256") else NA_character_)
  }
  config
}

cfg_get <- function(x, keys, default = NULL, required = FALSE) {
  value <- x
  for (key in strsplit(keys, "\\.", fixed = FALSE)[[1]]) {
    if (!is.list(value) || is.null(value[[key]])) {
      if (required) stop("Missing config field: ", keys)
      return(default)
    }
    value <- value[[key]]
  }
  value
}

load_scrna_object <- function(path, format = "auto", sample_id = NULL) {
  path <- normalizePath(path, mustWork = TRUE)
  if (format == "auto") {
    if (dir.exists(path)) format <- "10x_dir" else format <- sub("^.*\\.", "", tolower(path))
  }
  if (format %in% c("qs", "seurat_qs")) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required")
    return(qs::qread(path))
  }
  if (format %in% c("rds", "seurat_rds")) return(readRDS(path))
  if (format == "10x_dir") {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")
    uncompressed <- file.path(path, c("matrix.mtx", "barcodes.tsv", "features.tsv"))
    if (all(file.exists(uncompressed)) && !file.exists(file.path(path, "barcodes.tsv.gz"))) {
      # Modern Seurat Read10X expects gzipped v3 files; STARsolo and unpacked
      # bundles can contain the same three files without compression.
      counts <- methods::as(Matrix::readMM(uncompressed[1]), "CsparseMatrix")
      barcodes <- utils::read.delim(uncompressed[2], header = FALSE, colClasses = "character")
      features <- utils::read.delim(uncompressed[3], header = FALSE, colClasses = "character")
      if (nrow(features) != nrow(counts) || nrow(barcodes) != ncol(counts)) stop("10x feature/barcode dimensions do not match matrix")
      symbols <- features[[if (ncol(features) >= 2L) 2L else 1L]]
      rownames(counts) <- make.unique(symbols)
      colnames(counts) <- barcodes[[1]]
      if (ncol(features) >= 3L) {
        rna <- features[[3]] == "Gene Expression"
        if (!any(rna)) stop("10x input has no Gene Expression features")
        counts <- counts[rna, , drop = FALSE]
      }
    } else {
      counts <- Seurat::Read10X(path)
      if (is.list(counts)) {
        if (is.null(counts[["Gene Expression"]])) stop("10x input has no Gene Expression matrix")
        counts <- counts[["Gene Expression"]]
      }
    }
    obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_id %||% "scrna")
    if (!is.null(sample_id)) obj$sample_id <- sample_id
    return(obj)
  }
  if (format %in% c("h5", "10x_h5")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")
    counts <- Seurat::Read10X_h5(path)
    obj <- Seurat::CreateSeuratObject(counts = counts, project = sample_id %||% "scrna")
    if (!is.null(sample_id)) obj$sample_id <- sample_id
    return(obj)
  }
  stop("Unsupported input format: ", format)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

assert_metadata <- function(obj, columns) {
  missing <- setdiff(columns, colnames(obj[[]]))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

prepare_output <- function(config) {
  out <- cfg_get(config, "output_dir", required = TRUE)
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  normalizePath(out, mustWork = TRUE)
}

save_scrna_object <- function(obj, path) {
  if (grepl("\\.qs$", path, ignore.case = TRUE)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required")
    qs::qsave(obj, path)
  } else saveRDS(obj, path)
}

get_raw_counts <- function(obj, assay = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Package 'SeuratObject' is required")
  assay <- assay %||% Seurat::DefaultAssay(obj)
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    selected <- obj[[assay]]
    layers <- SeuratObject::Layers(selected, search = "^counts($|[.])")
    if (!length(layers)) stop("No raw counts layer in assay: ", assay)
    # 只在内存中合并，保证导出的矩阵覆盖全部对象细胞。
    if (length(layers) > 1L) selected <- SeuratObject::JoinLayers(selected, layers = "counts", new = "counts")
    counts <- SeuratObject::LayerData(selected, layer = if (length(layers) > 1L) "counts" else layers[[1]])
  } else {
    counts <- SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
  cells <- colnames(obj)
  if (!nrow(counts) || !ncol(counts)) stop("Raw counts matrix is empty: ", assay)
  if (anyDuplicated(colnames(counts)) || !setequal(cells, colnames(counts))) {
    stop("Raw counts cells do not match all object cells in assay: ", assay)
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
    stop("Raw counts must be finite non-negative integers: ", assay)
  }
  counts[, cells, drop = FALSE]
}

technical_path <- function(out, name) {
  directory <- file.path(out, "_provenance")
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  file.path(directory, name)
}

write_run_manifest <- function(config, skill, out, artifacts, notes = character(), exit_status = 0L) {
  figure_records <- file.path(out, "_provenance", c("figure_status.tsv", "figure_colors.tsv"))
  artifacts <- c(artifacts, figure_records[file.exists(figure_records)])
  # A rerun must not include an older copy of its own manifest as an artifact.
  artifacts <- artifacts[!grepl("^run_manifest.*\\.json$", basename(artifacts))]
  artifact_records <- lapply(unique(artifacts), function(path) {
    info <- file.info(path)
    # The supervising runner appends its exit record after this function returns.
    mutable <- grepl("\\.log$", path)
    list(path = path, bytes = if (isTRUE(info$isdir) || mutable) NA_real_ else unname(info$size),
         sha256 = if (mutable) NA_character_ else .scrna_sha256(path), mutable = mutable)
  })
  config_path <- attr(config, "config_path")
  input_record <- attr(config, "input_record")
  finished_at <- Sys.time()
  manifest <- list(
    schema_version = 2L,
    run_id = Sys.getenv("SCRNA_RUN_ID", unset = ""),
    skill = skill,
    project_id = cfg_get(config, "project.id", required = TRUE),
    started_at = format(.scrna_runtime_started_at, tz = "UTC", usetz = TRUE),
    finished_at = format(finished_at, tz = "UTC", usetz = TRUE),
    duration_seconds = as.numeric(difftime(finished_at, .scrna_runtime_started_at, units = "secs")),
    command = commandArgs(FALSE),
    config = list(path = config_path, sha256 = .scrna_sha256(config_path)),
    input = input_record,
    output_dir = out,
    artifacts = artifact_records,
    random_seed = attr(config, "resolved_random_seed") %||% cfg_get(config, "preprocessing.seed", cfg_get(config, "random_seed")),
    exit_status = as.integer(exit_status),
    status = if (exit_status == 0L) "completed" else "failed",
    notes = unname(as.list(notes))
  )
  action <- cfg_get(config, "workflow.action")
  manifest_name <- if (skill == "06-scrna-preprocess-and-cluster" && identical(action, "finalize_resolution")) {
    "run_manifest_finalize.json"
  } else if (skill == "06-scrna-preprocess-and-cluster") {
    "run_manifest_preprocess.json"
  } else {
    "run_manifest.json"
  }
  jsonlite::write_json(manifest, technical_path(out, manifest_name), auto_unbox = TRUE, pretty = TRUE)
}
