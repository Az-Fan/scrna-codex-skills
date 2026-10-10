.figure_script <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if (is.null(.figure_script)) .figure_script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
if (length(.figure_script) && !is.na(.figure_script)) {
  source(file.path(dirname(normalizePath(.figure_script)), "figure_output.R"))
} else if (file.exists("toolkit/R/figure_output.R")) source("toolkit/R/figure_output.R")
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

scrna_expand_path <- function(path) {
  matches <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}|\\$[A-Za-z_][A-Za-z0-9_]*", path)[[1]]
  widths <- attr(matches, "match.length")
  for (i in rev(seq_along(matches))) {
    if (matches[i] < 0L) next
    token <- substr(path, matches[i], matches[i] + widths[i] - 1L)
    name <- gsub("[${}]", "", token)
    path <- paste0(substr(path, 1L, matches[i] - 1L), Sys.getenv(name, unset = token), substring(path, matches[i] + widths[i]))
  }
  path.expand(path)
}

scrna_environment_root <- function(config) {
  root <- cfg_get(config, "runtime.pixi_root")
  if (is.null(root) || !nzchar(root)) root <- Sys.getenv("SCRNA_PIXI_ROOT")
  if (!nzchar(root)) root <- "~/projects/scrna_envs"
  normalizePath(scrna_expand_path(root), mustWork = FALSE)
}

scrna_python_prefix <- function(config, method = "integration") {
  if (method == "integration") {
    prefix <- cfg_get(config, "benchmark.python_argv_prefix")
    if (!is.null(prefix)) {
      prefix <- unlist(prefix, use.names = FALSE)
      prefix[[1]] <- scrna_expand_path(prefix[[1]])
      return(prefix)
    }
    return(file.path(scrna_environment_root(config), "03-integration", ".pixi", "envs", "scvi", "bin", "python"))
  }
  configured <- cfg_get(config, "runtime.sccoda_python")
  if (!is.null(configured) && nzchar(configured)) return(scrna_expand_path(configured))
  file.path(scrna_environment_root(config), "07-cell-abundance", ".pixi", "envs", "sccoda", "bin", "python")
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
  selected <- obj[[assay]]
  if (inherits(selected, "Assay5")) {
    layers <- SeuratObject::Layers(selected, search = "^counts($|[.])")
    if (!length(layers)) stop("No raw counts layer in assay: ", assay)
    # 只在内存中合并，保证导出的矩阵覆盖全部对象细胞。
    if (length(layers) > 1L) selected <- SeuratObject::JoinLayers(selected, layers = "counts", new = "counts")
    counts <- SeuratObject::LayerData(selected, layer = if (length(layers) > 1L) "counts" else layers[[1]])
  } else {
    counts <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      SeuratObject::GetAssayData(selected, layer = "counts")
    } else SeuratObject::GetAssayData(selected, slot = "counts")
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

assert_pseudobulk_counts_source <- function(obj, assay, config) {
  source <- cfg_get(config, "analysis.counts_source", list())
  if (!is.list(source)) stop("analysis.counts_source must be an object")
  kind <- source$kind %||% ""
  if (length(kind) != 1L || is.na(kind) || !kind %in% c("raw_umi", "raw_read")) stop("Formal pseudobulk requires analysis.counts_source.kind = raw_umi or raw_read")
  if (inherits(obj[[assay]], "SCTAssay") || grepl("(^|[._-])(integrated|sct|harmony|corrected|scaled|normalized)([._-]|$)", assay, ignore.case = TRUE)) {
    stop("Formal pseudobulk cannot use an integrated, SCT, corrected or normalized assay: ", assay)
  }
  counts <- get_raw_counts(obj, assay)
  reference_path <- source$reference_object
  if (!is.null(reference_path) && (length(reference_path) != 1L || is.na(reference_path) || !is.character(reference_path))) stop("counts_source.reference_object must be a path string")
  evidence <- "user_declared; original source is not independently authenticated"
  reference_hash <- NA_character_
  if (!is.null(reference_path) && nzchar(reference_path)) {
    reference_path <- normalizePath(scrna_expand_path(reference_path), mustWork = TRUE)
    reference <- load_scrna_object(reference_path)
    reference_assay <- source$reference_assay %||% "RNA"
    if (inherits(reference[[reference_assay]], "SCTAssay") || grepl("(^|[._-])(integrated|sct|harmony|corrected|scaled|normalized)([._-]|$)", reference_assay, ignore.case = TRUE)) stop("Raw-count reference must use an uncorrected assay")
    reference_counts <- get_raw_counts(reference, reference_assay)
    if (!all(rownames(counts) %in% rownames(reference_counts)) || !all(colnames(counts) %in% colnames(reference_counts))) stop("Raw-count reference does not cover every selected gene and cell")
    if (any(counts != reference_counts[rownames(counts), colnames(counts), drop = FALSE])) stop("Selected counts differ from the supplied raw-count reference")
    reference_hash <- .scrna_sha256(reference_path)
    evidence <- "matched_supplied_reference; reference origin is user_declared"
  }
  data.frame(assay = assay, assay_class = class(obj[[assay]])[[1]], counts_kind = kind,
    evidence = evidence, reference_object = reference_path %||% "", reference_sha256 = reference_hash,
    n_genes = nrow(counts), n_cells = ncol(counts), stringsAsFactors = FALSE)
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
