.scrna_runtime_started_at <- Sys.time()
.scrna_sha256 <- function(path) {
  if (is.null(path) || !file.exists(path) || !requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(path, algo = "sha256", file = TRUE)
}

read_skill_config <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required")
  config <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  config_path <- normalizePath(path, mustWork = TRUE)
  input <- config$input$object %||% config$input$path
  attr(config, "config_path") <- config_path
  attr(config, "input_record") <- list(
    path = input,
    bytes = if (!is.null(input) && file.exists(input)) unname(file.info(input)$size) else NA_real_,
    sha256 = .scrna_sha256(input)
  )
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
    counts <- Seurat::Read10X(path)
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
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  } else {
    SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
}

write_run_manifest <- function(config, skill, out, artifacts, notes = character()) {
  artifact_records <- lapply(unique(artifacts), function(path) {
    info <- file.info(path)
    list(path = path, bytes = if (isTRUE(info$isdir)) NA_real_ else unname(info$size), sha256 = .scrna_sha256(path))
  })
  config_path <- attr(config, "config_path")
  input_record <- attr(config, "input_record")
  finished_at <- Sys.time()
  manifest <- list(
    schema_version = 2L,
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
    exit_status = 0L,
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
  jsonlite::write_json(manifest, file.path(out, manifest_name), auto_unbox = TRUE, pretty = TRUE)
}
