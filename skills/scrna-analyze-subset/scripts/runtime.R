read_skill_config <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
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
  manifest <- list(
    schema_version = 1L,
    skill = skill,
    project_id = cfg_get(config, "project.id", required = TRUE),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input = cfg_get(config, "input.object", cfg_get(config, "input.path")),
    output_dir = out,
    artifacts = unname(as.list(artifacts)),
    notes = unname(as.list(notes))
  )
  jsonlite::write_json(manifest, file.path(out, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
}
