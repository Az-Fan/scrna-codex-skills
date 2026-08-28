args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript integration_benchmark.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
config_path <- normalizePath(args[[1]], mustWork = TRUE)
config <- read_skill_config(config_path)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required")

as_chr <- function(x) unlist(x %||% list(), use.names = FALSE)
safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
grid_rows <- function(grid) {
  if (is.null(grid) || !length(grid)) return(list(list()))
  values <- lapply(grid, function(x) unlist(x, use.names = FALSE))
  frame <- expand.grid(values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  lapply(seq_len(nrow(frame)), function(i) as.list(frame[i, , drop = FALSE]))
}
scenario_id <- function(method, params) {
  if (!length(params)) return(method)
  suffix <- paste(vapply(names(params), function(nm) paste0(nm, "_", params[[nm]]), character(1)), collapse = "__")
  safe_id(paste(method, suffix, sep = "__"))
}
add_reduction <- function(target, embeddings, name, assay, key) {
  embeddings <- embeddings[colnames(target), , drop = FALSE]
  colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))
  target[[name]] <- SeuratObject::CreateDimReducObject(embeddings = embeddings, assay = assay, key = key)
  target
}

obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE), "auto")
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
batch_cols <- as_chr(cfg_get(config, "metadata.batch_variables", required = TRUE))
label_cols <- as_chr(cfg_get(config, "metadata.biological_labels", list()))
condition_col <- cfg_get(config, "metadata.condition")
required_meta <- unique(c(sample_col, batch_cols, label_cols, condition_col))
required_meta <- required_meta[!is.na(required_meta) & nzchar(required_meta)]
assert_metadata(obj, required_meta)
for (column in batch_cols) if (length(unique(stats::na.omit(obj[[]][[column]]))) < 2L) stop("Batch variable has fewer than two levels: ", column)

assay <- cfg_get(config, "input.assay", if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj))
Seurat::DefaultAssay(obj) <- assay
dims <- as.integer(as_chr(cfg_get(config, "benchmark.dims", as.list(1:30))))
seed <- as.integer(cfg_get(config, "benchmark.seed", 1234)); set.seed(seed)
if (!"pca" %in% names(obj@reductions)) {
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = max(dims), seed.use = seed, verbose = FALSE)
}

methods <- cfg_get(config, "benchmark.methods", required = TRUE)
method_names <- vapply(methods, function(x) tolower(x$name), character(1))
if (!"none" %in% method_names) methods <- c(list(list(name = "none")), methods)
scenarios <- list()
for (method in methods) for (params in grid_rows(method$parameter_grid)) {
  name <- tolower(method$name)
  scenarios[[length(scenarios) + 1L]] <- list(name = name, id = scenario_id(name, params), params = params, spec = method)
}

out <- prepare_output(config)
exchange <- file.path(out, "exchange"); embeddings_dir <- file.path(exchange, "embeddings")
dir.create(embeddings_dir, recursive = TRUE, showWarnings = FALSE)
run_rows <- list(); embedding_rows <- list()
for (scenario in scenarios[vapply(scenarios, function(x) x$name %in% c("none", "harmony", "rpca", "precomputed"), logical(1))]) {
  status <- "completed"; note <- ""; reduction <- NA_character_; representation_type <- "embedding"
  tryCatch({
    if (scenario$name == "none") {
      reduction <- "pca"
    } else if (scenario$name == "harmony") {
      if (!requireNamespace("harmony", quietly = TRUE)) stop("missing dependency: harmony")
      reduction <- paste0("bench_", scenario$id)
      obj <- harmony::RunHarmony(obj, group.by.vars = batch_cols, reduction = "pca", dims.use = dims,
        theta = as.numeric(scenario$params$theta %||% 2), reduction.save = reduction, verbose = FALSE)
    } else if (scenario$name == "rpca") {
      split_key <- ".benchmark_batch"
      obj[[split_key]] <- interaction(obj[[]][batch_cols], drop = TRUE, lex.order = TRUE)
      parts <- Seurat::SplitObject(obj, split.by = split_key)
      parts <- lapply(parts, function(x) Seurat::NormalizeData(x, verbose = FALSE) |>
        Seurat::FindVariableFeatures(nfeatures = as.integer(scenario$params$nfeatures %||% 2000), verbose = FALSE))
      features <- Seurat::SelectIntegrationFeatures(parts, nfeatures = as.integer(scenario$params$nfeatures %||% 2000))
      parts <- lapply(parts, function(x) Seurat::ScaleData(x, features = features, verbose = FALSE) |>
        Seurat::RunPCA(features = features, npcs = max(dims), seed.use = seed, verbose = FALSE))
      anchors <- Seurat::FindIntegrationAnchors(parts, anchor.features = features, reduction = "rpca", dims = dims,
        k.anchor = as.integer(scenario$params$k_anchor %||% 5), verbose = FALSE)
      integrated <- Seurat::IntegrateData(anchors, dims = dims,
        k.weight = as.integer(scenario$params$k_weight %||% 50), verbose = FALSE)
      Seurat::DefaultAssay(integrated) <- "integrated"
      integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
      integrated <- Seurat::RunPCA(integrated, npcs = max(dims), reduction.name = "rpca", seed.use = seed, verbose = FALSE)
      reduction <- paste0("bench_", scenario$id)
      obj <- add_reduction(obj, Seurat::Embeddings(integrated, "rpca"), reduction, assay, paste0("RPCA", length(embedding_rows) + 1L, "_"))
      obj[[split_key]] <- NULL
    } else {
      reduction <- scenario$spec$reduction %||% stop("precomputed method requires reduction")
      if (!reduction %in% names(obj@reductions)) stop("missing precomputed reduction: ", reduction)
    }
    embedding_path <- file.path(embeddings_dir, paste0(scenario$id, ".tsv.gz"))
    emb <- as.data.frame(Seurat::Embeddings(obj, reduction)[, dims, drop = FALSE]); emb$cell_id <- rownames(emb)
    con <- gzfile(embedding_path, "wt"); utils::write.table(emb[, c("cell_id", setdiff(names(emb), "cell_id"))], con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
    embedding_rows[[length(embedding_rows) + 1L]] <- data.frame(scenario = scenario$id, method = scenario$name,
      parameters = jsonlite::toJSON(scenario$params, auto_unbox = TRUE),
      path = normalizePath(embedding_path, winslash = "/", mustWork = TRUE), representation_type = representation_type)
  }, error = function(e) {
    status <<- if (grepl("missing dependency", conditionMessage(e), fixed = TRUE)) "missing_dependency" else "failed"
    note <<- conditionMessage(e)
  })
  run_rows[[length(run_rows) + 1L]] <- data.frame(scenario = scenario$id, method = scenario$name,
    parameters = jsonlite::toJSON(scenario$params, auto_unbox = TRUE), representation = reduction, status = status, notes = note)
}

metadata <- obj[[]]; metadata$cell_id <- rownames(metadata)
write_tsv(metadata[, c("cell_id", setdiff(names(metadata), "cell_id"))], file.path(exchange, "metadata.tsv"))
counts <- get_raw_counts(obj, assay)
Matrix::writeMM(Matrix::t(counts), file.path(exchange, "counts.mtx"))
writeLines(rownames(counts), file.path(exchange, "features.tsv")); writeLines(colnames(counts), file.path(exchange, "barcodes.tsv"))
embedding_manifest <- if (length(embedding_rows)) do.call(rbind, embedding_rows) else data.frame(scenario=character(), method=character(), parameters=character(), path=character(), representation_type=character())
write_tsv(embedding_manifest, file.path(exchange, "embedding_manifest.tsv"))

runs <- if (length(run_rows)) do.call(rbind, run_rows) else data.frame()
write_tsv(runs, file.path(out, "method_runs_r.tsv"))
python_methods <- any(vapply(scenarios, function(x) x$name %in% c("scvi", "scanvi", "bbknn"), logical(1)))
requested_metrics <- length(as_chr(cfg_get(config, "metrics.batch_removal", list()))) + length(as_chr(cfg_get(config, "metrics.biological_conservation", list()))) > 0L
requested_plots <- length(as_chr(cfg_get(config, "plots", list()))) > 0L
if (python_methods || requested_metrics || requested_plots) {
  prefix <- as_chr(cfg_get(config, "benchmark.python_argv_prefix", list("python3")))
  py_candidates <- c(
    file.path(dirname(script_file), "integration_python.py"),
    file.path(dirname(script_file), "..", "..", "skills", "05-scrna-benchmark-integration", "scripts", "integration_python.py")
  )
  py_script <- normalizePath(py_candidates[file.exists(py_candidates)][[1]], mustWork = TRUE)
  command <- prefix[[1]]; py_args <- c(prefix[-1], py_script, "--config", config_path, "--exchange", normalizePath(exchange, mustWork = TRUE))
  status <- system2(command, py_args)
  if (!identical(status, 0L)) stop("Python integration/metric stage failed with exit code ", status)
}

object_path <- file.path(out, "integration_benchmark_object.qs"); save_scrna_object(obj, object_path)
artifacts <- c(file.path(out, "method_runs_r.tsv"), file.path(exchange, "embedding_manifest.tsv"), object_path)
for (candidate in c("method_runs.tsv", "metric_results_long.tsv", "method_summary.tsv", "method_ranking.tsv", "design_confounding.tsv", "recommendation.md", "benchmark_embeddings.h5ad")) {
  path <- file.path(out, candidate); if (file.exists(path)) artifacts <- c(artifacts, path)
}
write_run_manifest(config, "05-scrna-benchmark-integration", out, artifacts)
