args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript score_programs.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
config <- read_skill_config(args[[1]])

require_pkg <- function(package, reason = "") {
  if (!requireNamespace(package, quietly = TRUE)) {
    suffix <- if (nzchar(reason)) paste0(" ", reason) else ""
    stop("Package '", package, "' is required.", suffix)
  }
}
require_pkg("Seurat")
require_pkg("SeuratObject")

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^[_.-]+|[_.-]+$", "", x)
  if (!nzchar(x)) stop("A task or signature name becomes empty after sanitization")
  x
}

as_character_vector <- function(x) unique(as.character(unlist(x, use.names = FALSE)))

digest_value <- function(x) {
  path <- tempfile("scrna-program-digest-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3)
  unname(tools::md5sum(path))
}

get_layer <- function(obj, assay, layer) {
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    layers <- SeuratObject::Layers(obj[[assay]])
    exact <- layers[layers == layer]
    split <- layers[grepl(paste0("^", layer, "[.]"), layers)]
    if (!length(exact) && length(split) > 1L) {
      obj <- SeuratObject::JoinLayers(obj, assay = assay, layers = split, new = layer)
      exact <- layer
    }
    if (!length(exact)) stop("Layer not found: ", assay, "/", layer)
    return(list(object = obj, matrix = SeuratObject::LayerData(obj, assay = assay, layer = layer)))
  }
  slot <- if (layer %in% c("counts", "data", "scale.data")) layer else stop("Unsupported Seurat v4 slot: ", layer)
  list(object = obj, matrix = SeuratObject::GetAssayData(obj, assay = assay, slot = slot))
}

has_layer <- function(obj, assay, layer) {
  if (!assay %in% names(obj@assays)) return(FALSE)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    layers <- SeuratObject::Layers(obj[[assay]])
    return(any(layers == layer | grepl(paste0("^", layer, "[.]"), layers)))
  }
  nrow(SeuratObject::GetAssayData(obj, assay = assay, slot = layer)) > 0L
}

read_gmt <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines[nzchar(lines)], "\t", fixed = FALSE)
  fields <- fields[lengths(fields) >= 3L]
  sets <- lapply(fields, function(x) unique(x[-c(1, 2)]))
  names(sets) <- vapply(fields, `[[`, character(1), 1L)
  sets
}

write_gmt <- function(sets, path) {
  lines <- vapply(names(sets), function(nm) paste(c(nm, "generated-by-scrna-score-programs", sets[[nm]]), collapse = "\t"), character(1))
  writeLines(lines, path)
}

resolve_gene_sets <- function(spec, species, resource_dir, task_name) {
  source <- tolower(cfg_get(spec, "source", required = TRUE))
  provenance <- list(source = source)
  if (source == "inline") {
    sets <- spec$sets
    if (!is.list(sets) || !length(sets) || is.null(names(sets))) stop("inline gene_sets.sets must be a named object")
    sets <- lapply(sets, as_character_vector)
  } else if (source == "gmt") {
    path <- cfg_get(spec, "path", required = TRUE)
    sets <- read_gmt(path)
    provenance$path <- normalizePath(path, mustWork = TRUE)
    provenance$md5 <- unname(tools::md5sum(provenance$path))
  } else if (source == "msigdb") {
    require_pkg("msigdbr", "Provide a cached GMT with source=gmt if msigdbr is unavailable.")
    collection <- cfg_get(spec, "collection", if (species == "mouse") "MH" else "H")
    subcollection <- cfg_get(spec, "subcollection")
    db_species <- if (species == "mouse") "MM" else "HS"
    organism <- if (species == "mouse") "Mus musculus" else "Homo sapiens"
    cache_name <- paste(c("msigdb", db_species, collection, subcollection), collapse = "_")
    cache_name <- paste0(sanitize_name(cache_name), ".gmt")
    cache_path <- file.path(resource_dir, cache_name)
    if (file.exists(cache_path)) {
      sets <- read_gmt(cache_path)
    } else {
      call <- list(db_species = db_species, species = organism, collection = collection)
      if (!is.null(subcollection) && nzchar(subcollection)) call$subcollection <- subcollection
      tab <- do.call(msigdbr::msigdbr, call)
      if (!nrow(tab)) stop("msigdbr returned no genes for collection ", collection)
      sets <- split(tab$gene_symbol, tab$gs_name)
      sets <- lapply(sets, unique)
      write_gmt(sets, cache_path)
    }
    provenance$collection <- collection
    provenance$subcollection <- subcollection
    provenance$path <- normalizePath(cache_path, mustWork = TRUE)
    provenance$md5 <- unname(tools::md5sum(cache_path))
    provenance$package_version <- as.character(utils::packageVersion("msigdbr"))
  } else if (source == "scmetabolism") {
    require_pkg("scMetabolism", "Provide the package GMT explicitly with source=gmt if needed.")
    if (species != "human") stop("The bundled scMetabolism GMT is human-symbol based; provide a reviewed species-matched GMT for mouse")
    database <- toupper(cfg_get(spec, "database", "KEGG"))
    filename <- if (database == "KEGG") "KEGG_metabolism_nc.gmt" else if (database == "REACTOME") "REACTOME_metabolism_nc.gmt" else stop("scmetabolism database must be KEGG or REACTOME")
    path <- system.file("data", filename, package = "scMetabolism")
    if (!nzchar(path)) stop("scMetabolism resource not found: ", filename)
    sets <- read_gmt(path)
    provenance$database <- database
    provenance$path <- normalizePath(path, mustWork = TRUE)
    provenance$md5 <- unname(tools::md5sum(path))
    provenance$package_version <- as.character(utils::packageVersion("scMetabolism"))
  } else stop("Unsupported gene set source: ", source)
  names(sets) <- make.unique(vapply(names(sets), sanitize_name, character(1)))
  sets <- lapply(sets, function(x) unique(x[!is.na(x) & nzchar(x)]))
  sets <- sets[lengths(sets) > 0L]
  if (!length(sets)) stop("No non-empty gene sets resolved for task ", task_name)
  list(sets = sets, provenance = provenance)
}

coverage_table <- function(sets, features, task_name, method, min_genes, min_fraction) {
  rows <- lapply(names(sets), function(nm) {
    input <- unique(sets[[nm]])
    matched <- intersect(input, features)
    missing <- setdiff(input, features)
    fraction <- if (length(input)) length(matched) / length(input) else 0
    status <- if (length(matched) < min_genes) "insufficient_genes" else if (fraction < min_fraction) "low_fraction" else "ok"
    data.frame(task = task_name, method = method, signature = nm, input_genes = length(input),
               matched_genes = length(matched), coverage_fraction = fraction, status = status,
               missing_genes = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

score_vision <- function(mat, sets, cores) {
  require_pkg("VISION")
  totals <- Matrix::colSums(mat)
  if (any(!is.finite(totals) | totals <= 0)) stop("VISION requires positive finite library sizes for every cell")
  scaled <- mat %*% Matrix::Diagonal(x = stats::median(totals) / totals)
  dimnames(scaled) <- dimnames(mat)
  gmt <- tempfile(fileext = ".gmt")
  on.exit(unlink(gmt), add = TRUE)
  write_gmt(sets, gmt)
  old_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_cores), add = TRUE)
  options(mc.cores = cores)
  model <- VISION::Vision(scaled, signatures = gmt)
  model <- VISION::calcSignatureScores(model)
  scores <- t(as.matrix(model@SigScores))
  if (ncol(scores) != ncol(mat)) stop("Unexpected VISION score dimensions: ", paste(dim(scores), collapse = "x"), "; expected signatures x ", ncol(mat))
  retained <- intersect(names(sets), rownames(scores))
  if (!length(retained)) stop("VISION signature names could not be reconciled with the requested sets")
  scores <- scores[retained, , drop = FALSE]
  colnames(scores) <- colnames(mat)
  scores
}

score_aucell <- function(mat, sets, cores, auc_max_rank = NULL) {
  require_pkg("AUCell")
  rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1L, plotStats = FALSE, verbose = FALSE)
  args <- list(geneSets = sets, rankings = rankings, nCores = 1L)
  if (!is.null(auc_max_rank)) args$aucMaxRank <- as.integer(auc_max_rank)
  auc <- do.call(AUCell::AUCell_calcAUC, args)
  as.matrix(SummarizedExperiment::assay(auc))
}

score_ucell <- function(obj, assay, layer, sets, name) {
  require_pkg("UCell")
  before <- colnames(obj[[]])
  obj <- UCell::AddModuleScore_UCell(obj, features = sets, assay = assay, slot = layer, name = paste0(name, "__"))
  added <- setdiff(colnames(obj[[]]), before)
  if (length(added) != length(sets)) stop("UCell returned an unexpected number of score columns")
  scores <- t(as.matrix(obj[[]][, added, drop = FALSE]))
  rownames(scores) <- names(sets)
  list(object = obj, scores = scores)
}

score_addmodule <- function(obj, assay, layer, sets, name, seed, nbin, ctrl) {
  before <- colnames(obj[[]])
  obj <- Seurat::AddModuleScore(obj, features = unname(sets), assay = assay, slot = layer,
                                name = paste0(name, "__"), seed = seed, nbin = nbin, ctrl = ctrl)
  added <- setdiff(colnames(obj[[]]), before)
  if (length(added) != length(sets)) stop("AddModuleScore returned an unexpected number of score columns")
  scores <- t(as.matrix(obj[[]][, added, drop = FALSE]))
  rownames(scores) <- names(sets)
  list(object = obj, scores = scores)
}

score_progeny <- function(mat, species, top, scale) {
  require_pkg("progeny")
  organism <- if (species == "mouse") "Mouse" else "Human"
  acts <- progeny::progeny(as.matrix(mat), organism = organism, top = top, scale = scale, perm = 1)
  acts <- as.matrix(acts)
  if (ncol(acts) == ncol(mat)) return(acts[, colnames(mat), drop = FALSE])
  if (nrow(acts) == ncol(mat)) return(t(acts[colnames(mat), , drop = FALSE]))
  stop("Unexpected PROGENy output dimensions")
}

summarize_scores <- function(scores, meta, groups, task_name) {
  if (!length(groups)) return(data.frame())
  missing <- setdiff(groups, colnames(meta))
  if (length(missing)) stop("Summary metadata columns not found: ", paste(missing, collapse = ", "))
  group_frame <- meta[colnames(scores), groups, drop = FALSE]
  if (anyNA(group_frame)) stop("Summary metadata columns contain missing values")
  key <- interaction(group_frame, drop = TRUE, lex.order = TRUE)
  means <- t(vapply(split(seq_len(ncol(scores)), key), function(idx) Matrix::rowMeans(scores[, idx, drop = FALSE]), numeric(nrow(scores))))
  keys <- do.call(rbind, lapply(split(seq_len(nrow(group_frame)), key), function(idx) group_frame[idx[1], , drop = FALSE]))
  long <- data.frame(keys[rep(seq_len(nrow(keys)), each = nrow(scores)), , drop = FALSE],
                     task = task_name, signature = rep(rownames(scores), times = nrow(keys)),
                     mean_score = as.vector(t(means)), n_cells = rep(as.integer(table(key)), each = nrow(scores)),
                     row.names = NULL, check.names = FALSE)
  long
}

obj_path <- cfg_get(config, "input.object", required = TRUE)
obj <- load_scrna_object(obj_path, "auto")
out <- prepare_output(config)
dir.create(file.path(out, "scores"), showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)
resource_dir <- cfg_get(config, "resources.cache_dir", file.path(out, "resource_cache"))
dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)

species <- tolower(cfg_get(config, "species", required = TRUE))
if (!species %in% c("human", "mouse")) stop("species must be human or mouse")
assay_default <- cfg_get(config, "input.assay", if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj))
layer_default <- cfg_get(config, "input.layer", "counts")
tasks <- config$tasks
if (!is.list(tasks) || !length(tasks)) stop("tasks must be a non-empty array")
task_names <- vapply(tasks, function(x) sanitize_name(cfg_get(x, "name", required = TRUE)), character(1))
if (anyDuplicated(task_names)) stop("Task names must be unique after sanitization")
summary_groups <- as_character_vector(cfg_get(config, "summarize_by", list()))
if (length(summary_groups)) assert_metadata(obj, summary_groups)

input_fingerprint <- unname(tools::md5sum(normalizePath(obj_path, mustWork = TRUE)))
coverage_all <- list()
summary_all <- list()
task_manifest <- list()
output_assays <- character()
feature_mappings <- list()
artifacts <- character()
seed <- as.integer(cfg_get(config, "random_seed", 1L))
cores_default <- as.integer(cfg_get(config, "cores", 1L))

for (i in seq_along(tasks)) {
  task <- tasks[[i]]
  task_name <- task_names[[i]]
  method <- tolower(cfg_get(task, "method", required = TRUE))
  assay <- cfg_get(task, "assay", assay_default)
  layer <- cfg_get(task, "layer", if (method == "progeny") "data" else layer_default)
  normalized_in_memory <- FALSE
  if (layer == "data" && !has_layer(obj, assay, "data")) {
    if (!isTRUE(cfg_get(task, "normalize_if_missing", cfg_get(config, "input.normalize_if_missing", TRUE)))) {
      stop("Normalized data layer is missing for task ", task_name)
    }
    obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
    normalized_in_memory <- TRUE
  }
  loaded <- get_layer(obj, assay, layer)
  obj <- loaded$object
  mat <- loaded$matrix
  if (!identical(colnames(mat), colnames(obj))) mat <- mat[, colnames(obj), drop = FALSE]
  matrix_values <- if (methods::is(mat, "sparseMatrix")) mat@x else as.vector(mat)
  if (any(!is.finite(matrix_values))) stop("Expression matrix contains non-finite values for task ", task_name)
  cores <- as.integer(cfg_get(task, "cores", cores_default))
  if (is.na(cores) || cores < 1L) stop("cores must be positive")
  min_genes <- as.integer(cfg_get(task, "coverage.min_genes", 3L))
  min_fraction <- as.numeric(cfg_get(task, "coverage.min_fraction", 0.1))
  on_insufficient <- tolower(cfg_get(task, "coverage.on_insufficient", "skip"))
  gene_sets <- NULL
  provenance <- list(source = "progeny_model")
  coverage <- NULL
  if (method != "progeny") {
    resolved <- resolve_gene_sets(task$gene_sets, species, resource_dir, task_name)
    gene_sets <- resolved$sets
    provenance <- resolved$provenance
    coverage <- coverage_table(gene_sets, rownames(mat), task_name, method, min_genes, min_fraction)
    bad <- coverage$status == "insufficient_genes"
    if (any(bad) && on_insufficient == "error") stop("Insufficient matched genes in task ", task_name, ": ", paste(coverage$signature[bad], collapse = ", "))
    if (!on_insufficient %in% c("skip", "error")) stop("coverage.on_insufficient must be skip or error")
    keep <- coverage$status != "insufficient_genes"
    gene_sets <- Map(intersect, gene_sets[coverage$signature[keep]], MoreArgs = list(y = rownames(mat)))
    if (!length(gene_sets)) stop("No scoreable signatures remain for task ", task_name)
  }
  method_package <- c(vision = "VISION", aucell = "AUCell", ucell = "UCell", addmodulescore = "Seurat", progeny = "progeny")[[method]]
  method_version <- if (!is.null(method_package) && requireNamespace(method_package, quietly = TRUE)) as.character(utils::packageVersion(method_package)) else NA_character_
  task_key <- digest_value(list(input = input_fingerprint, cells = colnames(mat), features = rownames(mat),
                                task = task, provenance = provenance, package = method_package, package_version = method_version))
  cache_file <- file.path(out, "scores", paste0(task_name, "_", task_key, ".rds"))
  cache_hit <- isTRUE(cfg_get(config, "cache.enabled", TRUE)) && file.exists(cache_file)
  if (cache_hit) {
    scores <- readRDS(cache_file)
    if (!identical(colnames(scores), colnames(obj))) stop("Cached score cells do not match the input object")
  } else {
    if (method == "vision") scores <- score_vision(mat, gene_sets, cores)
    else if (method == "aucell") scores <- score_aucell(mat, gene_sets, cores, cfg_get(task, "parameters.auc_max_rank"))
    else if (method == "ucell") {
      scored <- score_ucell(obj, assay, layer, gene_sets, task_name)
      obj <- scored$object; scores <- scored$scores
    } else if (method == "addmodulescore") {
      scored <- score_addmodule(obj, assay, layer, gene_sets, task_name, seed,
                                as.integer(cfg_get(task, "parameters.nbin", 24L)), as.integer(cfg_get(task, "parameters.ctrl", 100L)))
      obj <- scored$object; scores <- scored$scores
    } else if (method == "progeny") {
      scores <- score_progeny(mat, species, as.integer(cfg_get(task, "parameters.top", 500L)), isTRUE(cfg_get(task, "parameters.scale", FALSE)))
    } else stop("Unsupported method: ", method)
    scores <- as.matrix(scores)
    if (!identical(colnames(scores), colnames(obj))) scores <- scores[, colnames(obj), drop = FALSE]
    if (any(!is.finite(scores))) stop("Non-finite scores produced for task ", task_name)
    saveRDS(scores, cache_file, compress = FALSE)
  }
  if (!is.null(coverage)) {
    dropped <- !coverage$signature %in% rownames(scores) & coverage$status != "insufficient_genes"
    coverage$status[dropped] <- "method_dropped"
    coverage_all[[length(coverage_all) + 1L]] <- coverage
  }
  assay_name <- tail(make.unique(c(names(obj@assays), paste0("program_", task_name))), 1L)
  assay_scores <- scores
  rownames(assay_scores) <- make.unique(gsub("_", "-", rownames(scores), fixed = TRUE))
  feature_mappings[[length(feature_mappings) + 1L]] <- data.frame(
    task = task_name, matrix_signature = rownames(scores), assay_feature = rownames(assay_scores), stringsAsFactors = FALSE
  )
  score_assay <- SeuratObject::CreateAssayObject(data = assay_scores)
  extra_cells <- setdiff(colnames(score_assay), colnames(obj))
  if (length(extra_cells)) stop("Score assay contains cells absent from Seurat object: ", paste(head(extra_cells), collapse = ","))
  if (!identical(colnames(score_assay), colnames(obj))) stop("Score assay cell order does not match the Seurat object")
  obj[[assay_name]] <- score_assay
  output_assays <- c(output_assays, assay_name)
  score_file <- file.path(out, "scores", paste0(task_name, "_scores.tsv.gz"))
  con <- gzfile(score_file, "wt")
  utils::write.table(data.frame(signature = rownames(scores), scores, check.names = FALSE), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  artifacts <- c(artifacts, score_file, cache_file)
  if (length(summary_groups)) summary_all[[length(summary_all) + 1L]] <- summarize_scores(scores, obj[[]], summary_groups, task_name)
  task_manifest[[task_name]] <- list(method = method, assay = assay, layer = layer, assay_output = assay_name,
                                     n_signatures = nrow(scores), n_cells = ncol(scores), cache_key = task_key,
                                     cache_hit = cache_hit, normalized_in_memory = normalized_in_memory,
                                     package = method_package, package_version = method_version, resource = provenance)
}

coverage_path <- file.path(out, "signature_coverage.tsv")
missing_output_assays <- setdiff(output_assays, names(obj@assays))
if (length(missing_output_assays)) stop("Score assays were lost before output: ", paste(missing_output_assays, collapse = ","))
coverage_out <- if (length(coverage_all)) do.call(rbind, coverage_all) else data.frame(task = character(), method = character(), signature = character(), input_genes = integer(), matched_genes = integer(), coverage_fraction = numeric(), status = character(), missing_genes = character())
utils::write.table(coverage_out, coverage_path, sep = "\t", quote = FALSE, row.names = FALSE)
summary_path <- file.path(out, "score_summary.tsv")
summary_out <- if (length(summary_all)) do.call(rbind, summary_all) else data.frame()
utils::write.table(summary_out, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
mapping_path <- file.path(out, "assay_feature_mapping.tsv")
mapping_out <- do.call(rbind, feature_mappings)
utils::write.table(mapping_out, mapping_path, sep = "\t", quote = FALSE, row.names = FALSE)

plot_paths <- character()
if (nrow(summary_out) && length(summary_groups)) {
  for (task_name in unique(summary_out$task)) {
    tab <- summary_out[summary_out$task == task_name, , drop = FALSE]
    labels <- apply(tab[, summary_groups, drop = FALSE], 1L, paste, collapse = " | ")
    mat_plot <- stats::xtabs(tab$mean_score ~ tab$signature + labels)
    if (nrow(mat_plot) > 1L && ncol(mat_plot) > 1L) {
      plot_path <- file.path(out, "figures", paste0(task_name, "_group_mean_heatmap.png"))
      grDevices::png(plot_path, width = min(6000, max(1200, 45 * ncol(mat_plot))), height = min(5000, max(900, 24 * nrow(mat_plot))), res = 150)
      stats::heatmap(mat_plot, scale = "row", margins = c(10, 12), main = paste(task_name, "group mean scores"))
      grDevices::dev.off()
      plot_paths <- c(plot_paths, plot_path)
    }
  }
}

format <- tolower(cfg_get(config, "output.object_format", "qs"))
if (!format %in% c("qs", "rds")) stop("output.object_format must be qs or rds")
object_path <- file.path(out, paste0("program_scored_object.", format))
save_scrna_object(obj, object_path)
session_path <- file.path(out, "session_info.txt")
writeLines(capture.output(utils::sessionInfo()), session_path)
task_manifest_path <- file.path(out, "task_manifest.json")
jsonlite::write_json(task_manifest, task_manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
artifacts <- unique(c(artifacts, coverage_path, summary_path, mapping_path, plot_paths, object_path, session_path, task_manifest_path))
write_run_manifest(config, "scrna-score-programs", out, artifacts,
                   notes = c(paste0("species=", species), paste0("tasks=", paste(task_names, collapse = ",")),
                             "Scores are descriptive at cell level; use biological samples for condition-level inference"))
