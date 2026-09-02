args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: preprocess_cluster.R CONFIG.json")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "runtime.R"))

for (pkg in c("Seurat", "SeuratObject", "jsonlite", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required")
}

config <- read_skill_config(args[[1]])
out <- prepare_output(config)
obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE))
if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object")

qc_status <- tolower(cfg_get(config, "input.qc_status", required = TRUE))
if (!qc_status %in% c("filtered", "unfiltered")) stop("input.qc_status must be filtered or unfiltered")
if (qc_status == "unfiltered" && !isTRUE(cfg_get(config, "input.allow_unfiltered", FALSE))) {
  stop("Unfiltered input requires input.allow_unfiltered=true")
}
analysis_label <- if (qc_status == "unfiltered") "exploratory_unfiltered" else "filtered"

assay <- cfg_get(config, "input.assay", Seurat::DefaultAssay(obj))
if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
Seurat::DefaultAssay(obj) <- assay
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
condition_col <- cfg_get(config, "metadata.condition")
batch_col <- cfg_get(config, "metadata.batch")
required_meta <- unique(Filter(Negate(is.null), c(sample_col, condition_col, batch_col)))
assert_metadata(obj, required_meta)
if (anyNA(obj[[sample_col, drop = TRUE]]) || any(obj[[sample_col, drop = TRUE]] == "")) stop("Sample identifiers contain missing values")

sanitize_name <- function(x) {
  if (length(x) != 1L || !grepl("^[A-Za-z0-9_-]+$", x)) stop("Invalid scenario name: ", x)
  gsub("-", "_", x, fixed = TRUE)
}
as_num <- function(x, default) as.numeric(x %||% default)
as_int <- function(x, default) as.integer(x %||% default)
resolution_text <- function(x) format(round(as.numeric(x), 3), trim = TRUE, scientific = FALSE)

save_object_atomic <- function(object, path) {
  extension <- sub("^.*\\.", "", path)
  stem <- sub(paste0("\\.", extension, "$"), "", basename(path))
  temporary <- file.path(dirname(path), paste0(".", stem, ".tmp-", Sys.getpid(), ".", extension))
  backup <- paste0(path, ".previous-", Sys.getpid())
  on.exit(unlink(c(temporary, backup), force = TRUE), add = TRUE)
  save_scrna_object(object, temporary)
  check <- load_scrna_object(temporary)
  if (!inherits(check, "Seurat") || !identical(dim(check), dim(object)) || !identical(colnames(check), colnames(object))) {
    stop("Atomic object readback validation failed")
  }
  if (file.exists(path) && !file.rename(path, backup)) stop("Could not stage the previous output object for replacement")
  if (!file.rename(temporary, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("Could not atomically publish the validated output object")
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}

write_cluster_artifacts <- function(object, scenario, cluster_col, output_dir) {
  cluster <- object[[cluster_col, drop = TRUE]]
  sizes <- as.data.frame(table(cluster = cluster), stringsAsFactors = FALSE)
  sizes$scenario <- scenario
  sizes <- sizes[c("scenario", "cluster", "Freq")]
  names(sizes)[3] <- "cells"
  size_file <- file.path(output_dir, paste0(scenario, "_cluster_sizes.tsv"))
  utils::write.table(sizes, size_file, sep = "\t", quote = FALSE, row.names = FALSE)

  composition <- as.data.frame(table(sample = object[[sample_col, drop = TRUE]], cluster = cluster), stringsAsFactors = FALSE)
  composition$scenario <- scenario
  composition <- composition[c("scenario", "sample", "cluster", "Freq")]
  names(composition)[4] <- "cells"
  comp_file <- file.path(output_dir, paste0(scenario, "_sample_cluster_counts.tsv"))
  utils::write.table(composition, comp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  reduction <- paste0(scenario, "_umap")
  if (!reduction %in% names(object@reductions)) stop("Final UMAP reduction not found: ", reduction)
  groups <- unique(c(cluster_col, intersect(unlist(cfg_get(config, "plots.group_by", list(sample_col))), colnames(object[[]]))))
  point_size <- as_num(cfg_get(config, "plots.pt_size"), 0.03)
  plot <- Seurat::DimPlot(object, reduction = reduction, group.by = groups, ncol = 1, pt.size = point_size)
  pdf_file <- file.path(output_dir, paste0(scenario, "_umap_diagnostics.pdf"))
  ggplot2::ggsave(pdf_file, plot, width = 7, height = max(5, 4 * length(groups)), units = "in", bg = "white")
  files <- c(size_file, comp_file, pdf_file)
  if (isTRUE(cfg_get(config, "plots.preview_png", FALSE))) {
    png_file <- file.path(output_dir, paste0(scenario, "_umap_diagnostics.png"))
    ggplot2::ggsave(png_file, plot, width = 7, height = max(5, 4 * length(groups)), units = "in", dpi = 200, bg = "white")
    files <- c(files, png_file)
  }
  files
}

action <- cfg_get(config, "workflow.action", "run")
if (action == "finalize_resolution") {
  scenario_cfg <- cfg_get(config, "finalize", required = TRUE)
  nm <- sanitize_name(scenario_cfg$scenario %||% "")
  resolution <- as_num(scenario_cfg$resolution, NA_real_)
  if (!is.finite(resolution) || resolution <= 0) stop("finalize.resolution must be positive")
  candidate_col <- paste0(nm, "_res.", resolution_text(resolution))
  if (!candidate_col %in% colnames(obj[[]])) stop("Candidate resolution column not found: ", candidate_col)
  final_col <- paste0(nm, "_clusters")
  obj[[final_col]] <- factor(obj[[candidate_col, drop = TRUE]])
  prior_seed <- obj@misc$scrna_preprocess$seed %||% cfg_get(config, "finalize.seed")
  attr(config, "resolved_random_seed") <- prior_seed
  obj$preprocess_qc_status <- analysis_label
  obj@misc$scrna_preprocess <- list(qc_status = qc_status, analysis_label = analysis_label,
                                    confirmed_resolution = resolution, scenario = nm, seed = prior_seed)
  cell_file <- file.path(out, "cell_assignments.tsv")
  utils::write.table(data.frame(cell_id = colnames(obj), cluster = obj[[final_col, drop = TRUE]]),
                     cell_file, sep = "\t", quote = FALSE, row.names = FALSE)
  format <- tolower(cfg_get(config, "output.object_format", "qs"))
  if (!format %in% c("qs", "rds")) stop("output.object_format must be qs or rds")
  object_file <- file.path(out, paste0("preprocessed_clustered_object.", format))
  cluster_artifacts <- write_cluster_artifacts(obj, nm, final_col, out)
  summary_file <- file.path(out, "scenario_summary.tsv")
  if (file.exists(summary_file)) {
    summary <- utils::read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
    match_row <- which(summary$scenario == nm)
    if (length(match_row) != 1L) stop("Existing scenario_summary.tsv does not contain exactly one matching scenario")
    if (!"confirmed_resolution" %in% names(summary)) summary$confirmed_resolution <- NA_real_
    if (!"qc_status" %in% names(summary)) summary$qc_status <- qc_status
    summary$confirmed_resolution[[match_row]] <- resolution
    summary$selected_resolution[[match_row]] <- resolution
    summary$status[[match_row]] <- "complete"
    summary$clusters[[match_row]] <- length(unique(obj[[final_col, drop = TRUE]]))
    summary$qc_status[[match_row]] <- qc_status
  } else {
    summary <- data.frame(
      scenario = nm, action = action, confirmed_resolution = resolution,
      selected_resolution = resolution, status = "complete",
      clusters = length(unique(obj[[final_col, drop = TRUE]])), qc_status = qc_status,
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  session_file <- file.path(out, "session_info.txt")
  writeLines(capture.output(utils::sessionInfo()), session_file)
  if (isTRUE(cfg_get(config, "output.drop_scale_data", TRUE))) {
    if ("scale.data" %in% SeuratObject::Layers(obj[[assay]])) obj[[assay]]@layers[["scale.data"]] <- NULL
  }
  save_object_atomic(obj, object_file)
  state <- list(status = "complete", action = action, scenario = nm,
                confirmed_resolution = resolution, source_column = candidate_col,
                qc_status = qc_status, analysis_label = analysis_label)
  jsonlite::write_json(state, file.path(out, "workflow_state.json"), auto_unbox = TRUE, pretty = TRUE)
  write_run_manifest(config, "06-scrna-preprocess-and-cluster", out,
                     c(object_file, cell_file, cluster_artifacts, summary_file, session_file,
                       file.path(out, "workflow_state.json")),
                     notes = c(paste0("Confirmed resolution ", resolution, " for scenario ", nm),
                               paste0("qc_status=", qc_status), paste0("analysis_label=", analysis_label)))
  quit(save = "no", status = 0)
}
if (action != "run") stop("workflow.action must be run or finalize_resolution")

seed <- as_int(cfg_get(config, "preprocessing.seed"), 1234L)
attr(config, "resolved_random_seed") <- seed
set.seed(seed)
assay_layers <- SeuratObject::Layers(obj[[assay]])
if (sum(grepl("^counts([.$]|$)", assay_layers)) > 1L) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay)
}
obj <- Seurat::NormalizeData(
  obj, assay = assay,
  normalization.method = cfg_get(config, "preprocessing.normalization_method", "LogNormalize"),
  scale.factor = as_num(cfg_get(config, "preprocessing.scale_factor"), 10000), verbose = FALSE
)

counts <- get_raw_counts(obj, assay)
if (inherits(counts, "sparseMatrix")) {
  count_values <- counts@x
  if (length(count_values) && (any(count_values < 0) || any(abs(count_values - round(count_values)) > 1e-8))) {
    stop("Raw counts must be non-negative integers")
  }
} else {
  if (min(counts) < 0) stop("Raw counts must be non-negative integers")
  check_index <- if (length(counts) > 1000000L) sample.int(length(counts), 1000000L) else seq_len(length(counts))
  if (any(abs(counts[check_index] - round(counts[check_index])) > 1e-8)) stop("Raw counts must be non-negative integers")
}

if (isTRUE(cfg_get(config, "cell_cycle.score", FALSE))) {
  species <- tolower(cfg_get(config, "cell_cycle.species", "human"))
  cc <- Seurat::cc.genes.updated.2019
  if (species == "mouse") {
    convert <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    cc$s.genes <- convert(cc$s.genes)
    cc$g2m.genes <- convert(cc$g2m.genes)
  } else if (species != "human") stop("cell_cycle.species must be human or mouse")
  present_s <- intersect(cc$s.genes, rownames(obj[[assay]]))
  present_g2m <- intersect(cc$g2m.genes, rownames(obj[[assay]]))
  if (length(present_s) < 5L || length(present_g2m) < 5L) stop("Too few cell-cycle genes found; check species and feature symbols")
  obj <- Seurat::CellCycleScoring(obj, s.features = present_s, g2m.features = present_g2m, assay = assay, set.ident = FALSE)
}

obj <- Seurat::FindVariableFeatures(
  obj, assay = assay,
  selection.method = cfg_get(config, "preprocessing.variable_feature_method", "vst"),
  nfeatures = as_int(cfg_get(config, "preprocessing.n_variable_features"), 2000L), verbose = FALSE
)

scenarios <- cfg_get(config, "scenarios", list())
if (!is.list(scenarios)) stop("scenarios must be a JSON array")
if (!length(scenarios)) stop("At least one selected scenario is required")
names_seen <- vapply(scenarios, function(x) sanitize_name(x$name %||% ""), character(1))
if (anyDuplicated(names_seen)) stop("Scenario names must be unique")
if (length(scenarios) > 4L) stop("At most four explicit scenarios are supported")

dims <- as.integer(unlist(cfg_get(config, "preprocessing.dims", as.list(1:30))))
npcs <- as_int(cfg_get(config, "preprocessing.pca_npcs"), 50L)
if (!length(dims) || min(dims) < 1L || max(dims) > npcs) stop("preprocessing.dims must fall within computed PCs")
resolution_default <- as_num(cfg_get(config, "preprocessing.resolution"), 0.4)
plot_groups <- unique(unlist(cfg_get(config, "plots.group_by", list(sample_col))))
plot_groups <- intersect(plot_groups, colnames(obj[[]]))
pt_size <- as_num(cfg_get(config, "plots.pt_size"), 0.03)

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n_pairs <- choose2(sum(tab))
  index <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  expected <- row_pairs * col_pairs / n_pairs
  maximum <- (row_pairs + col_pairs) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

cluster_graph <- function(graph, resolution, random_seed) {
  labels <- Seurat::FindClusters(
    object = SeuratObject::as.Graph(as(graph, "dgCMatrix")), resolution = resolution,
    algorithm = 1, random.seed = random_seed, verbose = FALSE
  )
  stats::setNames(as.character(labels[[1]]), rownames(labels))
}

resolution_stability <- function(graph, resolutions, n_subsamples, subsample_frac, random_seed) {
  graph_cells <- colnames(graph)
  set.seed(random_seed)
  sub_n <- max(30L, as.integer(round(length(graph_cells) * subsample_frac)))
  subsamples <- lapply(seq_len(n_subsamples), function(i) sort(sample(graph_cells, sub_n)))
  references <- lapply(resolutions, function(resolution) cluster_graph(graph, resolution, random_seed))
  rows <- lapply(seq_along(resolutions), function(index) {
    resolution <- resolutions[[index]]
    reference <- references[[index]]
    aris <- vapply(subsamples, function(cells) {
      subgraph <- graph[cells, cells, drop = FALSE]
      labels <- cluster_graph(subgraph, resolution, random_seed)
      adjusted_rand(reference[cells], labels[cells])
    }, numeric(1))
    cluster_sizes <- table(reference)
    data.frame(resolution = resolution, n_clusters = length(cluster_sizes),
               stability_mean = mean(aris), stability_sd = stats::sd(aris),
               stability_conservative = mean(aris) - stats::sd(aris),
               min_cluster_cells = min(cluster_sizes), max_cluster_fraction = max(cluster_sizes) / length(reference),
               clusters_lt50 = sum(cluster_sizes < 50), clusters_lt100 = sum(cluster_sizes < 100))
  })
  result <- do.call(rbind, rows)
  result$ari_previous <- c(NA_real_, vapply(seq.int(2L, length(references)), function(index) {
    adjusted_rand(references[[index - 1L]], references[[index]])
  }, numeric(1)))
  result
}

scenario_rows <- list()
artifacts <- character()
for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  nm <- sanitize_name(scenario$name)
  key_stem <- gsub("[^A-Za-z0-9]", "", nm)
  regress <- unlist(scenario$regress_variables %||% list(), use.names = FALSE)
  if (length(regress)) {
    assert_metadata(obj, regress)
    bad <- regress[!vapply(obj[[]][regress], is.numeric, logical(1))]
    if (length(bad)) stop("Regression variables must be numeric: ", paste(bad, collapse = ", "))
    flat <- regress[vapply(obj[[]][regress], function(x) length(unique(x[is.finite(x)])) < 2L, logical(1))]
    if (length(flat)) stop("Regression variables lack variation: ", paste(flat, collapse = ", "))
  }
  pca_name <- paste0(nm, "_pca")
  pca_key <- paste0(key_stem, "PC_")
  umap_name <- paste0(nm, "_umap")
  reduction_use <- pca_name
  obj <- Seurat::ScaleData(obj, assay = assay, vars.to.regress = if (length(regress)) regress else NULL, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, assay = assay, npcs = npcs, reduction.name = pca_name,
                        reduction.key = pca_key, seed.use = seed, verbose = FALSE)

  harmony_enabled <- isTRUE(cfg_get(scenario, "harmony.enabled", FALSE))
  harmony_theta <- NA_real_
  harmony_groups <- character()
  if (harmony_enabled) {
    if (!requireNamespace("harmony", quietly = TRUE)) stop("Package 'harmony' is required for scenario ", nm)
    harmony_groups <- unlist(cfg_get(scenario, "harmony.group_by", list()), use.names = FALSE)
    if (!length(harmony_groups)) stop("Harmony group_by is required for scenario ", nm)
    assert_metadata(obj, harmony_groups)
    harmony_theta <- as_num(cfg_get(scenario, "harmony.theta"), 2)
    reduction_use <- paste0(nm, "_harmony")
    obj <- harmony::RunHarmony(
      object = obj, group.by.vars = harmony_groups, reduction.use = pca_name,
      dims.use = dims, theta = harmony_theta, reduction.save = reduction_use, verbose = FALSE
    )
  }

  nn_name <- paste0(nm, "_nn")
  snn_name <- paste0(nm, "_snn")
  cluster_name <- paste0(nm, "_clusters")
  clustering_mode <- cfg_get(scenario, "clustering.mode", "fixed")
  resolution <- as_num(cfg_get(scenario, "clustering.resolution", scenario$resolution), resolution_default)
  k_param <- as_int(cfg_get(scenario, "neighbors.k_param"), 20L)
  if (!is.finite(k_param) || k_param < 1L || k_param >= ncol(obj)) stop("neighbors.k_param must be between 1 and the number of cells minus one")
  umap_n_neighbors <- as_int(cfg_get(scenario, "umap.n_neighbors"), 30L)
  umap_min_dist <- as_num(cfg_get(scenario, "umap.min_dist"), 0.3)
  umap_metric <- cfg_get(scenario, "umap.metric", "cosine")
  umap_method <- cfg_get(scenario, "umap.method", "uwot")
  if (!is.finite(umap_n_neighbors) || umap_n_neighbors < 2L || umap_n_neighbors >= ncol(obj)) stop("umap.n_neighbors must be between 2 and the number of cells minus one")
  if (!is.finite(umap_min_dist) || umap_min_dist < 0 || umap_min_dist > 1) stop("umap.min_dist must be between 0 and 1")
  if (!is.character(umap_metric) || length(umap_metric) != 1L || !nzchar(umap_metric)) stop("umap.metric must be a non-empty string")
  if (!umap_method %in% c("uwot", "umap-learn")) stop("umap.method must be uwot or umap-learn")
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_use, dims = dims, k.param = k_param,
                               graph.name = c(nn_name, snn_name), verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = reduction_use, dims = dims, reduction.name = umap_name,
                         reduction.key = paste0(key_stem, "UMAP_"), seed.use = seed,
                         n.neighbors = umap_n_neighbors, min.dist = umap_min_dist,
                         metric = umap_metric, umap.method = umap_method, verbose = FALSE)

  recommended_resolution <- NA_real_
  selected_resolution <- resolution
  final_available <- TRUE
  if (clustering_mode == "fixed") {
    obj <- Seurat::FindClusters(obj, graph.name = snn_name, resolution = resolution,
                                cluster.name = cluster_name, random.seed = seed, verbose = FALSE)
  } else if (clustering_mode == "scan") {
    resolutions <- sort(unique(as.numeric(unlist(cfg_get(scenario, "clustering.resolutions", as.list(seq(0.1, 1, 0.1)))))))
    resolutions <- resolutions[is.finite(resolutions) & resolutions > 0]
    if (length(resolutions) < 2L) stop("Resolution scan requires at least two positive candidates")
    graph <- obj@graphs[[snn_name]]
    candidate_cols <- vapply(resolutions, function(candidate) {
      col <- paste0(nm, "_res.", resolution_text(candidate))
      labels <- cluster_graph(graph, candidate, seed)
      obj[[col]] <<- factor(labels[colnames(obj)])
      col
    }, character(1))
    stability <- resolution_stability(
      graph, resolutions,
      n_subsamples = as_int(cfg_get(scenario, "clustering.stability.n_subsamples"), 5L),
      subsample_frac = as_num(cfg_get(scenario, "clustering.stability.subsample_frac"), 0.8),
      random_seed = seed
    )
    min_clusters <- as_int(cfg_get(scenario, "clustering.stability.min_clusters"), 3L)
    eligible <- stability[stability$n_clusters >= min_clusters, , drop = FALSE]
    if (!nrow(eligible)) stop("No scanned resolution produced the minimum number of clusters")
    recommended_resolution <- eligible$resolution[[which.max(eligible$stability_conservative)]]
    stability_file <- file.path(out, paste0(nm, "_resolution_stability.tsv"))
    utils::write.table(stability, stability_file, sep = "\t", quote = FALSE, row.names = FALSE)
    p_stability <- ggplot2::ggplot(stability, ggplot2::aes(x = resolution, y = stability_mean)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = stability_mean - stability_sd, ymax = stability_mean + stability_sd), width = 0.02) +
      ggplot2::geom_vline(xintercept = recommended_resolution, linetype = 2, colour = "red") + ggplot2::theme_linedraw()
    stability_plot <- file.path(out, paste0(nm, "_resolution_stability.png"))
    ggplot2::ggsave(stability_plot, p_stability, width = 7, height = 4, dpi = 300, bg = "white")
    p_grid <- Seurat::DimPlot(obj, reduction = umap_name, group.by = candidate_cols,
                              ncol = min(5, length(candidate_cols)), label = TRUE, pt.size = pt_size) & Seurat::NoLegend()
    grid_file <- file.path(out, paste0(nm, "_umap_clusters_by_resolution.png"))
    ggplot2::ggsave(grid_file, p_grid, width = 20, height = ceiling(length(candidate_cols) / 5) * 4, dpi = 300, bg = "white")
    artifacts <- c(artifacts, stability_file, stability_plot, grid_file)
    if (requireNamespace("clustree", quietly = TRUE)) {
      p_tree <- clustree::clustree(obj[[]], prefix = paste0(nm, "_res.")) +
        ggplot2::guides(edge_colour = "none", edge_alpha = "none")
      tree_file <- file.path(out, paste0(nm, "_clustree_resolution.png"))
      tree_error <- tryCatch({
        ggplot2::ggsave(tree_file, p_tree, width = 13, height = 11, dpi = 300, bg = "white")
        NULL
      }, error = function(e) conditionMessage(e))
      if (is.null(tree_error)) artifacts <- c(artifacts, tree_file)
      else writeLines(tree_error, file.path(out, paste0(nm, "_clustree_error.txt")))
    } else {
      writeLines("Package 'clustree' is unavailable", file.path(out, paste0(nm, "_clustree_error.txt")))
    }
    selection <- cfg_get(scenario, "clustering.selection", "review")
    if (selection == "recommended") selected_resolution <- recommended_resolution
    else if (selection == "fixed") selected_resolution <- resolution
    else if (selection == "review") final_available <- FALSE
    else stop("clustering.selection must be review, recommended, or fixed")
    if (final_available) {
      selected_col <- paste0(nm, "_res.", resolution_text(selected_resolution))
      if (!selected_col %in% candidate_cols) stop("Selected resolution was not included in the scan")
      obj[[cluster_name]] <- factor(obj[[selected_col, drop = TRUE]])
    }
  } else stop("clustering.mode must be fixed or scan")

  cluster <- if (final_available) obj[[cluster_name, drop = TRUE]] else NULL
  if (final_available) {
    artifacts <- c(artifacts, write_cluster_artifacts(obj, nm, cluster_name, out))
  }
  scenario_rows[[i]] <- data.frame(
    scenario = nm, regress_variables = paste(regress, collapse = ","), harmony = harmony_enabled,
    harmony_group_by = paste(harmony_groups, collapse = ","), theta = harmony_theta,
    reduction = reduction_use, dims = paste(range(dims), collapse = "-"), clustering_mode = clustering_mode,
    recommended_resolution = recommended_resolution,
    selected_resolution = if (final_available) selected_resolution else NA_real_,
    k_param = k_param, umap_n_neighbors = umap_n_neighbors, umap_min_dist = umap_min_dist,
    umap_metric = umap_metric, umap_method = umap_method,
    status = if (final_available) "complete" else "awaiting_resolution_confirmation",
    clusters = if (final_available) length(unique(cluster)) else NA_integer_, stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, scenario_rows)
summary_file <- file.path(out, "scenario_summary.tsv")
utils::write.table(summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

scenario_names <- vapply(scenarios, function(x) sanitize_name(x$name), character(1))
complete_names <- scenario_names[paste0(scenario_names, "_clusters") %in% colnames(obj[[]])]
pairs <- if (length(complete_names) > 1L) utils::combn(complete_names, 2, simplify = FALSE) else list()
similarity <- if (length(pairs)) do.call(rbind, lapply(pairs, function(pair) {
  data.frame(
    scenario_a = pair[[1]], scenario_b = pair[[2]],
    adjusted_rand_index = adjusted_rand(
      obj[[paste0(pair[[1]], "_clusters"), drop = TRUE]],
      obj[[paste0(pair[[2]], "_clusters"), drop = TRUE]]
    ), stringsAsFactors = FALSE
  )
})) else data.frame(scenario_a = character(), scenario_b = character(), adjusted_rand_index = numeric())
similarity_file <- file.path(out, "scenario_cluster_similarity.tsv")
utils::write.table(similarity, similarity_file, sep = "\t", quote = FALSE, row.names = FALSE)
first_scenario <- sanitize_name(scenarios[[1]]$name)
elbow <- Seurat::ElbowPlot(obj, reduction = paste0(first_scenario, "_pca"), ndims = npcs)
elbow_file <- file.path(out, paste0(first_scenario, "_elbow.png"))
ggplot2::ggsave(elbow_file, elbow, width = 4, height = 4, units = "in", dpi = 300, bg = "white")

format <- tolower(cfg_get(config, "output.object_format", "qs"))
if (!format %in% c("qs", "rds")) stop("output.object_format must be qs or rds")
if (isTRUE(cfg_get(config, "output.drop_scale_data", TRUE))) {
  if ("scale.data" %in% SeuratObject::Layers(obj[[assay]])) obj[[assay]]@layers[["scale.data"]] <- NULL
}
object_file <- file.path(out, paste0("preprocessed_clustered_object.", format))
obj$preprocess_qc_status <- analysis_label
obj@misc$scrna_preprocess <- list(qc_status = qc_status, analysis_label = analysis_label,
                                  seed = seed, scenarios = names_seen)
save_object_atomic(obj, object_file)
cell_file <- file.path(out, "cell_assignments.tsv")
cluster_cols <- paste0(vapply(scenarios, function(x) sanitize_name(x$name), character(1)), "_clusters")
cluster_cols <- intersect(cluster_cols, colnames(obj[[]]))
if (length(cluster_cols)) {
  cell_table <- data.frame(cell_id = colnames(obj), obj[[]][cluster_cols], check.names = FALSE)
  utils::write.table(cell_table, cell_file, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cell_file <- character()
}
writeLines(capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"))

artifacts <- unique(c(artifacts, summary_file, similarity_file, elbow_file, object_file, cell_file, file.path(out, "session_info.txt")))
state_status <- if (any(summary$status == "awaiting_resolution_confirmation")) "awaiting_resolution_confirmation" else "complete"
state_scenarios <- lapply(seq_len(nrow(summary)), function(index) as.list(summary[index, , drop = FALSE]))
state <- list(status = state_status, action = action, scenarios = state_scenarios,
              qc_status = qc_status, analysis_label = analysis_label,
              next_action = if (state_status == "awaiting_resolution_confirmation") "Review resolution plots and run finalize_resolution" else "none")
jsonlite::write_json(state, file.path(out, "workflow_state.json"), auto_unbox = TRUE, pretty = TRUE)
artifacts <- unique(c(artifacts, file.path(out, "workflow_state.json")))
write_run_manifest(config, "06-scrna-preprocess-and-cluster", out, artifacts,
                   notes = c(
                     paste0(length(scenarios), " user-selected scenario(s) executed"),
                     "No additional scenario was inserted or selected automatically",
                     paste0("qc_status=", qc_status), paste0("analysis_label=", analysis_label)
                   ))
