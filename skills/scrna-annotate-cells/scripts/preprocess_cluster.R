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
  cell_file <- file.path(out, "cell_assignments.tsv")
  utils::write.table(data.frame(cell_id = colnames(obj), cluster = obj[[final_col, drop = TRUE]]),
                     cell_file, sep = "\t", quote = FALSE, row.names = FALSE)
  format <- tolower(cfg_get(config, "output.object_format", "qs"))
  object_file <- file.path(out, paste0("preprocessed_clustered_object.", format))
  save_scrna_object(obj, object_file)
  state <- list(status = "complete", action = action, scenario = nm,
                confirmed_resolution = resolution, source_column = candidate_col)
  jsonlite::write_json(state, file.path(out, "workflow_state.json"), auto_unbox = TRUE, pretty = TRUE)
  write_run_manifest(config, "preprocess-and-cluster-scrna", out,
                     c(object_file, cell_file, file.path(out, "workflow_state.json")),
                     notes = paste0("Confirmed resolution ", resolution, " for scenario ", nm))
  quit(save = "no", status = 0)
}
if (action != "run") stop("workflow.action must be run or finalize_resolution")

seed <- as_int(cfg_get(config, "preprocessing.seed"), 1234L)
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
  rows <- lapply(resolutions, function(resolution) {
    reference <- cluster_graph(graph, resolution, random_seed)
    aris <- vapply(subsamples, function(cells) {
      subgraph <- graph[cells, cells, drop = FALSE]
      labels <- cluster_graph(subgraph, resolution, random_seed)
      adjusted_rand(reference[cells], labels[cells])
    }, numeric(1))
    data.frame(resolution = resolution, n_clusters = length(unique(reference)),
               stability_mean = mean(aris), stability_sd = stats::sd(aris),
               stability_conservative = mean(aris) - stats::sd(aris))
  })
  do.call(rbind, rows)
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
  neighbor_k <- as_int(cfg_get(scenario, "neighbors.k_param"), 20L)
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_use, dims = dims, k.param = neighbor_k,
                               graph.name = c(nn_name, snn_name), verbose = FALSE)
  umap_n_neighbors <- as_int(cfg_get(scenario, "umap.n_neighbors"), 30L)
  umap_min_dist <- as_num(cfg_get(scenario, "umap.min_dist"), 0.3)
  umap_metric <- cfg_get(scenario, "umap.metric", "cosine")
  umap_method <- cfg_get(scenario, "umap.method", "uwot")
  obj <- Seurat::RunUMAP(obj, reduction = reduction_use, dims = dims, reduction.name = umap_name,
                         reduction.key = paste0(key_stem, "UMAP_"), n.neighbors = umap_n_neighbors,
                         min.dist = umap_min_dist, metric = umap_metric, umap.method = umap_method,
                         seed.use = seed, verbose = FALSE)

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
  sizes <- as.data.frame(table(cluster = cluster), stringsAsFactors = FALSE)
  sizes$scenario <- nm
  sizes <- sizes[c("scenario", "cluster", "Freq")]
  names(sizes)[3] <- "cells"
  size_file <- file.path(out, paste0(nm, "_cluster_sizes.tsv"))
  utils::write.table(sizes, size_file, sep = "\t", quote = FALSE, row.names = FALSE)
  artifacts <- c(artifacts, size_file)

  composition <- as.data.frame(table(sample = obj[[sample_col, drop = TRUE]], cluster = cluster), stringsAsFactors = FALSE)
  composition$scenario <- nm
  composition <- composition[c("scenario", "sample", "cluster", "Freq")]
  names(composition)[4] <- "cells"
  comp_file <- file.path(out, paste0(nm, "_sample_cluster_counts.tsv"))
  utils::write.table(composition, comp_file, sep = "\t", quote = FALSE, row.names = FALSE)
  artifacts <- c(artifacts, comp_file)

  groups <- unique(c(cluster_name, plot_groups))
  p <- Seurat::DimPlot(obj, reduction = umap_name, group.by = groups, ncol = 1, pt.size = pt_size)
  plot_file <- file.path(out, paste0(nm, "_umap_diagnostics.png"))
  ggplot2::ggsave(plot_file, p, width = 7, height = max(5, 4 * length(groups)), units = "in", dpi = 300, bg = "white")
  artifacts <- c(artifacts, plot_file)
  }
  scenario_rows[[i]] <- data.frame(
    scenario = nm, regress_variables = paste(regress, collapse = ","), harmony = harmony_enabled,
    harmony_group_by = paste(harmony_groups, collapse = ","), theta = harmony_theta,
    reduction = reduction_use, dims = paste(range(dims), collapse = "-"), clustering_mode = clustering_mode,
    neighbor_k = neighbor_k, umap_n_neighbors = umap_n_neighbors,
    umap_min_dist = umap_min_dist, umap_metric = umap_metric, umap_method = umap_method,
    recommended_resolution = recommended_resolution,
    selected_resolution = if (final_available) selected_resolution else NA_real_,
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
save_scrna_object(obj, object_file)
cell_file <- file.path(out, "cell_assignments.tsv")
cluster_cols <- paste0(vapply(scenarios, function(x) sanitize_name(x$name), character(1)), "_clusters")
cluster_cols <- intersect(cluster_cols, colnames(obj[[]]))
cell_table <- data.frame(cell_id = colnames(obj), obj[[]][cluster_cols], check.names = FALSE)
utils::write.table(cell_table, cell_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"))

artifacts <- unique(c(artifacts, summary_file, similarity_file, elbow_file, object_file, cell_file, file.path(out, "session_info.txt")))
state_status <- if (any(summary$status == "awaiting_resolution_confirmation")) "awaiting_resolution_confirmation" else "complete"
state <- list(status = state_status, action = action, scenarios = scenario_rows,
              next_action = if (state_status == "awaiting_resolution_confirmation") "Review resolution plots and run finalize_resolution" else "none")
jsonlite::write_json(state, file.path(out, "workflow_state.json"), auto_unbox = TRUE, pretty = TRUE)
artifacts <- unique(c(artifacts, file.path(out, "workflow_state.json")))
write_run_manifest(config, "preprocess-and-cluster-scrna", out, artifacts,
                   notes = c(
                     paste0(length(scenarios), " user-selected scenario(s) executed"),
                     "No additional scenario was inserted or selected automatically"
                   ))
