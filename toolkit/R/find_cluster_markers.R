args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript find_cluster_markers.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
config <- read_skill_config(args[[1]])
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")

obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE), "auto")
assay <- cfg_get(config, "analysis.assay", Seurat::DefaultAssay(obj))
if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
Seurat::DefaultAssay(obj) <- assay

cluster_col <- cfg_get(config, "metadata.cluster")
if (!is.null(cluster_col) && nzchar(cluster_col)) {
  assert_metadata(obj, cluster_col)
  cluster_values <- obj[[]][[cluster_col]]
  if (anyNA(cluster_values)) stop("Cluster column contains missing values: ", cluster_col)
  Seurat::Idents(obj) <- as.factor(cluster_values)
}
idents <- droplevels(Seurat::Idents(obj))
if (nlevels(idents) < 2L) stop("FindAllMarkers requires at least two identity groups")
cluster_counts <- as.data.frame(table(cluster = idents), stringsAsFactors = FALSE)
names(cluster_counts)[2] <- "n_cells"

if (isTRUE(cfg_get(config, "analysis.join_layers", TRUE)) &&
    requireNamespace("SeuratObject", quietly = TRUE) &&
    utils::packageVersion("SeuratObject") >= "5.0.0") {
  layers <- SeuratObject::Layers(obj[[assay]])
  if (sum(grepl("^data\\.", layers)) > 1L || sum(grepl("^counts\\.", layers)) > 1L) {
    obj <- Seurat::JoinLayers(obj, assay = assay)
  }
}
normalized_in_memory <- FALSE
if (requireNamespace("SeuratObject", quietly = TRUE) &&
    utils::packageVersion("SeuratObject") >= "5.0.0") {
  layers <- SeuratObject::Layers(obj[[assay]])
  has_data <- any(layers == "data" | grepl("^data\\.", layers))
} else {
  has_data <- nrow(Seurat::GetAssayData(obj, assay = assay, slot = "data")) > 0L
}
if (!has_data) {
  if (!isTRUE(cfg_get(config, "analysis.normalize_if_missing", TRUE))) {
    stop("Selected assay has no normalized data layer; normalize it first or set analysis.normalize_if_missing=true")
  }
  obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
  normalized_in_memory <- TRUE
}

max_cells <- cfg_get(config, "analysis.max_cells_per_ident")
if (is.null(max_cells)) max_cells <- Inf
set.seed(as.integer(cfg_get(config, "analysis.random_seed", 1)))
markers <- Seurat::FindAllMarkers(
  object = obj,
  assay = assay,
  only.pos = isTRUE(cfg_get(config, "analysis.only_pos", TRUE)),
  test.use = cfg_get(config, "analysis.test_use", "wilcox"),
  logfc.threshold = as.numeric(cfg_get(config, "analysis.logfc_threshold", 0.25)),
  min.pct = as.numeric(cfg_get(config, "analysis.min_pct", 0.1)),
  min.diff.pct = as.numeric(cfg_get(config, "analysis.min_diff_pct", -Inf)),
  return.thresh = as.numeric(cfg_get(config, "analysis.return_thresh", 0.05)),
  max.cells.per.ident = as.numeric(max_cells),
  random.seed = as.integer(cfg_get(config, "analysis.random_seed", 1)),
  verbose = FALSE
)
if (!nrow(markers)) stop("No markers passed the configured thresholds; review assay, identities, layers, and thresholds")
if (!"gene" %in% names(markers)) markers$gene <- rownames(markers)
fc_col <- intersect(c("avg_log2FC", "avg_logFC", "avg_diff"), names(markers))[1]
if (is.na(fc_col)) stop("FindAllMarkers result has no recognized effect-size column")
p_col <- if ("p_val_adj" %in% names(markers)) "p_val_adj" else "p_val"
ordering <- order(as.character(markers$cluster), markers[[p_col]], -markers[[fc_col]], markers$gene)
markers <- markers[ordering, , drop = FALSE]
markers$rank_within_cluster <- ave(seq_len(nrow(markers)), as.character(markers$cluster), FUN = seq_along)
top_n <- as.integer(cfg_get(config, "reporting.top_n", 20))
if (is.na(top_n) || top_n < 1L) stop("reporting.top_n must be a positive integer")
top_markers <- markers[markers$rank_within_cluster <= top_n, , drop = FALSE]

marker_counts <- as.data.frame(table(cluster = factor(as.character(markers$cluster), levels = as.character(cluster_counts$cluster))))
names(marker_counts)[2] <- "n_markers"
summary <- merge(cluster_counts, marker_counts, by = "cluster", all.x = TRUE, sort = FALSE)
summary$n_markers[is.na(summary$n_markers)] <- 0L

out <- prepare_output(config)
marker_path <- file.path(out, "cluster_markers.tsv")
top_path <- file.path(out, "top_cluster_markers.tsv")
summary_path <- file.path(out, "cluster_marker_summary.tsv")
plot_path <- file.path(out, "top_marker_dotplot.pdf")
utils::write.table(markers, marker_path, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(top_markers, top_path, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

dot_n <- as.integer(cfg_get(config, "reporting.dotplot_top_n", 5))
if (is.na(dot_n) || dot_n < 1L) stop("reporting.dotplot_top_n must be a positive integer")
dot_genes <- unique(top_markers$gene[top_markers$rank_within_cluster <= dot_n])
grDevices::pdf(plot_path, width = max(8, min(20, length(dot_genes) * 0.22 + 4)), height = max(5, nrow(cluster_counts) * 0.28 + 3))
print(Seurat::DotPlot(obj, features = dot_genes, assay = assay) + Seurat::RotatedAxis())
grDevices::dev.off()

write_run_manifest(
  config, "06-scrna-find-cluster-markers", out,
  c(marker_path, top_path, summary_path, plot_path),
  c(paste0("assay=", assay), paste0("grouping=", if (is.null(cluster_col)) "active identities" else cluster_col),
    paste0("normalized_in_memory=", normalized_in_memory), "Input object was not rewritten")
)
