args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript annotate_cells.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
source(file.path(dirname(normalizePath(script_file)), "figure_style.R"))
config <- read_skill_config(args[[1]])
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")

action <- cfg_get(config, "workflow.action", required = TRUE)
obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE), "auto")
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
cluster_col <- cfg_get(config, "metadata.cluster")
condition_col <- cfg_get(config, "metadata.condition")
reduction <- cfg_get(config, "metadata.reduction", "umap")
assert_metadata(obj, c(sample_col, condition_col %||% character()))
out <- prepare_output(config)

write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cluster_levels <- function(x) {
  values <- unique(as.character(x))
  numeric_values <- suppressWarnings(as.numeric(values))
  if (!anyNA(numeric_values)) values[order(numeric_values)] else sort(values)
}

if (identical(action, "prepare_review")) {
  markers_path <- cfg_get(config, "input.markers")
  compute_clustering <- isTRUE(cfg_get(config, "clustering.compute_if_missing", FALSE))
  assay <- cfg_get(config, "markers.assay", Seurat::DefaultAssay(obj))
  Seurat::DefaultAssay(obj) <- assay
  layers <- if (utils::packageVersion("SeuratObject") >= "5.0.0") SeuratObject::Layers(obj[[assay]]) else c("counts", "data")
  if (!any(layers == "data" | grepl("^data\\.", layers))) obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)

  if (!is.null(cluster_col) && nzchar(cluster_col)) {
    assert_metadata(obj, cluster_col)
    clusters <- as.character(obj[[]][[cluster_col]])
    if (anyNA(clusters)) stop("Cluster column contains missing values: ", cluster_col)
    Seurat::Idents(obj) <- factor(clusters, levels = cluster_levels(clusters))
  } else if (compute_clustering) {
    dims <- unlist(cfg_get(config, "clustering.dims", as.list(1:30)), use.names = FALSE)
    resolution <- as.numeric(cfg_get(config, "clustering.resolution", 0.5))
    if (!"pca" %in% names(obj@reductions)) obj <- Seurat::FindVariableFeatures(obj, assay = assay, verbose = FALSE) |> Seurat::ScaleData(assay = assay, verbose = FALSE) |> Seurat::RunPCA(assay = assay, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE) |> Seurat::FindClusters(resolution = resolution, verbose = FALSE)
    if (!reduction %in% names(obj@reductions)) obj <- Seurat::RunUMAP(obj, dims = dims, reduction.name = reduction, verbose = FALSE)
    cluster_col <- "seurat_clusters"
  } else {
    stop("prepare_review requires metadata.cluster or clustering.compute_if_missing=true")
  }

  if (!reduction %in% names(obj@reductions)) stop("Reduction not found: ", reduction)
  clusters <- as.character(Seurat::Idents(obj))
  if (!is.null(markers_path)) {
    markers_path <- normalizePath(markers_path, mustWork = TRUE)
    markers <- utils::read.delim(markers_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("cluster", "gene") %in% names(markers))) stop("Existing marker table must contain cluster and gene columns")
    marker_note <- paste0("markers_reused_from=", markers_path)
  } else {
    markers <- Seurat::FindAllMarkers(
      obj, assay = assay, only.pos = isTRUE(cfg_get(config, "markers.only_pos", TRUE)),
      test.use = cfg_get(config, "markers.test_use", "wilcox"),
      logfc.threshold = as.numeric(cfg_get(config, "markers.logfc_threshold", 0.25)),
      min.pct = as.numeric(cfg_get(config, "markers.min_pct", 0.1)), verbose = FALSE
    )
    if (!nrow(markers)) stop("No cluster markers were identified")
    if (!"gene" %in% names(markers)) markers$gene <- rownames(markers)
    marker_note <- "markers_computed_in_annotation_prepare_review"
  }

  marker_file <- file.path(out, "cluster_markers.tsv"); write_tsv(markers, marker_file)
  cluster_ids <- cluster_levels(clusters)
  review <- data.frame(
    cluster = cluster_ids, candidate_broad = "", candidate_fine = "", candidate_state = "",
    evidence = "", conflicts = "", sample_bias = "", qc_flag = "", confidence = "",
    decision = "pending", notes = "", stringsAsFactors = FALSE
  )
  review_file <- file.path(out, "annotation_review.tsv"); write_tsv(review, review_file)
  object_file <- file.path(out, "clustered_object.qs"); save_scrna_object(obj, object_file)

  plot_file <- paper_dimplot(obj, reduction, cluster_col, file.path(out, "cluster_umap"), config, out)
  sample_plot_file <- paper_dimplot(obj, reduction, c(cluster_col, sample_col), file.path(out, "cluster_sample_umap"), config, out)
  feature_files <- paper_features(obj, reduction, config, out)

  canonical <- cfg_get(config, "markers.canonical", list())
  canonical_file <- NULL
  if (length(canonical)) {
    genes <- unique(unlist(canonical, use.names = FALSE)); genes <- genes[genes %in% rownames(obj)]
    if (length(genes)) {
      canonical_file <- paper_dotplot(obj, genes, assay, cluster_col, file.path(out, "canonical_marker_dotplot"), config, out)
    }
  }
  write_run_manifest(config, "08-scrna-annotate-cells", out,
    c(marker_file, review_file, object_file, plot_file, sample_plot_file, canonical_file %||% character(), feature_files),
    c("action=prepare_review", marker_note, "Formal labels were not written"))
  quit(save = "no", status = 0L)
}

if (!identical(action, "apply_confirmed")) stop("Unknown workflow.action: ", action)
if (is.null(cluster_col) || !nzchar(cluster_col)) stop("apply_confirmed requires metadata.cluster")
assert_metadata(obj, cluster_col)
if (!reduction %in% names(obj@reductions)) stop("Reduction not found: ", reduction)

decision_path <- normalizePath(cfg_get(config, "input.decisions", required = TRUE), mustWork = TRUE)
decisions <- utils::read.delim(decision_path, check.names = FALSE, stringsAsFactors = FALSE)
broad_col <- cfg_get(config, "annotation.broad_column", required = TRUE)
fine_col <- cfg_get(config, "annotation.fine_column", required = TRUE)
state_col <- cfg_get(config, "annotation.state_column")
decision_col <- cfg_get(config, "annotation.decision_column", "decision")
required <- unique(c("cluster", broad_col, fine_col, state_col %||% character(), decision_col))
missing <- setdiff(required, names(decisions)); if (length(missing)) stop("Annotation decisions are missing columns: ", paste(missing, collapse = ", "))
decisions$cluster <- as.character(decisions$cluster)
if (anyDuplicated(decisions$cluster)) stop("Annotation decisions contain duplicated clusters")
clusters <- as.character(obj[[]][[cluster_col]])
object_clusters <- cluster_levels(clusters)
if (!setequal(object_clusters, decisions$cluster)) stop("Annotation decisions must cover exactly every object cluster")
confirmed_values <- tolower(unlist(cfg_get(config, "annotation.confirmed_values", list("confirmed")), use.names = FALSE))
if (any(!tolower(as.character(decisions[[decision_col]])) %in% confirmed_values)) stop("Every annotation decision must be confirmed before apply_confirmed")
if (any(!nzchar(as.character(decisions[[broad_col]]))) || any(!nzchar(as.character(decisions[[fine_col]])))) stop("Confirmed broad and fine labels cannot be empty")

broad_out <- cfg_get(config, "annotation.output_broad_column", broad_col)
fine_out <- cfg_get(config, "annotation.output_fine_column", fine_col)
state_out <- if (!is.null(state_col)) cfg_get(config, "annotation.output_state_column", state_col) else NULL
broad_map <- stats::setNames(as.character(decisions[[broad_col]]), decisions$cluster)
fine_map <- stats::setNames(as.character(decisions[[fine_col]]), decisions$cluster)
annotation_meta <- data.frame(row.names = colnames(obj))
annotation_meta[[broad_out]] <- factor(broad_map[clusters]); annotation_meta[[fine_out]] <- factor(fine_map[clusters])
if (!is.null(state_col)) {
  state_map <- stats::setNames(as.character(decisions[[state_col]]), decisions$cluster)
  annotation_meta[[state_out]] <- factor(state_map[clusters])
}
obj <- Seurat::AddMetaData(obj, metadata = annotation_meta)

cell_table <- data.frame(cell_id = colnames(obj), cluster = clusters, sample = as.character(obj[[]][[sample_col]]), stringsAsFactors = FALSE)
if (!is.null(condition_col)) cell_table$condition <- as.character(obj[[]][[condition_col]])
cell_table[[broad_out]] <- as.character(obj[[]][[broad_out]]); cell_table[[fine_out]] <- as.character(obj[[]][[fine_out]])
if (!is.null(state_out)) cell_table[[state_out]] <- as.character(obj[[]][[state_out]])
cell_file <- file.path(out, "cell_annotations.tsv"); write_tsv(cell_table, cell_file)

summary <- decisions[match(object_clusters, decisions$cluster), , drop = FALSE]
summary$n_cells <- as.integer(table(factor(clusters, levels = object_clusters)))
summary_file <- file.path(out, "annotation_summary.tsv"); write_tsv(summary, summary_file)
object_file <- file.path(out, "annotated_object.qs"); save_scrna_object(obj, object_file)

annotation_plot_file <- paper_dimplot(obj, reduction, c(broad_out, fine_out), file.path(out, "annotated_umap"), config, out)
audit_groups <- c(cluster_col, sample_col, condition_col %||% character())
audit_plot_file <- paper_dimplot(obj, reduction, audit_groups, file.path(out, "cluster_sample_condition_umap"), config, out)
feature_files <- paper_features(obj, reduction, config, out)
session_file <- technical_path(out, "session_info.txt"); writeLines(capture.output(utils::sessionInfo()), session_file)
write_run_manifest(config, "08-scrna-annotate-cells", out,
  c(object_file, cell_file, summary_file, annotation_plot_file, audit_plot_file, feature_files, session_file),
  c("action=apply_confirmed", paste0("decisions=", decision_path), "Cluster IDs remain separate from biological labels"))
