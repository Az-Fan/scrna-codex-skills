#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(jsonlite); library(ggplot2); library(Seurat); library(patchwork) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript visualize_gene.R config.json")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
source(file.path(script_dir, "runtime.R")); source(file.path(script_dir, "figure_style.R"))
config <- read_skill_config(args[[1]]); out <- prepare_output(config)
if (is.null(config[["output"]])) config[["output"]] <- list()
config[["output"]][["figure_format"]] <- cfg_get(config, "plots.figure_format", "png")
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", tolower(x))
enabled <- function(name) isTRUE(cfg_get(config, paste0("plots.", name), TRUE))
figure_stems <- c("target_gene_featureplots", "target_gene_dotplot", "target_gene_violinplot", "target_gene_sample_expression", "target_gene_sample_heatmap", "target_gene_de_effect", "target_gene_volcano_highlight")
stale <- unlist(lapply(figure_stems, function(stem) file.path(out, paste0(stem, c(".png", ".pdf")))))
unlink(stale[file.exists(stale)])

gene_cfg <- cfg_get(config, "genes", required = TRUE)
gene_info <- do.call(rbind, lapply(gene_cfg, function(x) {
  if (is.character(x)) data.frame(symbol = x, label = x) else data.frame(symbol = x$symbol, label = x$label %||% x$symbol)
}))
gene_info$symbol <- as.character(gene_info$symbol); gene_info$label <- as.character(gene_info$label)
obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE))
assay <- cfg_get(config, "expression.assay", required = TRUE)
if (!assay %in% names(obj@assays)) stop("Expression assay not found: ", assay)
Seurat::DefaultAssay(obj) <- assay
layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character())
if (length(layers) && !"data" %in% layers) stop("A normalized data layer is required; this visualization skill does not normalize the object")

md <- obj[[]]; sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
condition_col <- cfg_get(config, "metadata.condition"); population_col <- cfg_get(config, "metadata.population")
reduction <- cfg_get(config, "metadata.reduction", "umap")
required_md <- unique(c(sample_col, condition_col, population_col)); required_md <- required_md[!is.na(required_md) & nzchar(required_md)]
assert_metadata(obj, required_md)
md$.sample <- as.character(md[[sample_col]])
md$.condition <- if (!is.null(condition_col) && nzchar(condition_col)) as.character(md[[condition_col]]) else "all"
md$.population <- if (!is.null(population_col) && nzchar(population_col)) as.character(md[[population_col]]) else md$.condition

available <- intersect(gene_info$symbol, rownames(obj)); missing <- setdiff(gene_info$symbol, available)
if (!length(available)) stop("None of the requested genes are present in assay ", assay)
expr <- Seurat::FetchData(obj, vars = available, layer = "data")
expr$cell_id <- rownames(expr)
meta <- data.frame(cell_id = rownames(md), sample = md$.sample, condition = md$.condition, population = md$.population, stringsAsFactors = FALSE)
cell <- merge(meta, expr, by = "cell_id", all.x = FALSE, sort = FALSE)
long <- do.call(rbind, lapply(available, function(g) data.frame(cell[, c("cell_id", "sample", "condition", "population")], symbol = g, label = gene_info$label[match(g, gene_info$symbol)], expression = as.numeric(cell[[g]]), stringsAsFactors = FALSE)))

keys <- interaction(long$symbol, long$label, long$sample, long$condition, long$population, drop = TRUE, lex.order = TRUE)
summary_rows <- lapply(split(long, keys), function(x) data.frame(symbol = x$symbol[[1]], label = x$label[[1]], sample = x$sample[[1]], condition = x$condition[[1]], population = x$population[[1]], n_cells = nrow(x), mean_expression = mean(x$expression), median_expression = median(x$expression), detected_fraction = mean(x$expression > 0), stringsAsFactors = FALSE))
cell_summary <- do.call(rbind, summary_rows); rownames(cell_summary) <- NULL
write_tsv(cell_summary, file.path(out, "target_gene_cell_expression_summary.tsv"))

de <- NULL; de_path <- cfg_get(config, "input.differential_table")
de_gene <- cfg_get(config, "differential_columns.gene", "gene"); de_lfc <- cfg_get(config, "differential_columns.log2_fold_change", "log2FoldChange"); de_padj <- cfg_get(config, "differential_columns.adjusted_p_value", "padj")
if (!is.null(de_path) && nzchar(de_path)) {
  de <- if (grepl("\\.csv$", de_path, ignore.case = TRUE)) read.csv(de_path, check.names = FALSE) else read.delim(de_path, check.names = FALSE)
  miss_cols <- setdiff(c(de_gene, de_lfc, de_padj), names(de)); if (length(miss_cols)) stop("Missing differential table columns: ", paste(miss_cols, collapse = ", "))
  de[[de_gene]] <- as.character(de[[de_gene]])
}

pb <- NULL; pb_path <- cfg_get(config, "input.pseudobulk_data")
if (!is.null(pb_path) && nzchar(pb_path)) {
  z <- readRDS(pb_path)
  if (is.list(z) && !is.null(z$normalized_counts) && !is.null(z$coldata)) pb <- z else warning("pseudobulk_data lacks normalized_counts and coldata; falling back to sample means")
}

sample_values <- NULL
if (!is.null(pb)) {
  mat <- as.matrix(pb$normalized_counts); common <- intersect(available, rownames(mat))
  if (length(common)) {
    cd <- as.data.frame(pb$coldata); ids <- colnames(mat)
    pb_sample_col <- if (sample_col %in% names(cd)) sample_col else if ("sample" %in% names(cd)) "sample" else NULL
    pb_condition_col <- if (!is.null(condition_col) && condition_col %in% names(cd)) condition_col else if ("condition" %in% names(cd)) "condition" else NULL
    sample_values <- do.call(rbind, lapply(common, function(g) data.frame(symbol = g, label = gene_info$label[match(g, gene_info$symbol)], sample = if (is.null(pb_sample_col)) ids else as.character(cd[[pb_sample_col]]), condition = if (is.null(pb_condition_col)) "all" else as.character(cd[[pb_condition_col]]), value = log2(as.numeric(mat[g, ]) + 1), source = "pseudobulk_normalized_counts_log2p1", stringsAsFactors = FALSE)))
  }
}
if (is.null(sample_values)) {
  sample_values <- aggregate(expression ~ symbol + label + sample + condition, long, mean)
  names(sample_values)[names(sample_values) == "expression"] <- "value"; sample_values$source <- "mean_normalized_cell_expression"
}
write_tsv(sample_values, file.path(out, "target_gene_sample_expression.tsv"))

gene_status <- gene_info
gene_status$found_in_object <- gene_status$symbol %in% available
gene_status$found_in_de <- if (is.null(de)) NA else gene_status$symbol %in% de[[de_gene]]
gene_status$found_in_pseudobulk <- if (is.null(pb)) NA else gene_status$symbol %in% rownames(pb$normalized_counts)
write_tsv(gene_status, file.path(out, "gene_status.tsv"))

target_summary <- gene_status
target_summary$log2_fold_change <- NA_real_; target_summary$adjusted_p_value <- NA_real_
if (!is.null(de)) { m <- match(target_summary$symbol, de[[de_gene]]); target_summary$log2_fold_change <- as.numeric(de[[de_lfc]][m]); target_summary$adjusted_p_value <- as.numeric(de[[de_padj]][m]) }
write_tsv(target_summary, file.path(out, "target_gene_summary.tsv"))

status <- data.frame(family = character(), file = character(), status = character(), reason = character(), stringsAsFactors = FALSE)
record <- function(family, files = character(), state = "generated", reason = "") { if (!length(files)) files <- ""; status <<- rbind(status, data.frame(family = family, file = files, status = state, reason = reason, stringsAsFactors = FALSE)) }
artifacts <- c(file.path(out, c("gene_status.tsv", "target_gene_cell_expression_summary.tsv", "target_gene_sample_expression.tsv", "target_gene_summary.tsv")))
condition_colors <- paper_colors(sample_values$condition, "condition", config, out)

if (enabled("featureplot") && reduction %in% names(obj@reductions)) {
  pages <- split(available, ceiling(seq_along(available) / 6)); files <- character()
  for (i in seq_along(pages)) {
    ps <- Seurat::FeaturePlot(obj, features = pages[[i]], reduction = reduction, cols = c("#ECECEC", "#B2182B"), order = TRUE, raster = FALSE, combine = FALSE, keep.scale = "feature")
    ps <- lapply(seq_along(ps), function(j) ps[[j]] + paper_theme() + Seurat::NoAxes() + labs(title = gene_info$label[match(pages[[i]][j], gene_info$symbol)], subtitle = paste(assay, "normalized expression; per-gene scale")))
    stem <- file.path(out, if (length(pages) == 1) "target_gene_featureplots" else sprintf("target_gene_featureplots_page%02d", i))
    files <- c(files, paper_save(wrap_plots(ps, ncol = min(3, length(ps))), stem, 5 * min(3, length(ps)), 4.8 * ceiling(length(ps) / 3), config, out, "target_gene_featureplots"))
  }; record("target_gene_featureplots", files); artifacts <- c(artifacts, files)
} else record("target_gene_featureplots", state = "skipped", reason = if (!enabled("featureplot")) "disabled" else paste("reduction unavailable:", reduction))

if (enabled("dotplot") && length(unique(md$.population)) >= 2) {
  obj$.gene_plot_group <- md$.population
  files <- paper_dotplot(obj, available, assay, ".gene_plot_group", file.path(out, "target_gene_dotplot"), config, out)
  record("target_gene_dotplot", files); artifacts <- c(artifacts, files)
} else record("target_gene_dotplot", state = "skipped", reason = if (!enabled("dotplot")) "disabled" else "population field has fewer than two levels")

if (enabled("violin") && length(unique(long$condition)) >= 2) {
  p <- ggplot(long, aes(condition, expression, fill = condition)) + geom_violin(scale = "width", trim = TRUE, linewidth = .25) + geom_boxplot(width = .12, outlier.shape = NA, fill = "white", linewidth = .25) + facet_wrap(~ label, scales = "free_y") + scale_fill_manual(values = condition_colors) + paper_theme() + labs(title = "Target-gene expression by condition", subtitle = "Cell-level distributions are descriptive; cells are not biological replicates", x = NULL, y = "Normalized expression", fill = "Condition")
  files <- paper_save(p, file.path(out, "target_gene_violinplot"), 4 + min(4, length(available)) * 2, 4.8 * ceiling(length(available) / 4), config, out, "target_gene_violinplot"); record("target_gene_violinplot", files); artifacts <- c(artifacts, files)
} else record("target_gene_violinplot", state = "skipped", reason = if (!enabled("violin")) "disabled" else "condition field has fewer than two levels")

if (enabled("sample_expression") && length(unique(sample_values$sample)) >= 2) {
  n_sample_genes <- length(unique(sample_values$symbol)); sample_ncol <- min(3, n_sample_genes)
  p <- ggplot(sample_values, aes(condition, value, colour = condition)) + geom_jitter(width = .12, height = 0, size = 2.4, alpha = .9) + stat_summary(fun = median, geom = "point", shape = 95, size = 7, colour = "#303030") + facet_wrap(~ label, scales = "free_y", ncol = sample_ncol) + scale_colour_manual(values = condition_colors) + paper_theme() + theme(legend.position = "bottom") + labs(title = "Target-gene sample expression", subtitle = paste("Each point is a biological sample; source:", unique(sample_values$source)[1]), x = NULL, y = "Expression", colour = "Condition")
  files <- paper_save(p, file.path(out, "target_gene_sample_expression"), max(8, sample_ncol * 3.2), 3.6 * ceiling(n_sample_genes / sample_ncol) + 1.2, config, out, "target_gene_sample_expression"); record("target_gene_sample_expression", files); artifacts <- c(artifacts, files)
} else record("target_gene_sample_expression", state = "skipped", reason = if (!enabled("sample_expression")) "disabled" else "fewer than two biological samples")

if (enabled("sample_heatmap") && length(unique(sample_values$symbol)) >= 2 && length(unique(sample_values$sample)) >= 2) {
  hm <- sample_values; hm$z <- ave(hm$value, hm$symbol, FUN = function(x) { z <- as.numeric(scale(x)); ifelse(is.finite(z), z, 0) })
  p <- ggplot(hm, aes(sample, label, fill = z)) + geom_tile(colour = "white", linewidth = .25) + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) + paper_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Target-gene sample heatmap", subtitle = "Gene-wise z scores across biological samples", x = NULL, y = NULL, fill = "Row z score")
  files <- paper_save(p, file.path(out, "target_gene_sample_heatmap"), max(7, length(unique(hm$sample)) * .42), max(4, length(unique(hm$symbol)) * .35 + 2.5), config, out, "target_gene_sample_heatmap"); record("target_gene_sample_heatmap", files); artifacts <- c(artifacts, files)
} else record("target_gene_sample_heatmap", state = "skipped", reason = if (!enabled("sample_heatmap")) "disabled" else "requires at least two genes and two samples")

de_targets <- target_summary[is.finite(target_summary$log2_fold_change), , drop = FALSE]
min_effect <- as.numeric(cfg_get(config, "plots.min_abs_effect_to_plot", 0.001))
if (enabled("de_effect") && nrow(de_targets) && max(abs(de_targets$log2_fold_change), na.rm = TRUE) >= min_effect) {
  de_targets$label <- factor(de_targets$label, levels = rev(de_targets$label[order(de_targets$log2_fold_change)])); de_targets$evidence <- pmin(-log10(pmax(de_targets$adjusted_p_value, .Machine$double.xmin)), 20)
  p <- ggplot(de_targets, aes(log2_fold_change, label)) + geom_vline(xintercept = 0, colour = "#777777", linewidth = .35) + geom_segment(aes(x = 0, xend = log2_fold_change, yend = label), colour = "#BDBDBD") + geom_point(aes(size = evidence, colour = log2_fold_change > 0)) + scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#4477AA"), guide = "none") + scale_size_continuous(name = expression(-log[10](adjusted~P)), range = c(2.5, 6)) + paper_theme() + labs(title = "Target-gene differential effects", subtitle = "Values are imported from the supplied complete differential table", x = expression(log[2]~fold~change), y = NULL)
  files <- paper_save(p, file.path(out, "target_gene_de_effect"), 7, max(4, nrow(de_targets) * .4 + 2), config, out, "target_gene_de_effect"); record("target_gene_de_effect", files); artifacts <- c(artifacts, files)
} else record("target_gene_de_effect", state = "skipped", reason = if (!enabled("de_effect")) "disabled" else if (!nrow(de_targets)) "no finite target-gene differential effects" else paste0("all absolute target-gene effects are below ", min_effect))

if (enabled("volcano_highlight") && !is.null(de)) {
  vd <- data.frame(symbol = de[[de_gene]], lfc = as.numeric(de[[de_lfc]]), padj = as.numeric(de[[de_padj]]), stringsAsFactors = FALSE); vd <- vd[is.finite(vd$lfc) & is.finite(vd$padj), ]; vd$y <- -log10(pmax(vd$padj, .Machine$double.xmin)); vd$target <- vd$symbol %in% available
  if (nrow(vd) >= 3) {
    p <- ggplot(vd, aes(lfc, y)) + geom_point(colour = "#C8C8C8", size = .8, alpha = .55) + geom_point(data = vd[vd$target, ], colour = "#D55E00", size = 2.2) + geom_vline(xintercept = 0, colour = "#777777", linewidth = .3) + paper_theme() + labs(title = "Target genes in the complete differential result", subtitle = "Grey points are all finite tested genes; selected genes are highlighted", x = expression(log[2]~fold~change), y = expression(-log[10](adjusted~P)))
    labels <- transform(vd[vd$target, ], display_label = gene_info$label[match(symbol, gene_info$symbol)])
    if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(data = labels, aes(label = display_label), colour = "#8C2D04", size = 3.2, max.overlaps = Inf, seed = 1)
    else p <- p + geom_text(data = labels, aes(label = display_label), colour = "#8C2D04", vjust = -0.8, check_overlap = FALSE)
    files <- paper_save(p, file.path(out, "target_gene_volcano_highlight"), 7, 5.5, config, out, "target_gene_volcano_highlight"); record("target_gene_volcano_highlight", files); artifacts <- c(artifacts, files)
  } else record("target_gene_volcano_highlight", state = "skipped", reason = "fewer than three finite DE rows")
} else record("target_gene_volcano_highlight", state = "skipped", reason = if (!enabled("volcano_highlight")) "disabled" else "no differential table supplied")

for (g in missing) paper_record(out, "requested_gene", g, "skipped", "gene absent from expression assay")
write_tsv(status, file.path(out, "plot_status.tsv")); artifacts <- c(artifacts, file.path(out, "plot_status.tsv"))
writeLines(capture.output(sessionInfo()), technical_path(out, "session_info.txt")); artifacts <- c(artifacts, technical_path(out, "session_info.txt"))
write_run_manifest(config, "15-scrna-visualize-gene", out, artifacts, notes = c("paper_v1 fixed target-gene figure set", paste("sample expression source:", unique(sample_values$source)[1])))
