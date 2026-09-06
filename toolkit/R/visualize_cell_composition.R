#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite); library(ggplot2); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript visualize_cell_composition.R config.json")
cfg <- fromJSON(args[[1]], simplifyVector = FALSE)
getv <- function(x, path, default = NULL) { for (p in strsplit(path, "\\.")[[1]]) { if (is.null(x) || is.null(x[[p]])) return(default); x <- x[[p]] }; x }
out <- normalizePath(getv(cfg, "output_dir"), mustWork = FALSE); dir.create(out, recursive = TRUE, showWarnings = FALSE)
tech <- file.path(out, "_provenance"); dir.create(tech, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", tolower(x))
fmt <- tolower(getv(cfg, "plots.figure_format", "png")); if (fmt == "both") formats <- c("png", "pdf") else formats <- fmt
if (any(!formats %in% c("png", "pdf"))) stop("plots.figure_format must be png, pdf or both")
theme_paper <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9), strip.background = element_rect(fill = "#F4F4F4", colour = NA), strip.text = element_text(face = "bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 9), plot.margin = margin(8, 12, 8, 8))
save_plot <- function(p, stem, width = 9, height = 6, family = basename(stem)) {
  files <- character(); for (f in formats) { path <- paste0(stem, ".", f); ggsave(path, p, width = width, height = height, units = "in", dpi = 300, bg = "white"); files <- c(files, path) }
  files
}
status <- data.frame(family = character(), file = character(), status = character(), reason = character(), stringsAsFactors = FALSE)
record <- function(family, file, state = "generated", reason = "") status <<- rbind(status, data.frame(family = family, file = file, status = state, reason = reason, stringsAsFactors = FALSE))

input_object <- getv(cfg, "input.object")
counts_path <- getv(cfg, "input.counts_table")
if ((is.null(input_object) || !nzchar(input_object)) == (is.null(counts_path) || !nzchar(counts_path))) stop("Provide exactly one of input.object or input.counts_table")
sample_col <- getv(cfg, "metadata.sample"); condition_col <- getv(cfg, "metadata.condition"); batch_col <- getv(cfg, "metadata.batch"); reduction <- getv(cfg, "metadata.reduction", "umap")
vars <- getv(cfg, "composition.variables", list())
if (!length(vars)) stop("composition.variables must contain at least one variable")
`%||%` <- function(a, b) if (is.null(a) || !length(a) || !nzchar(as.character(a))) b else a
var_info <- lapply(vars, function(v) list(column = v$column, label = v$label %||% v$column, kind = v$kind %||% "metadata"))
denom_mode <- getv(cfg, "composition.denominator.mode", "all_input_cells"); denom_desc <- getv(cfg, "composition.denominator.description")
if (is.null(denom_desc) || !nzchar(denom_desc)) stop("composition.denominator.description is required")

obj <- NULL; md <- NULL
if (!is.null(input_object) && nzchar(input_object)) {
  suppressPackageStartupMessages(library(Seurat))
  if (grepl("\\.qs$", input_object, ignore.case = TRUE)) { suppressPackageStartupMessages(library(qs)); obj <- qread(input_object) } else obj <- readRDS(input_object)
  md <- obj[[]]
} else {
  md <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
}
required <- c(sample_col, condition_col, batch_col, vapply(var_info, `[[`, character(1), "column")); required <- unique(required[nzchar(required)])
missing <- setdiff(required, names(md)); if (length(missing)) stop("Missing metadata/count columns: ", paste(missing, collapse = ", "))
if (is.null(condition_col) || !nzchar(condition_col)) md$.condition <- "all" else md$.condition <- as.character(md[[condition_col]])
if (is.null(batch_col) || !nzchar(batch_col)) md$.batch <- "not_configured" else md$.batch <- as.character(md[[batch_col]])
md$.sample <- as.character(md[[sample_col]])
count_col <- getv(cfg, "input.count_column", "n_cells")
parent_col <- getv(cfg, "composition.parent_column")
if (!is.null(parent_col) && nzchar(parent_col) && !parent_col %in% names(md)) stop("Missing parent column: ", parent_col)

make_long <- function(v) {
  cat <- as.character(md[[v$column]]); cat[is.na(cat) | !nzchar(cat)] <- "Missing"
  parent <- if (!is.null(parent_col) && nzchar(parent_col)) as.character(md[[parent_col]]) else "all"
  parent[is.na(parent) | !nzchar(parent)] <- "Missing"
  if (denom_mode %in% c("selected_parent", "selected_cell_types")) {
    include <- unlist(getv(cfg, "composition.denominator.include", character()))
    if (length(include)) keep <- parent %in% include else keep <- rep(TRUE, length(parent))
  } else keep <- rep(TRUE, length(parent))
  if (is.null(obj)) {
    n <- as.numeric(md[[count_col]]); if (anyNA(n) || any(n < 0)) stop("Count column must be non-negative numeric")
  } else n <- rep(1, nrow(md))
  d <- data.frame(sample = md$.sample, condition = md$.condition, batch = md$.batch, category = cat, parent = parent, n_cells = n, stringsAsFactors = FALSE)
  d <- d[keep, , drop = FALSE]
  agg <- aggregate(n_cells ~ sample + condition + batch + parent + category, d, sum)
  samples <- unique(d[c("sample", "condition", "batch", "parent")]); cats <- unique(agg$category)
  grid <- merge(samples, data.frame(category = cats, stringsAsFactors = FALSE), all = TRUE)
  grid <- merge(grid, agg, by = c("sample", "condition", "batch", "parent", "category"), all.x = TRUE); grid$n_cells[is.na(grid$n_cells)] <- 0
  den <- aggregate(n_cells ~ sample + parent, grid, sum); names(den)[3] <- "denominator_cells"
  grid <- merge(grid, den, by = c("sample", "parent"), all.x = TRUE); grid$proportion <- ifelse(grid$denominator_cells > 0, grid$n_cells / grid$denominator_cells, NA_real_)
  grid$variable <- v$column; grid$variable_label <- v$label; grid$kind <- v$kind
  grid[order(grid$sample, grid$parent, grid$category), c("variable", "variable_label", "kind", "sample", "condition", "batch", "parent", "category", "n_cells", "denominator_cells", "proportion")]
}

all_long <- do.call(rbind, lapply(var_info, make_long)); write_tsv(all_long, file.path(out, "composition_counts.tsv")); write_tsv(all_long, file.path(out, "composition_proportions.tsv"))
coverage <- aggregate(cbind(n_cells, denominator_cells) ~ variable + sample + condition + batch + parent, all_long, function(x) x[[1]])
coverage$n_categories_observed <- vapply(seq_len(nrow(coverage)), function(i) sum(all_long$variable == coverage$variable[i] & all_long$sample == coverage$sample[i] & all_long$parent == coverage$parent[i] & all_long$n_cells > 0), integer(1))
coverage$low_coverage <- coverage$denominator_cells < as.numeric(getv(cfg, "quality.min_cells_per_sample", 20))
write_tsv(coverage, file.path(out, "sample_coverage.tsv"))

audit <- do.call(rbind, lapply(var_info, function(v) {
  x <- all_long[all_long$variable == v$column, ]; data.frame(variable = v$column, label = v$label, kind = v$kind, n_samples = length(unique(x$sample)), n_categories = length(unique(x$category)), condition_levels = length(unique(x$condition)), batch_levels = length(unique(x$batch)), denominator_mode = denom_mode, denominator_description = denom_desc, stringsAsFactors = FALSE)
}))
audit$batch_diagnostic <- ifelse(audit$batch_levels < 2, "unavailable_or_single_level", ifelse(audit$condition_levels < 2, "composition_only_no_condition_comparison", "available_descriptive_diagnostic")); write_tsv(audit, file.path(out, "composition_audit.tsv"))

colmap <- function(values, field) { configured <- getv(cfg, paste0("figure_style.colors.", field), list()); labels <- sort(unique(as.character(values))); pal <- grDevices::hcl.colors(length(labels), "Dynamic"); names(pal) <- labels; if (length(configured)) { z <- unlist(configured); pal[names(z)[names(z) %in% labels]] <- z[names(z) %in% labels] }; data.frame(field = rep(field, length(labels)), label = labels, color = unname(pal), stringsAsFactors = FALSE) }
colors <- rbind(colmap(all_long$category, "category"), colmap(all_long$condition, "condition"), colmap(all_long$batch, "batch")); write_tsv(colors, file.path(tech, "figure_colors.tsv"))
cat_colors <- setNames(colors$color[colors$field == "category"], colors$label[colors$field == "category"])
condition_colors <- setNames(colors$color[colors$field == "condition"], colors$label[colors$field == "condition"])

for (v in var_info) {
  x <- all_long[all_long$variable == v$column, , drop = FALSE]; x$category <- factor(x$category, levels = sort(unique(x$category)))
  if (length(unique(x$category)) < 2L) {
    for (family in c("composition_overview", "composition_dotplot", "composition_heatmap", "composition_counts", "embedding_diagnostics")) record(family, v$column, "skipped", "composition variable has fewer than two observed categories")
    next
  }
  title <- paste(v$label, "composition")
  p1 <- ggplot(x, aes(sample, proportion, fill = category)) + geom_col(width = .82, colour = "white", linewidth = .15) + scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) + scale_fill_manual(values = cat_colors, drop = FALSE) + labs(title = title, subtitle = paste0("One bar per sample; denominator: ", denom_desc), x = NULL, y = "Proportion") + theme_paper + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  n_condition <- length(unique(x$condition)); n_batch <- length(unique(x$batch))
  if (n_condition > 1 && n_batch > 1) p1 <- p1 + facet_grid(condition ~ batch, scales = "free_x", space = "free_x")
  else if (n_condition > 1) p1 <- p1 + facet_grid(. ~ condition, scales = "free_x", space = "free_x")
  else if (n_batch > 1) p1 <- p1 + facet_grid(. ~ batch, scales = "free_x", space = "free_x")
  save_plot(p1, file.path(out, paste0("composition_overview_", safe(v$column))), max(8, length(unique(x$sample)) * .35), 6); record("composition_overview", paste0(v$column))
  p2 <- ggplot(x, aes(condition, proportion)) + geom_jitter(aes(colour = condition), width = .12, height = 0, size = 2.2, alpha = .85) + stat_summary(fun = mean, geom = "point", shape = 95, size = 7, colour = "#303030") + facet_wrap(~ category, ncol = min(4, max(1, length(unique(x$category)))), scales = "free_y") + scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = condition_colors, drop = FALSE) + labs(title = paste(title, "by group"), subtitle = "Each point is a sample; black marker is the descriptive group mean", x = NULL, y = "Proportion", colour = "Condition") + theme_paper
  save_plot(p2, file.path(out, paste0("composition_dotplot_", safe(v$column))), 9, 6); record("composition_dotplot", v$column)
  hm <- aggregate(proportion ~ sample + category, x, mean); p3 <- ggplot(hm, aes(sample, category, fill = proportion)) + geom_tile(colour = "white", linewidth = .2) + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = .5, labels = scales::percent) + labs(title = paste(title, "heatmap"), subtitle = "Sample-level proportions", x = NULL, y = NULL, fill = "Proportion") + theme_paper + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p3, file.path(out, paste0("composition_heatmap_", safe(v$column))), max(8, length(unique(x$sample)) * .35), max(5, length(unique(x$category)) * .25)); record("composition_heatmap", v$column)
  p4 <- ggplot(x, aes(condition, n_cells)) + geom_jitter(aes(colour = condition), width = .12, height = 0, size = 2.2, alpha = .85) + facet_wrap(~ category, ncol = min(4, max(1, length(unique(x$category)))), scales = "free_y") + scale_y_continuous(labels = scales::comma) + scale_colour_manual(values = condition_colors, drop = FALSE) + labs(title = paste(title, "counts"), subtitle = "Counts are shown alongside relative composition", x = NULL, y = "Cells", colour = "Condition") + theme_paper
  save_plot(p4, file.path(out, paste0("composition_counts_", safe(v$column))), 9, 6); record("composition_counts", v$column)
  if (!is.null(obj) && isTRUE(getv(cfg, "plots.umap", TRUE)) && reduction %in% names(obj@reductions)) {
    fields <- unique(c(sample_col, condition_col, batch_col, v$column)); fields <- fields[!is.null(fields) & nzchar(fields) & fields %in% names(md)]
    fields <- fields[vapply(fields, function(f) length(unique(as.character(md[[f]]))) > 1L, logical(1))]
    if (!length(fields)) { record("embedding_diagnostics", v$column, "skipped", "no requested metadata field has at least two observed levels"); next }
    plots <- lapply(fields, function(f) DimPlot(obj, reduction = reduction, group.by = f, raster = FALSE, pt.size = .08) + theme_paper + labs(title = f))
    p5 <- wrap_plots(plots, ncol = 2); save_plot(p5, file.path(out, paste0("embedding_diagnostics_", safe(v$column))), 12, 5 * ceiling(length(plots) / 2)); record("embedding_diagnostics", v$column)
  } else record("embedding_diagnostics", v$column, "skipped", if (is.null(obj)) "requires Seurat object" else paste("reduction unavailable:", reduction))
}
write_tsv(status, file.path(out, "plot_status.tsv"));
artifact_paths <- list.files(out, recursive = FALSE, full.names = FALSE)
manifest <- list(schema_version = 1, skill = "14-scrna-visualize-cell-composition", project_id = getv(cfg, "project.id", "unknown"), input = input_object %||% counts_path, output_dir = out, denominator = list(mode = denom_mode, description = denom_desc), composition_variables = lapply(var_info, function(v) list(column = v$column, label = v$label, kind = v$kind)), artifacts = artifact_paths, finished_at = as.character(Sys.time()), status = "completed")
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(tech, "run_manifest.json")); writeLines(capture.output(sessionInfo()), file.path(tech, "session_info.txt"))
