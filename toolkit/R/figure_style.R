# Fixed paper-style plotting helpers. No analysis or object mutation.
paper_theme <- function() {
  ggplot2::theme_classic(base_size = 11, base_family = "sans") +
    ggplot2::theme(text = ggplot2::element_text(colour = "#303030"),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(colour = "#303030", size = 9),
      axis.line = ggplot2::element_line(linewidth = .35),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      strip.background = ggplot2::element_rect(fill = "#F4F4F4", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 12, 8, 8))
}

paper_formats <- function(config) {
  formats <- tolower(unlist(cfg_get(config, "output.figure_format", "png"), use.names = FALSE))
  if (identical(formats, "both")) formats <- c("png", "pdf")
  if (!length(formats) || any(!formats %in% c("png", "pdf"))) stop("output.figure_format must be png, pdf or both")
  unique(formats)
}

paper_record <- function(out, family, file = "", status = "generated", reason = "") {
  path <- technical_path(out, "figure_status.tsv")
  row <- data.frame(family = family, file = file, status = status, reason = reason,
                    style_version = "paper_v1", stringsAsFactors = FALSE)
  previous <- if (file.exists(path)) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE) else row[FALSE, ]
  previous <- previous[!(previous$family == family & previous$file == file), , drop = FALSE]
  utils::write.table(rbind(previous, row), path, sep = "\t", quote = TRUE, row.names = FALSE)
}

paper_colors <- function(values, field, config, out) {
  labels <- sort(unique(as.character(values[!is.na(values)])), method = "radix")
  # Per-label hashing is independent of row order and which other labels are present.
  fallback <- vapply(labels, function(label) {
    hash <- 0
    for (code in utf8ToInt(enc2utf8(label))) hash <- (hash * 31 + code) %% 104729
    grDevices::hcl(h = (hash * 137.508) %% 360, c = 65, l = 58, fixup = TRUE)
  }, character(1))
  names(fallback) <- labels
  configured <- cfg_get(config, "figure_style.colors", list())[[field]]
  if (!is.null(configured)) {
    configured <- unlist(configured)
    if (is.null(names(configured)) || any(!nzchar(names(configured))) || anyDuplicated(names(configured))) stop("figure_style.colors maps must have unique label names")
    tryCatch(grDevices::col2rgb(configured), error = function(e) stop("Invalid configured figure colour: ", conditionMessage(e)))
    common <- intersect(labels, names(configured))
    fallback[common] <- configured[common]
  }
  path <- technical_path(out, "figure_colors.tsv")
  rows <- data.frame(field = field, label = labels, color = unname(fallback), stringsAsFactors = FALSE)
  if (file.exists(path)) {
    previous <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    previous <- previous[!(previous$field == field & previous$label %in% labels), , drop = FALSE]
    rows <- rbind(previous, rows)
  }
  utils::write.table(rows, path, sep = "\t", quote = TRUE, row.names = FALSE)
  fallback
}

paper_save <- function(plot, stem, width, height, config, out, family = basename(stem)) {
  files <- character()
  for (format in paper_formats(config)) {
    file <- paste0(stem, ".", format)
    withr::with_seed(1, ggplot2::ggsave(file, plot, width = width, height = height, units = "in", dpi = 300, bg = "white"))
    paper_record(out, family, file)
    files <- c(files, file)
  }
  files
}

paper_dimplot <- function(obj, reduction, groups, stem, config, out, point_size = .1, label = TRUE, legends = TRUE) {
  groups <- unique(groups)
  audit_fields <- unlist(cfg_get(config, "metadata", list())[c("sample", "condition", "batch")], use.names = FALSE)
  plots <- lapply(groups, function(field) {
    colors <- paper_colors(obj[[]][[field]], field, config, out)
    p <- Seurat::DimPlot(obj, reduction = reduction, group.by = field, cols = colors,
                         label = label && !field %in% audit_fields, repel = TRUE, pt.size = point_size, raster = FALSE,
                         label.size = 3, seed = 1) + paper_theme() + ggplot2::labs(title = field) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = max(1, ceiling(length(colors) / 20)), override.aes = list(size = 2)))
    if (!legends) p <- p + Seurat::NoLegend()
    p
  })
  files <- character()
  pages <- split(seq_along(plots), ceiling(seq_along(plots) / 4))
  for (i in seq_along(pages)) {
    index <- pages[[i]]
    page_stem <- if (length(pages) == 1L) stem else paste0(stem, sprintf("_page%02d", i))
    combined <- patchwork::wrap_plots(plots[index], ncol = min(2L, length(index)))
    files <- c(files, paper_save(combined, page_stem, 7 * min(2L, length(index)), 6 * ceiling(length(index) / 2), config, out, basename(stem)))
  }
  files
}

paper_dotplot <- function(obj, genes, assay, group, stem, config, out) {
  # Compute all scaled expression once; pagination must not change the colour scale.
  p <- Seurat::DotPlot(obj, features = genes, assay = assay, group.by = group) + paper_theme() +
    Seurat::RotatedAxis() + ggplot2::labs(x = "Gene", y = "Cluster", title = "Marker expression")
  data <- p$data
  color_range <- range(data$avg.exp.scaled, finite = TRUE)
  if (!all(is.finite(color_range))) color_range <- c(-1, 1)
  p <- p + ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = color_range) +
    ggplot2::scale_size(limits = c(0, 100), range = c(0, 6))
  gene_pages <- split(genes, ceiling(seq_along(genes) / 30))
  identities <- levels(data$id)
  if (!length(identities)) identities <- unique(as.character(data$id))
  group_pages <- split(identities, ceiling(seq_along(identities) / 25))
  files <- character(); page <- 0L; total <- length(gene_pages) * length(group_pages)
  for (selected_groups in group_pages) for (selected_genes in gene_pages) {
    page <- page + 1L
    q <- p
    q$data <- data[as.character(data$features.plot) %in% selected_genes & as.character(data$id) %in% selected_groups, , drop = FALSE]
    q <- q + ggplot2::scale_x_discrete(limits = selected_genes) + ggplot2::scale_y_discrete(limits = selected_groups) +
      ggplot2::labs(subtitle = paste0("Dot size: detected fraction; colour: scaled mean expression | page ", page, "/", total))
    page_stem <- if (total == 1L) stem else paste0(stem, sprintf("_page%02d", page))
    files <- c(files, paper_save(q, page_stem, max(8, length(selected_genes) * .28 + 3), max(5, length(selected_groups) * .24 + 2.5), config, out, basename(stem)))
  }
  files
}

paper_features <- function(obj, reduction, config, out) {
  genes <- unique(as.character(unlist(cfg_get(config, "figure_style.feature_genes", list()), use.names = FALSE)))
  if (!length(genes)) return(character())
  assay <- cfg_get(config, "figure_style.feature_assay", Seurat::DefaultAssay(obj))
  if (!assay %in% names(obj@assays)) stop("Feature plot assay not found: ", assay)
  Seurat::DefaultAssay(obj) <- assay
  layers <- SeuratObject::Layers(obj[[assay]])
  if (!"data" %in% layers) stop("Feature plots require a normalized data layer in figure_style.feature_assay; input was not normalized or modified by the plotting helper")
  missing <- setdiff(genes, rownames(obj))
  for (gene in missing) paper_record(out, "target_gene_featureplots", status = "skipped", reason = paste("Gene not found:", gene), file = gene)
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(character())
  files <- character(); pages <- split(genes, ceiling(seq_along(genes) / 6))
  for (i in seq_along(pages)) {
    plots <- Seurat::FeaturePlot(obj, features = pages[[i]], reduction = reduction,
      cols = c("#E8E8E8", "#B2182B"), order = TRUE, raster = FALSE, combine = FALSE, keep.scale = "feature")
    plots <- lapply(plots, function(p) p + paper_theme() + ggplot2::labs(subtitle = paste(assay, "normalized expression; per-gene scale")))
    stem <- file.path(out, if (length(pages) == 1) "target_gene_featureplots" else sprintf("target_gene_featureplots_page%02d", i))
    files <- c(files, paper_save(patchwork::wrap_plots(plots, ncol = min(3, length(plots))), stem,
      5 * min(3, length(plots)), 4.8 * ceiling(length(plots) / 3), config, out, "target_gene_featureplots"))
  }
  files
}
