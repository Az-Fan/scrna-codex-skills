# Shared-coordinate, plotting-only review of an existing Seurat embedding.
review_cluster_umap <- function(obj, config, out, prefix="") {
  suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
  fmt <- figure_format(config)
  red <- cfg_get(config,"metadata.reduction",required=TRUE); sm <- cfg_get(config,"metadata.sample",required=TRUE); cl <- cfg_get(config,"metadata.cluster",required=TRUE)
  if (!red %in% names(obj@reductions)) stop("Reduction not found: ", red)
  if (!all(c(sm, cl) %in% colnames(obj[[]]))) stop("Sample or cluster column missing")
  embedding <- Seurat::Embeddings(obj, red)
  if (ncol(embedding)<2L || anyDuplicated(rownames(embedding)) || !setequal(rownames(embedding),colnames(obj))) stop("Embedding must cover every input cell exactly once with at least two dimensions")
  xy <- embedding[colnames(obj), 1:2, drop = FALSE]
  md <- obj[[]][rownames(xy), , drop = FALSE]
  d <- data.frame(x = xy[,1], y = xy[,2], sample = as.character(md[[sm]]), cluster = as.character(md[[cl]]))
  if (anyNA(d) || any(!is.finite(xy)) || any(!nzchar(trimws(d$sample))) || any(!nzchar(trimws(d$cluster)))) stop("Missing/invalid embedding or labels")
  clusters <- unique(d$cluster)
  clusters <- if (all(!is.na(suppressWarnings(as.numeric(clusters))))) clusters[order(as.numeric(clusters))] else sort(clusters)
  samples <- unlist(cfg_get(config,"review_umap.sample_order",list()), use.names = FALSE)
  if (!length(samples)) samples <- unique(d$sample)
  if (anyDuplicated(samples) || !setequal(samples, unique(d$sample))) stop("sample_order must contain every sample exactly once")
  d$cluster <- factor(d$cluster, levels = clusters)
  configured <- cfg_get(config,"review_umap.cluster_colors",list())
  if(length(configured) && !setequal(names(configured),clusters)) stop("cluster_colors must cover every cluster")
  colors <- umap_palette(clusters,configured)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data.frame(cluster = clusters, color = colors[clusters]), file.path(out, paste0(prefix,"cluster_colors.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  tab <- as.data.frame(table(sample = factor(d$sample, levels = samples), cluster = d$cluster), stringsAsFactors = FALSE)
  names(tab)[3] <- "cells"
  tab$sample_fraction <- tab$cells / as.numeric(table(factor(d$sample, levels = samples)))[match(tab$sample, samples)]
  tab$cluster_fraction <- tab$cells / as.numeric(table(d$cluster))[match(tab$cluster, clusters)]
  stopifnot(sum(tab$cells) == nrow(d), all(abs(tapply(tab$sample_fraction, tab$sample, sum)-1) < 1e-10), all(abs(tapply(tab$cluster_fraction, tab$cluster, sum)-1) < 1e-10))
  utils::write.table(tab, file.path(out, paste0(prefix,"sample_cluster_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  xlim <- range(d$x); ylim <- range(d$y)
  base <- function() ggplot() + coord_equal(xlim = xlim, ylim = ylim) + theme_classic(base_size = 11) + labs(x = "UMAP 1", y = "UMAP 2")
  save_plot <- function(p, stem, width, height) ggplot2::ggsave(file.path(out, paste0(prefix,stem, ".", fmt)), p, width = width, height = height, units = "in", dpi = 300, bg = "white")
  point_size <- umap_point_size(nrow(d),cfg_get(config,"review_umap.point_size",cfg_get(config,"figure_style.umap_point_size")))
  if(!is.numeric(point_size) || length(point_size)!=1L || !is.finite(point_size) || point_size<=0) stop("review_umap.point_size must be a positive finite number")
  centres <- aggregate(cbind(x,y) ~ cluster, d, median)
  dd <- withr::with_seed(1,d[sample.int(nrow(d)), ])
  p <- base() + geom_point(data = dd, aes(x,y,color=cluster), size=point_size) + scale_color_manual(values=colors, drop=FALSE) + guides(color=guide_legend(override.aes=list(size=2))) + geom_label(data=centres, aes(x,y,label=cluster), size=3, fill="white", linewidth=.15) + labs(title="Cluster UMAP", color="Cluster")
  save_plot(p, "cluster_umap", 10, 7)
  plots <- lapply(samples, function(s) base() + geom_point(data=d, aes(x,y), color="grey85", size=.85*point_size, alpha=.5) + geom_point(data=d[d$sample==s,], aes(x,y), color="#0072B2", size=1.2*point_size, alpha=.88) + labs(title=paste0(s," (n=",sum(d$sample==s),")")))
  pages <- split(seq_along(plots), ceiling(seq_along(plots)/4))
  if(fmt=="pdf") {
    grDevices::pdf(file.path(out,paste0(prefix,"sample_highlights.pdf")),width=12,height=9,onefile=TRUE)
    on.exit(grDevices::dev.off(),add=TRUE)
  }
  for (i in seq_along(pages)) {
    p <- wrap_plots(plots[pages[[i]]], ncol=2) + plot_annotation(title="Sample representation across the shared UMAP", subtitle="Grey: all cells; blue: titled sample | identical coordinates and axes")
    if(fmt=="pdf") print(p) else save_plot(p, sprintf("sample_highlight_page%02d",i), 12, 9)
  }

  files <- c(paste0(prefix,"cluster_umap.",fmt),
    if(fmt=="pdf") paste0(prefix,"sample_highlights.pdf") else sprintf("%ssample_highlight_page%02d.png",prefix,seq_along(pages)),
    paste0(prefix,c("cluster_colors.tsv","sample_cluster_counts.tsv")))
  jsonlite::write_json(list(reduction=red,sample_column=sm,cluster_column=cl,cells=nrow(d),
    samples=samples,clusters=clusters,figure_format=fmt,palette=as.list(colors),
    sample_pages=length(pages),panels_per_page=4,point_size=point_size,shared_xlim=xlim,shared_ylim=ylim),
    technical_path(out,paste0(prefix,"umap_review_methods.json")),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  file.path(out,files)
}
