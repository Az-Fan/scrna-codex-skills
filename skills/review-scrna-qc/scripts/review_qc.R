#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(jsonlite)})
`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: review_qc.R CONFIG.json")
cfg <- fromJSON(args[[1]], simplifyVector = FALSE)
out <- normalizePath(cfg$output_dir, mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
obj_path <- normalizePath(cfg$input$object, mustWork = TRUE)
obj <- if (grepl("\\.qs$", obj_path, ignore.case = TRUE)) {
  if (!requireNamespace("qs", quietly = TRUE)) stop("QS input requires qs in the selected pixi environment")
  qs::qread(obj_path)
} else readRDS(obj_path)
if (!inherits(obj, "Seurat")) stop("input.object is not a Seurat object")
md <- obj[[]]
sample_col <- cfg$metadata$sample
if (!sample_col %in% names(md)) stop("metadata.sample column not found: ", sample_col)
md$.sample <- as.character(md[[sample_col]])
optional_group <- c(condition = cfg$metadata$condition %||% NULL, batch = cfg$metadata$batch %||% NULL)
pick_column <- function(configured, candidates) {
  if (!is.null(configured) && configured %in% names(md)) return(configured)
  hit <- candidates[candidates %in% names(md)]
  if (length(hit)) hit[[1]] else NA_character_
}
cluster_col <- pick_column(cfg$metadata$cluster %||% NULL, c("seurat_clusters", "cluster", "clusters"))
annotation_col <- pick_column(cfg$metadata$annotation %||% NULL,
                              c("annotation_manual", "annotation_marker_pred", "annotation", "cell_type", "celltype"))

aliases <- list(
  n_genes = c("n_genes", "nFeature_RNA"), n_UMIs = c("n_UMIs", "nCount_RNA"),
  mito_frac = c("mito_frac", "percent.mt"), nuclear_frac = "nuclear_frac",
  ambient_frac = c("ambient_frac_decontx", "ambient_frac"), doublet_score = "doublet_score",
  hbb_score = "hbb_score", s_score = c("s_score", "S.Score"),
  g2m_score = c("g2m_score", "G2M.Score"))
resolve <- function(x) { hit <- x[x %in% names(md)]; if (length(hit)) hit[[1]] else NA_character_ }
cols <- vapply(aliases, resolve, character(1))
availability <- data.frame(metric = names(cols), source_column = unname(cols),
                           available = !is.na(cols), stringsAsFactors = FALSE)
write.table(availability, file.path(out, "metric_availability.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
available <- names(cols)[!is.na(cols)]
if (!length(available)) stop("No recognized QC metric columns were found")
for (m in available) md[[m]] <- suppressWarnings(as.numeric(md[[cols[[m]]]]))

qfun <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(n=0, mean=NA, median=NA, mad=NA, q01=NA, q05=NA, q25=NA, q75=NA, q95=NA, q99=NA))
  c(n=length(x), mean=mean(x), median=median(x), mad=mad(x),
    setNames(as.numeric(quantile(x, c(.01,.05,.25,.75,.95,.99), na.rm=TRUE)), c("q01","q05","q25","q75","q95","q99")))
}
long_stats <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$.sample), function(ii)
  do.call(rbind, lapply(available, function(m) data.frame(sample=md$.sample[ii[1]], metric=m,
                                                        as.list(qfun(md[[m]][ii])), check.names=FALSE)))))
write.table(long_stats, file.path(out, "qc_quantiles_by_sample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
summary_wide <- reshape(long_stats[,c("sample","metric","n","median","q05","q95")], idvar="sample", timevar="metric", direction="wide")
write.table(summary_wide, file.path(out, "qc_summary_by_sample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

directions <- c(n_genes="both", n_UMIs="both", mito_frac="upper", nuclear_frac="upper",
                ambient_frac="upper", doublet_score="upper", hbb_score="upper",
                s_score="descriptive", g2m_score="descriptive")
thresholds <- do.call(rbind, lapply(available, function(m) {
  z <- qfun(md[[m]]); direction <- directions[[m]]
  lower <- if (direction == "both") max(z[["q01"]], z[["median"]] - 3*z[["mad"]], na.rm=TRUE) else NA_real_
  upper <- if (direction %in% c("both","upper")) min(z[["q99"]], z[["median"]] + 3*z[["mad"]], na.rm=TRUE) else NA_real_
  override <- cfg$thresholds[[m]] %||% list()
  if (!is.null(override$lower)) lower <- as.numeric(override$lower)
  if (!is.null(override$upper)) upper <- as.numeric(override$upper)
  data.frame(metric=m, source_column=cols[[m]], direction=direction, candidate_lower=lower,
             candidate_upper=upper, method=if(length(override)) "configured_override" else "q01/q99_and_median_3MAD",
             approval="", decision="", notes="", stringsAsFactors=FALSE)
}))
write.table(thresholds, file.path(out, "threshold_review.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")

keep <- rep(TRUE, nrow(md))
for (i in seq_len(nrow(thresholds))) {
  if (thresholds$direction[i] == "descriptive") next
  x <- md[[thresholds$metric[i]]]
  if (!is.na(thresholds$candidate_lower[i])) keep <- keep & (is.na(x) | x >= thresholds$candidate_lower[i])
  if (!is.na(thresholds$candidate_upper[i])) keep <- keep & (is.na(x) | x <= thresholds$candidate_upper[i])
}
retention <- function(group, label) {
  a <- aggregate(cbind(total=rep(1,length(keep)), retained=as.integer(keep)), list(group=group), sum)
  a$retention_fraction <- a$retained/a$total; names(a)[1] <- label; a
}
rs <- retention(md$.sample, "sample")
write.table(rs, file.path(out, "candidate_retention_by_sample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
groups <- list()
for (nm in names(optional_group)) if (!is.null(optional_group[[nm]]) && optional_group[[nm]] %in% names(md))
  groups[[nm]] <- retention(as.character(md[[optional_group[[nm]]]]), nm)
rg <- if(length(groups)) do.call(rbind, lapply(names(groups), function(nm) {
  x <- groups[[nm]]; names(x)[1] <- "group"; x$group_type <- nm; x[,c("group_type","group","total","retained","retention_fraction")]
})) else data.frame(group_type=character(),group=character(),total=integer(),retained=integer(),retention_fraction=numeric())
write.table(rg, file.path(out, "candidate_retention_by_group.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

plot_long <- do.call(rbind, lapply(available, function(m) data.frame(sample=md$.sample, metric=m, value=md[[m]])))
atlas <- list(); plot_log <- list()
save_plot <- function(name, plot, width=12, height=8, reason="") {
  ggsave(file.path(out, name), plot, width=width, height=height, dpi=300, bg="white", limitsize=FALSE)
  atlas[[length(atlas)+1]] <<- plot
  plot_log[[length(plot_log)+1]] <<- data.frame(file=name, status="generated", reason=reason)
}
skip_plot <- function(name, reason) plot_log[[length(plot_log)+1]] <<-
  data.frame(file=name, status="skipped", reason=reason)

p_sample <- ggplot(plot_long, aes(sample, value)) +
  geom_violin(scale="width", trim=TRUE, linewidth=.1, na.rm=TRUE, fill="#7DB7D8") +
  geom_boxplot(width=.12, outlier.shape=NA, na.rm=TRUE, fill="white") +
  facet_wrap(~metric, scales="free_y", ncol=2) + theme_classic(base_size=9) +
  theme(axis.text.x=element_text(angle=60,hjust=1)) + labs(x="Sample",y="QC value")
save_plot("qc_distribution_by_sample.png", p_sample, 12, max(8,3*ceiling(length(available)/2)))

p_dist <- ggplot(plot_long,aes(value)) + geom_histogram(bins=60,fill="#4C78A8",color="white",na.rm=TRUE) +
  facet_wrap(~metric,scales="free",ncol=3) + theme_classic(base_size=9) + labs(x=NULL,y="Cells")
save_plot("qc_distributions.png",p_dist,12,max(6,3*ceiling(length(available)/3)))

group_violin <- function(column, filename, xlab) {
  z <- do.call(rbind, lapply(available, function(m) data.frame(group=as.character(md[[column]]),metric=m,value=md[[m]])))
  p <- ggplot(z,aes(group,value)) + geom_violin(scale="width",trim=TRUE,linewidth=.1,fill="#9CCB86",na.rm=TRUE) +
    geom_boxplot(width=.12,outlier.shape=NA,fill="white",na.rm=TRUE) + facet_wrap(~metric,scales="free_y",ncol=2) +
    theme_classic(base_size=9) + theme(axis.text.x=element_text(angle=60,hjust=1)) + labs(x=xlab,y="QC value")
  save_plot(filename,p,12,max(8,3*ceiling(length(available)/2)),paste("group column:",column))
}
if (!is.na(cluster_col)) group_violin(cluster_col,"qc_metrics_per_cluster.png","Cluster") else
  skip_plot("qc_metrics_per_cluster.png","No cluster metadata column")
if (!is.na(annotation_col)) group_violin(annotation_col,"qc_distribution_by_celltype.png","Annotation") else
  skip_plot("qc_distribution_by_celltype.png","No annotation metadata column")

if (all(c("n_UMIs","n_genes") %in% available)) {
  p_scatter <- ggplot(md,aes(n_UMIs,n_genes)) + geom_point(size=.2,alpha=.35,color="#4C78A8") +
    facet_wrap(~.sample,scales="free") + scale_x_log10() + scale_y_log10() + theme_classic(base_size=9)
  save_plot("n_UMIs_vs_n_genes.png",p_scatter,12,8)
} else skip_plot("n_UMIs_vs_n_genes.png","Requires n_UMIs and n_genes")
if (all(c("nuclear_frac","n_UMIs") %in% available)) {
  p_nuclear <- ggplot(md,aes(nuclear_frac,n_UMIs)) + geom_point(size=.25,alpha=.4,color="#4C78A8") +
    scale_y_log10() + facet_wrap(~.sample,scales="free") + theme_classic(base_size=9)
  save_plot("nuclear_frac_vs_n_UMIs.png",p_nuclear,12,8)
} else skip_plot("nuclear_frac_vs_n_UMIs.png","Requires nuclear_frac and n_UMIs")

if ("mito_frac" %in% available && !is.na(annotation_col)) {
  mt <- thresholds[thresholds$metric=="mito_frac",,drop=FALSE]
  md$.annotation <- as.character(md[[annotation_col]])
  md$.candidate_keep <- if(nrow(mt) && !is.na(mt$candidate_upper)) is.na(md$mito_frac) | md$mito_frac <= mt$candidate_upper else TRUE
  p_mito <- ggplot(md,aes(.annotation,mito_frac,color=.candidate_keep)) + geom_jitter(size=.15,alpha=.35,width=.25) +
    {if(nrow(mt) && !is.na(mt$candidate_upper)) geom_hline(yintercept=mt$candidate_upper,linetype="dashed") else NULL} +
    coord_flip() + theme_classic(base_size=9) + labs(x="Annotation",y="mito_frac",color="Candidate keep")
  save_plot("filter_mito_by_annotation.png",p_mito,8,max(6,.22*length(unique(md$.annotation))))
} else skip_plot("filter_mito_by_annotation.png","Requires mito_frac and annotation")

umap_xy <- NULL
if ("umap" %in% names(obj@reductions)) umap_xy <- Embeddings(obj,"umap")[rownames(md),1:2,drop=FALSE]
if (is.null(umap_xy) && all(c("UMAP_1","UMAP_2") %in% names(md))) umap_xy <- as.matrix(md[,c("UMAP_1","UMAP_2")])
if (!is.null(umap_xy)) {
  md$.UMAP1 <- umap_xy[,1]; md$.UMAP2 <- umap_xy[,2]
  if (!is.na(annotation_col)) {
    p_ann <- ggplot(md,aes(.UMAP1,.UMAP2,color=.data[[annotation_col]])) + geom_point(size=.08,alpha=.7) +
      theme_void() + guides(color=guide_legend(override.aes=list(size=2))) + labs(color="Annotation")
    save_plot("umap_annotation_qc_context.png",p_ann,8,6)
  } else skip_plot("umap_annotation_qc_context.png","No annotation metadata column")
  for (m in available) {
    p_umap <- ggplot(md,aes(.UMAP1,.UMAP2,color=.data[[m]])) + geom_point(size=.08,alpha=.7) +
      scale_color_gradient(low="grey90",high="#B2182B",na.value="grey80") + theme_void() + labs(title=paste("UMAP:",m),color=m)
    save_plot(paste0("umap_qc_",m,".png"),p_umap,8,6)
  }
} else {
  skip_plot("umap_annotation_qc_context.png","No UMAP reduction or UMAP_1/UMAP_2 metadata")
  for(m in available) skip_plot(paste0("umap_qc_",m,".png"),"No UMAP coordinates")
}

p_ret <- ggplot(rs,aes(sample,retention_fraction)) + geom_col(fill="#59A14F") + geom_hline(yintercept=.8,lty=2) +
  coord_cartesian(ylim=c(0,1)) + theme_classic() + theme(axis.text.x=element_text(angle=60,hjust=1)) +
  labs(x="Sample",y="Candidate retention")
save_plot("candidate_retention_by_sample.png",p_ret,10,5)
write.table(do.call(rbind,plot_log),file.path(out,"plot_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
pdf(file.path(out,"qc_atlas.pdf"),width=12,height=8,onefile=TRUE)
for(p in atlas) print(p)
dev.off()
message("QC review completed: ", out)
