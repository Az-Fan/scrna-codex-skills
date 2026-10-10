#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(jsonlite)})
`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: review_qc.R CONFIG.json")
cfg <- fromJSON(args[[1]], simplifyVector = FALSE)
out <- normalizePath(cfg$output_dir, mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
detail_level <- tolower(cfg$output$detail_level %||% "compact")
if (!detail_level %in% c("compact", "full")) stop("output.detail_level must be compact or full")
full_output <- identical(detail_level, "full")
preview_png <- cfg$output$preview_png %||% FALSE
if (!is.logical(preview_png) || length(preview_png) != 1L || is.na(preview_png)) stop("output.preview_png must be true or false")
detail_out <- file.path(out, "details")
if (full_output) dir.create(detail_out, recursive = TRUE, showWarnings = FALSE)
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
if (anyNA(md$.sample) || any(!nzchar(trimws(md$.sample)))) stop("Sample IDs must be present for every cell")
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
numeric_values <- function(x) suppressWarnings(as.numeric(if (is.factor(x)) as.character(x) else x))
resolve <- function(x) { hit <- x[x %in% names(md)]; if (length(hit)) hit[[1]] else NA_character_ }
cols <- vapply(aliases, resolve, character(1))
finite_n <- vapply(cols, function(column) {
  if (is.na(column)) return(0L)
  sum(is.finite(numeric_values(md[[column]])))
}, integer(1))
availability <- data.frame(metric = names(cols), source_column = unname(cols), finite_cells = finite_n,
                           available = !is.na(cols) & finite_n > 0L,
                           reason = ifelse(is.na(cols), "column_not_found", ifelse(finite_n == 0L, "no_finite_values", "")),
                           stringsAsFactors = FALSE)
if (full_output) write.table(availability, file.path(detail_out, "metric_availability.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
available <- availability$metric[availability$available]
if (!length(available)) stop("No recognized QC metric columns were found")
for (m in available) {
  md[[m]] <- numeric_values(md[[cols[[m]]]])
  md[[m]][!is.finite(md[[m]])] <- NA_real_
}
for (column in c(cluster_col, annotation_col)) if (!is.na(column)) {
  values <- as.character(md[[column]])
  missing <- is.na(values) | !nzchar(trimws(values))
  if (any(missing)) {
    label <- "(missing)"
    while (label %in% values[!missing]) label <- paste0(label, "*")
    values[missing] <- label
  }
  md[[column]] <- values
}

qfun <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(n=0, mean=NA, median=NA, mad=NA, q01=NA, q05=NA, q25=NA, q75=NA, q95=NA, q99=NA))
  c(n=length(x), mean=mean(x), median=median(x), mad=mad(x),
    setNames(as.numeric(quantile(x, c(.01,.05,.25,.75,.95,.99), na.rm=TRUE)), c("q01","q05","q25","q75","q95","q99")))
}
long_stats <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$.sample), function(ii)
  do.call(rbind, lapply(available, function(m) data.frame(sample=md$.sample[ii[1]], metric=m,
                                                        as.list(qfun(md[[m]][ii])), check.names=FALSE)))))
if (full_output) write.table(long_stats, file.path(detail_out, "qc_quantiles_by_sample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
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
if (full_output) write.table(rs, file.path(detail_out, "candidate_retention_by_sample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
groups <- list()
for (nm in names(optional_group)) if (!is.null(optional_group[[nm]]) && optional_group[[nm]] %in% names(md))
  groups[[nm]] <- retention(as.character(md[[optional_group[[nm]]]]), nm)
rg <- if(length(groups)) do.call(rbind, lapply(names(groups), function(nm) {
  x <- groups[[nm]]; names(x)[1] <- "group"; x$group_type <- nm; x[,c("group_type","group","total","retained","retention_fraction")]
})) else data.frame(group_type=character(),group=character(),total=integer(),retained=integer(),retention_fraction=numeric())
if (full_output) write.table(rg, file.path(detail_out, "candidate_retention_by_group.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

suppressPackageStartupMessages(library(patchwork))
base_theme <- theme_classic(base_size=10) + theme(plot.title=element_text(face="bold",size=11),
  plot.subtitle=element_text(size=9,color="#555555"),legend.title=element_text(size=9),
  legend.text=element_text(size=8),plot.margin=margin(5,7,5,5))
theme_set(base_theme)
blank <- function(title, note) ggplot()+annotate("text",x=0,y=0,label=note,size=3)+theme_void()+labs(title=title)
ordered <- function(x) { z<-unique(as.character(x)); if(all(grepl("^[0-9]+$",z))) z[order(as.numeric(z))] else sort(z) }
samples <- unique(md$.sample)
labels <- vapply(samples,function(s) cfg$display$sample_labels[[s]] %||% s,character(1))
if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Sample display labels cannot be empty")
if(anyDuplicated(labels)) stop("Sample display labels must be unique")
md$.label <- factor(labels[match(md$.sample,samples)],levels=labels)
pal <- setNames(vapply(seq_along(labels),function(i) cfg$display$sample_colors[[labels[i]]] %||% grDevices::hcl.colors(length(labels),"Dark 3")[i],character(1)),labels)
mt_unit <- cfg$display$mito_unit %||% if(identical(cols[["mito_frac"]],"percent.mt")) "percent" else "fraction"
if(!mt_unit %in% c("percent","fraction")) stop("display.mito_unit must be percent or fraction")
if (length(cfg$approved_rules)) {
  if (!is.character(cfg$approved_rules_source) || length(cfg$approved_rules_source) != 1L || !nzchar(trimws(cfg$approved_rules_source))) stop("approved_rules_source is required for approved-rule overlays")
  for (m in names(cfg$approved_rules)) {
    if (!m %in% names(aliases)) stop("Unknown approved-rule metric: ", m)
    bounds <- cfg$approved_rules[[m]]
    if (!length(bounds) || any(!names(bounds) %in% c("lower", "upper")) || any(!vapply(bounds, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1)))) stop("Approved bounds must be finite numeric lower/upper values")
    if (!is.null(bounds$lower) && !is.null(bounds$upper) && bounds$lower > bounds$upper) stop("Approved lower bound exceeds upper bound")
  }
}
core <- intersect(c("n_genes","n_UMIs","mito_frac","doublet_score"),available)
metric_labels <- c(n_genes="Detected genes",n_UMIs="UMIs",mito_frac="Mitochondrial (%)",doublet_score="Doublet score")
display_value <- function(m,x) if(m=="mito_frac" && mt_unit=="fraction") x*100 else x
positive_log <- function(x) { x[!is.finite(x) | x <= 0] <- NA_real_; log10(x) }
rule <- function(m) {
  a<-cfg$approved_rules[[m]] %||% list()
  v<-unlist(a[c("lower","upper")],use.names=FALSE)
  display_value(m,v)
}
violin <- function(m,group=".label",data=md) {
  z<-data.frame(group=data[[group]],value=display_value(m,data[[m]]))
  p<-ggplot(z,aes(group,value,fill=group))+geom_violin(scale="width",linewidth=.25,trim=TRUE,na.rm=TRUE)+
    geom_boxplot(width=.12,outlier.shape=NA,fill="white",linewidth=.25,na.rm=TRUE)+
    labs(title=metric_labels[[m]],x=NULL,y=NULL)+theme(legend.position="none")
  if(group==".label") p<-p+scale_fill_manual(values=pal)
  else p<-p+scale_fill_manual(values=setNames(rep("#B5C8D1",nlevels(z$group)),levels(z$group)))
  if(m %in% c("n_genes","n_UMIs")) p<-p+scale_y_log10(labels=scales::label_number(big.mark=","))
  if(length(rule(m))) p<-p+geom_hline(yintercept=rule(m),linetype=2,color="#B34444",linewidth=.35)
  p
}
summary_display <- data.frame(sample=samples,label=labels,cells=vapply(samples,function(s)sum(md$.sample==s),integer(1)))
for(m in core) summary_display[[m]]<-vapply(samples,function(s)median(display_value(m,md[[m]][md$.sample==s]),na.rm=TRUE),numeric(1))
write.table(summary_display,file.path(out,"sample_display_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
sample_page <- function(selected) {
md <- md[md$.sample %in% selected, , drop=FALSE]
samples <- selected
labels <- labels[match(selected, summary_display$sample)]
md$.label <- factor(md$.label, levels=labels)
summary_display <- summary_display[match(selected, summary_display$sample), , drop=FALSE]
# A compact vector table shares the first row with the sample size bar chart.
headers<-c("Sample","Cells",unname(c(n_genes="Genes",n_UMIs="UMIs",mito_frac="Mito (%)",doublet_score="Doublet")[core]))
vals<-cbind(labels,format(summary_display$cells,big.mark=",",trim=TRUE),do.call(cbind,lapply(core,function(m)format(round(summary_display[[m]],if(m %in% c("n_genes","n_UMIs"))0 else 2),big.mark=",",trim=TRUE))))
text_data<-rbind(data.frame(x=seq_along(headers),y=length(samples)+1,text=headers,weight="bold"),
  data.frame(x=rep(seq_along(headers),each=length(samples)),y=rep(rev(seq_along(samples)),length(headers)),text=as.vector(vals),weight="plain"))
tab<-ggplot(text_data,aes(x,y,label=text,fontface=weight))+geom_text(size=2.9)+
  geom_hline(yintercept=length(samples)+.5,linewidth=.3,color="#999999")+theme_void()+
  coord_cartesian(xlim=c(.4,length(headers)+.6),ylim=c(.5,length(samples)+1.5))+labs(title="Sample summary · medians")
bar<-ggplot(summary_display,aes(factor(label,levels=labels),cells,fill=label))+geom_col(width=.65)+
  geom_text(aes(label=scales::comma(cells)),vjust=-.4,size=3)+scale_fill_manual(values=pal)+
  scale_y_continuous(expand=expansion(mult=c(0,.16)),labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  labs(title="Cells per sample",x=NULL,y="Cells")+theme(legend.position="none")
vp<-if(length(core)) wrap_plots(lapply(core,function(m)violin(m,data=md)),nrow=1) else blank("Core QC", "Core metrics unavailable; see supplementary diagnostics")
density <- if(all(c("n_UMIs","n_genes") %in% available)) {
  ggplot(md,aes(n_UMIs,n_genes))+geom_bin_2d(bins=65)+facet_wrap(~.label,nrow=1)+
    scale_x_log10(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+scale_y_log10(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
    scale_fill_viridis_c(trans="log10",name="Cells / bin")+labs(title="UMI-gene density · shared axes",x="UMIs (log scale)",y="Detected genes (log scale)")+
    geom_vline(xintercept=rule("n_UMIs"),linetype=2,color="#B34444",linewidth=.3)+
    geom_hline(yintercept=rule("n_genes"),linetype=2,color="#B34444",linewidth=.3)
} else blank("UMI-gene density","Required metrics unavailable")
rule_text<-if(length(cfg$approved_rules)) "Dashed red lines: previously approved rules; shown for context only. No filtering in this review." else "No approved rules configured. Candidate thresholds are available separately for review."
((bar|tab)/vp/density)+plot_layout(heights=c(.75,1.15,1.3))+
  plot_annotation(title=sprintf("QC overview  |  %s cells · %s genes · %s samples",scales::comma(ncol(obj)),scales::comma(nrow(obj)),length(samples)),
    subtitle="Counts use log scales; mitochondrial values use percent; boxes show median and interquartile range.",caption=rule_text,
    theme=theme(plot.title=element_text(face="bold",size=16),plot.caption=element_text(size=8)))
}
sample_pages <- lapply(split(samples, ceiling(seq_along(samples)/8)), sample_page)
cluster_pages<-list(blank("Cluster QC","Cluster metadata unavailable"))
if(!is.na(cluster_col)) {
  cl<-ordered(md[[cluster_col]]); md$.cluster<-factor(as.character(md[[cluster_col]]),levels=cl)
  cs<-data.frame(cluster=cl,cells=as.integer(table(md$.cluster)))
  for(m in available) cs[[paste0("median_",m)]]<-vapply(cl,function(k)median(md[[m]][md$.cluster==k],na.rm=TRUE),numeric(1))
  write.table(cs,file.path(out,"qc_summary_by_cluster.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  composition<-as.data.frame(table(cluster=md$.cluster,sample=md$.label));names(composition)[3]<-"cells"
  composition$fraction<-composition$cells/cs$cells[match(composition$cluster,cl)]
  composition$sample_id <- samples[match(as.character(composition$sample), labels)]
  write.table(composition,file.path(out,"cluster_sample_composition.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  cluster_page <- function(selected) {
  cl <- selected
  cs <- cs[match(cl, cs$cluster), , drop=FALSE]
  composition <- composition[composition$cluster %in% cl, , drop=FALSE]
  hm<-do.call(rbind,lapply(core,function(m) {
    raw<-display_value(m,cs[[paste0("median_",m)]]); z<-if(m %in% c("n_genes","n_UMIs"))positive_log(raw) else raw
    finite <- is.finite(z)
    sdv<-sd(z[finite]); z<-if(is.finite(sdv)&&sdv>0)(z-mean(z[finite]))/sdv else ifelse(finite,0,NA_real_)
    data.frame(cluster=factor(cl,levels=rev(cl)),metric=factor(metric_labels[[m]],levels=metric_labels[core]),z=z,
      text=if(m %in% c("n_genes","n_UMIs"))scales::comma(round(raw)) else sprintf("%.2f",raw))
  }))
  heat<-if(length(core)) ggplot(hm,aes(metric,cluster,fill=z))+geom_tile(color="white",linewidth=.4)+geom_text(aes(label=text),size=2.8)+
    scale_fill_gradient2(low="#4287AC",mid="white",high="#CE6B54",midpoint=0,na.value="grey85",name="Column z-score")+
    labs(title="QC medians",x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=25,hjust=1),axis.ticks=element_blank(),legend.position="bottom")
    else blank("QC medians", "Core metrics unavailable")
  cb<-ggplot(cs,aes(cells,factor(cluster,levels=rev(cl))))+geom_col(fill="#778C9B",width=.75)+
    scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+labs(title="Cluster size",x="Cells",y="Cluster")
  comp<-ggplot(composition,aes(fraction,factor(cluster,levels=rev(cl)),fill=sample))+geom_col(width=.8)+
    scale_fill_manual(values=pal)+scale_x_continuous(labels=scales::label_percent(),expand=c(0,0))+
    labs(title="Sample composition",x="Fraction of cluster",y=NULL,fill=NULL)+theme(legend.position="bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank())
  focus<-intersect(unlist(cfg$display$focus_clusters %||% list()),cl)
  focusplots<-if(length(focus) && length(core)) {zz<-md[md$.cluster %in% focus,,drop=FALSE];zz$.cluster<-factor(zz$.cluster,levels=focus);wrap_plots(lapply(core,function(m)violin(m,".cluster",zz)),nrow=1)} else blank("Focused review","No focus clusters or core metrics configured")
  ((cb+heat+comp+plot_layout(widths=c(.7,1.8,1.2)))/focusplots)+plot_layout(heights=c(2,1))+
    plot_annotation(title="Cluster QC  |  complexity, composition and focused review",
      subtitle=paste0("Numeric cluster order; heatmap prints raw medians. Colors standardize each metric (log10 counts before z-score). Focus: ",paste(focus,collapse=", "),"."),
      caption="QC and sample imbalance describe clusters; neither alone establishes a low-quality population or a batch effect.",
      theme=theme(plot.title=element_text(face="bold",size=16),plot.caption=element_text(size=8)))
  }
  cluster_pages <- lapply(split(cl, ceiling(seq_along(cl)/24)), cluster_page)
}
# Rasterize only dense points; PDF titles, legends, axes and cluster labels remain vector.
red<-cfg$display$umap_reduction %||% intersect(c("umap","standard_umap"),names(obj@reductions))[1]
page3<-blank("UMAP QC","UMAP coordinates unavailable")
xy <- NULL
if(length(red)&&!is.na(red)&&red %in% names(obj@reductions)) {
  embedding <- Embeddings(obj,red)
  if(ncol(embedding) >= 2L && all(rownames(md) %in% rownames(embedding))) xy <- embedding[rownames(md),1:2,drop=FALSE]
} else if(is.null(cfg$display$umap_reduction) && all(c("UMAP_1","UMAP_2") %in% names(md))) {
  xy <- cbind(numeric_values(md$UMAP_1), numeric_values(md$UMAP_2)); red <- "metadata:UMAP_1/UMAP_2"
}
if(!is.null(xy) && all(is.finite(xy))) {
  xr<-range(xy[,1]);yr<-range(xy[,2]);xr<-xr+c(-1,1)*max(diff(xr)*.025,.001);yr<-yr+c(-1,1)*max(diff(yr)*.025,.001)
  raster_umap<-function(values,title,discrete=FALSE,colors=NULL,cluster_labels=FALSE) {
    if(discrete) {lev<-levels(values);cc<-unname(colors[as.character(values)]);leg<-data.frame(x=mean(xr),y=mean(yr),v=factor(lev,levels=lev))}
    else {
      if(!any(is.finite(values))) return(blank(title,"No finite values on the displayed scale"))
      lim<-range(values[is.finite(values)]);ramp<-viridisLite::viridis(256); idx<-if(diff(lim)>0)pmax(1,pmin(256,1+floor(255*(values-lim[1])/diff(lim)))) else rep(128,length(values));cc<-ramp[idx];cc[!is.finite(values)]<-"#BBBBBB";leg<-data.frame(x=mean(xr),y=mean(yr),v=lim)
    }
    tf<-tempfile(fileext=".png");png(tf,width=1100,height=1100,bg="white")
    par(mar=c(0,0,0,0),xaxs="i",yaxs="i");plot.new();plot.window(xlim=xr,ylim=yr)
    set.seed(777); ii<-sample(seq_len(nrow(xy)))
    point_cex <- min(2, max(.28, 20/sqrt(nrow(xy))))
    points(xy[ii,1],xy[ii,2],pch=16,cex=point_cex,col=cc[ii]);dev.off()
    raster<-png::readPNG(tf);unlink(tf)
    p<-ggplot(leg,aes(x,y,color=v))+annotation_raster(raster,xr[1],xr[2],yr[1],yr[2])+geom_point(size=0)+
      coord_fixed(xlim=xr,ylim=yr,expand=FALSE)+labs(title=title,x="UMAP 1",y="UMAP 2",color=NULL)+
      theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=grid::unit(.28,"cm"),legend.key.width=grid::unit(.35,"cm"))
    if(discrete) p<-p+scale_color_manual(values=colors,drop=FALSE)+guides(color=guide_legend(override.aes=list(size=2),ncol=if(length(lev)>5)2 else 1))
    else p<-p+scale_color_viridis_c(limits=lim)+guides(color=guide_colorbar(barheight=grid::unit(2.2,"cm"),barwidth=grid::unit(.3,"cm")))
    if(cluster_labels) {cent<-aggregate(xy,list(cluster=values),median);p<-p+geom_label(data=cent,aes(x=.data[[names(cent)[2]]],y=.data[[names(cent)[3]]],label=cluster),inherit.aes=FALSE,size=2.5,linewidth=.1,fill="white")}
    p
  }
  uc<-if(!is.na(cluster_col))raster_umap(md$.cluster,"Cluster",TRUE,setNames(hcl.colors(length(cl),"Dark 3"),cl),TRUE) else blank("Cluster","No cluster metadata")
  us<-raster_umap(md$.label,"Sample",TRUE,pal)
  um<-lapply(core,function(m)raster_umap(if(m %in% c("n_genes","n_UMIs"))positive_log(md[[m]]) else display_value(m,md[[m]]),
    paste0(metric_labels[[m]],if(m %in% c("n_genes","n_UMIs"))" · log10" else "")))
  page3<-wrap_plots(c(list(uc,us),um),ncol=3)+plot_annotation(title="UMAP QC  |  shared coordinates across all panels",
    subtitle="Full cell set; identical axes and point order. Continuous colors span all finite values without percentile clipping.",
    caption="Doublet score is a diagnostic score, not a doublet probability or an approved exclusion rule.",theme=theme(plot.title=element_text(face="bold",size=16),plot.caption=element_text(size=8)))
}
pages<-c(sample_pages,cluster_pages,list(page3))
pdf(file.path(out,"qc_atlas.pdf"),width=12,height=9,onefile=TRUE)
for(p in pages) print(p)
dev.off()
page_names <- function(stem,n) if(n==1L) paste0(stem,".png") else paste0(stem,"_page",seq_len(n),".png")
preview_names <- if(preview_png) c(page_names("qc_01_samples",length(sample_pages)),page_names("qc_02_clusters",length(cluster_pages)),"qc_03_umap.png") else character()
if(preview_png) for(i in seq_along(pages)) ggsave(file.path(out,preview_names[i]),pages[[i]],width=12,height=9,dpi=300,bg="white")
# Preserve all established annotation, secondary metric and candidate-retention
# diagnostics in a separate atlas; full mode retains their individual PNGs.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"review_qc_diagnostics.R"),local=TRUE)
pdf(file.path(out,"qc_supplement.pdf"),width=12,height=8,onefile=TRUE)
for(p in atlas) print(p)
dev.off()
dir.create(file.path(out,"_provenance"),showWarnings=FALSE)
write_json(list(cells=ncol(obj),genes=nrow(obj),clusters=if(!is.na(cluster_col))length(cl) else NULL,
  source_object=obj_path,new_seurat_objects=0,cells_removed=0,mitochondrial_source_unit=mt_unit,mitochondrial_display_unit="percent",
  umap_reduction=red,sample_labels=as.list(setNames(labels,samples)),sample_palette=as.list(pal),
  approved_rules=cfg$approved_rules,approved_rules_source=cfg$approved_rules_source,
  threshold_candidates_applied=FALSE,point_raster_pixels=1100,point_order_seed=777,png_dpi=if(preview_png)300 else NULL,preview_png=preview_png,
  atlas_pages=length(pages),sample_pages=length(sample_pages),cluster_pages=length(cluster_pages),previews=preview_names,
  supplemental_diagnostics=vapply(Filter(function(x)x$status=="generated",plot_log),function(x)x$file,character(1))),file.path(out,"_provenance","figure_methods.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"_provenance","session_info.txt"))
message("Compact QC overview and supplemental diagnostics complete: ",out)
