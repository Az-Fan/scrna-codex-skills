# One-page PDF overview; preserves all input tables and the source Seurat object.
plot_marker_overview <- function(obj, markers, assay, cluster_col, config, out) {
  suppressPackageStartupMessages({library(ggplot2);library(Seurat)})
  groups <- unique(as.character(obj[[]][[cluster_col]]))
  if(anyNA(groups) || any(!nzchar(groups))) stop("Missing cluster identities")
  groups <- if(all(grepl("^[0-9]+$",groups))) groups[order(as.numeric(groups))] else sort(groups)
  obj[[cluster_col]] <- factor(as.character(obj[[]][[cluster_col]]),levels=groups)
  required <- c("cluster","gene","p_val_adj")
  if(length(setdiff(required,names(markers)))) stop("Existing marker table requires cluster, gene and p_val_adj")
  if(length(setdiff(as.character(markers$cluster),groups))) stop("Marker clusters do not match object clusters")
  fc <- intersect(c("avg_log2FC","avg_logFC","avg_diff"),names(markers))[1]
  if(is.na(fc)) stop("Marker table lacks an effect size")
  top_n <- cfg_get(config,"reporting.dotplot_top_n",4)
  if(!is.numeric(top_n) || length(top_n)!=1 || !is.finite(top_n) || top_n<1 || top_n!=floor(top_n)) stop("dotplot_top_n must be a positive integer")
  padj <- cfg_get(config,"reporting.adjusted_p_threshold",.05)
  if(!is.numeric(padj) || length(padj)!=1 || !is.finite(padj) || padj<=0 || padj>1) stop("adjusted_p_threshold must be in (0, 1]")
  if(!is.numeric(markers$p_val_adj) || !is.numeric(markers[[fc]])) stop("Marker statistics must be numeric")
  if(anyNA(markers$gene) || any(!nzchar(as.character(markers$gene)))) stop("Missing marker genes")
  passing <- markers[is.finite(markers$p_val_adj)&markers$p_val_adj<padj&is.finite(markers[[fc]]),,drop=FALSE]
  passing <- passing[order(match(as.character(passing$cluster),groups),passing$p_val_adj,-passing[[fc]],passing$gene),,drop=FALSE]
  passing <- passing[!duplicated(paste(passing$cluster,passing$gene,sep="\r")),,drop=FALSE]
  passing$overview_rank <- if(nrow(passing)) ave(seq_len(nrow(passing)),as.character(passing$cluster),FUN=seq_along) else integer()
  chosen <- passing[passing$overview_rank<=top_n,,drop=FALSE]
  if(any(!chosen$gene %in% rownames(obj[[assay]]))) stop("Selected markers missing from assay")
  chosen$row_id <- paste(chosen$cluster,chosen$gene,sep="::")
  utils::write.table(chosen,file.path(out,"marker_overview_selection.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  # Compute one background with all clusters and unique genes; repeated source blocks
  # reuse those identical values instead of rescaling expression within each block.
  dp <- if(nrow(chosen)) Seurat::DotPlot(obj,features=unique(chosen$gene),assay=assay,group.by=cluster_col)$data else data.frame(id=character(),features.plot=character(),avg.exp.scaled=numeric(),pct.exp=numeric())
  keys <- paste(as.character(dp$id),as.character(dp$features.plot),sep="\r")
  cell_counts <- table(obj[[]][[cluster_col]])
  entries <- lapply(seq_len(nrow(chosen)),function(i) {
    z<-dp[match(paste(groups,chosen$gene[i],sep="\r"),keys),,drop=FALSE]
    z$row_id<-chosen$row_id[i];z$source_cluster<-as.character(chosen$cluster[i]);z$marker_gene<-chosen$gene[i];z
  })
  z<-if(length(entries)) do.call(rbind,entries) else data.frame(id=character(),avg.exp.scaled=numeric(),pct.exp=numeric(),row_id=character(),source_cluster=character(),marker_gene=character())
  # Preserve a visible block for every cluster, including ones with no passing markers.
  missing_groups <- setdiff(groups,as.character(chosen$cluster))
  extra_ids <- character()
  if(length(missing_groups)) for(g in missing_groups) {
    zz<-data.frame(id=groups,avg.exp.scaled=NA_real_,pct.exp=NA_real_,row_id="",source_cluster=g,marker_gene="",stringsAsFactors=FALSE)
    for(nm in setdiff(names(z),names(zz))) zz[[nm]]<-NA
    zz<-zz[,names(z),drop=FALSE]
    zz$row_id<-paste0(g,"::no-marker");zz$source_cluster<-g;zz$marker_gene<-"No passing marker"
    extra_ids<-c(extra_ids,zz$row_id[1]);z<-rbind(z,zz)
  }
  z$id<-factor(as.character(z$id),levels=groups)
  z$source_cluster<-factor(z$source_cluster,levels=groups)
  row_order <- c(chosen$row_id,extra_ids)
  z$row_id<-factor(z$row_id,levels=rev(row_order))
  labels<-setNames(c(chosen$gene,rep("No passing marker",length(extra_ids))),row_order)
  strips<-setNames(paste0("C",groups,"\nn=",as.integer(cell_counts[groups])),groups)
  p<-ggplot(z,aes(id,row_id,color=avg.exp.scaled,size=pct.exp))+geom_point(na.rm=TRUE)+
    facet_grid(source_cluster~.,scales="free_y",space="free_y",labeller=labeller(source_cluster=strips))+
    scale_y_discrete(labels=function(v)unname(labels[v]))+
    scale_color_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0,name="Scaled mean")+
    scale_size(limits=c(0,100),range=c(0,4.2),breaks=c(25,50,75,100),name="Detected (%)")+
    theme_classic(base_size=10)+theme(axis.text.y=element_text(size=8.5),axis.text.x=element_text(size=9),
      strip.background=element_rect(fill="#F0F2F4",color=NA),strip.text.y=element_text(angle=0,size=8),
      panel.spacing.y=grid::unit(.10,"lines"),legend.position="bottom",legend.box="horizontal",
      plot.title=element_text(face="bold",size=14),plot.subtitle=element_text(size=9),plot.caption=element_text(size=8),
      axis.ticks.y=element_blank(),panel.grid.major.x=element_line(color="#F0F0F0",linewidth=.2))+
    guides(color=guide_colorbar(barwidth=grid::unit(3,"cm"),barheight=grid::unit(.3,"cm")),
      size=guide_legend(nrow=1))+
    labs(title="Cluster marker overview",subtitle=paste0("Top ",top_n," adjusted-significant markers per cluster | ",length(groups)," clusters | ",ncol(obj)," cells"),
      x="Cluster",y=NULL,caption="Right: marker-source cluster and cell count. Repeated genes retain identical expression values.\nDot size: detected fraction; color: gene-wise scaled mean RNA expression across all clusters.")
  width<-as.numeric(cfg_get(config,"plots.overview_width",max(7.5,3.2+.4*length(groups))))
  height<-as.numeric(cfg_get(config,"plots.overview_height",max(6,2.7+.185*length(row_order)+.08*length(groups))))
  if(length(width)!=1 || length(height)!=1 || !is.finite(width) || !is.finite(height) || width<=0 || height<=0) stop("Overview dimensions must be positive finite inches")
  format<-figure_format(config)
  file<-file.path(out,paste0("top_marker_dotplot.",format))
  ggsave(file,p,width=width,height=height,units="in",dpi=300,bg="white",limitsize=FALSE)
  paper_record(out,"top_marker_dotplot",file)
  dir.create(file.path(out,"_provenance"),showWarnings=FALSE)
  jsonlite::write_json(list(file=file,pages=1,format=format,layout="genes_by_clusters_source_blocks",top_n=top_n,
    marker_rows=nrow(chosen),unique_genes=length(unique(chosen$gene)),clusters=groups,
    zero_marker_clusters=missing_groups,assay=assay,scale_background="all clusters, unique selected genes",
    width_inches=width,height_inches=height),file.path(out,"_provenance","marker_overview_methods.json"),auto_unbox=TRUE,pretty=TRUE)
  file
}
