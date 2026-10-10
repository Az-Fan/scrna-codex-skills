suppressPackageStartupMessages({library(Seurat);library(jsonlite)})
source('toolkit/R/runtime.R');source('toolkit/R/figure_style.R');source('toolkit/R/plot_marker_overview.R')
set.seed(42)
counts<-matrix(rpois(16*130,4),16,130,dimnames=list(paste0('Gene',1:16),paste0('Cell',1:130)))
obj<-NormalizeData(CreateSeuratObject(counts),verbose=FALSE)
obj$cluster<-rep(as.character(0:12),each=10)
markers<-data.frame(cluster=rep(as.character(0:11),each=4),gene=rep(paste0('Gene',1:4),12),p_val_adj=.001,avg_log2FC=rep(4:1,12))
original<-serialize(obj,NULL)
tmp<-tempfile('marker-overview-');dir.create(tmp)
for(kind in c('pdf','png','none','labels')) {
  out<-file.path(tmp,kind);dir.create(out)
  cfg<-list(output=list(figure_format=if(kind=='png') 'png' else 'pdf'))
  input<-markers;target<-obj
  if(kind=='none') input$p_val_adj<-1
  if(kind=='labels') {
    target$cluster<-paste0('Label',target$cluster)
    input$cluster<-paste0('Label',input$cluster)
  }
  file<-plot_marker_overview(target,input,'RNA','cluster',cfg,out)
  stopifnot(file.exists(file),identical(serialize(obj,NULL),original))
  info<-read_json(file.path(out,'_provenance/marker_overview_methods.json'),simplifyVector=TRUE)
  selected<-read.delim(file.path(out,'marker_overview_selection.tsv'))
  stopifnot(info$pages==1,length(info$clusters)==13,info$marker_rows==if(kind=='none')0 else 48)
  stopifnot(nrow(selected)==if(kind=='none')0 else 48)
  stopifnot(length(info$zero_marker_clusters)==if(kind=='none')13 else 1)
  stopifnot(length(list.files(out,pattern='[.](pdf|png)$'))==1)
  if(kind!='png') {
    bytes<-readBin(file,'raw',n=file.info(file)$size)
    text<-paste(rawToChar(bytes,multiple=TRUE),collapse='')
    pages<-gregexpr('/Type /Page\\b',text,perl=TRUE,useBytes=TRUE)[[1]]
    stopifnot(length(pages)==1,pages[[1]]>0)
  }
  if(kind=='pdf') stopifnot(identical(info$clusters,as.character(0:12)))
}
cat('PASS single-page marker overview: 13 clusters, shared genes, zero-marker blocks, numeric and text IDs, PDF/PNG and unchanged object\n')
