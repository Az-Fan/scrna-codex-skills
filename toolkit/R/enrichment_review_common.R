# Shared preflight for presentation-only enrichment reviews.
source(file.path(dirname(renderer_path), 'figure_output.R'))
review_read <- function(path, cfg) {
  x <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE, colClasses='character', na.strings=c('', 'NA'))
  for (field in c('population','comparison_id','task_id','input_table')) {
    value <- if(field=="input_table")cfg[["source_table"]]else cfg[[field]]
    if (!is.null(value) && field %in% names(x)) x <- x[!is.na(x[[field]]) & x[[field]] == value,,drop=FALSE]
  }
  for (field in intersect(c('task_id','input_table'),names(x))) {
    if (length(unique(x[[field]][!is.na(x[[field]])])) > 1L)
      stop('Ambiguous input: explicitly select ',field,' in review configuration')
  }
  for (field in intersect(c('p.adjust','NES','Count','RichFactor','FoldEnrichment','setSize','log2FoldChange'),names(x))) {
    value <- suppressWarnings(as.numeric(x[[field]]))
    if (any(!is.na(x[[field]]) & is.na(value))) stop('Nonnumeric scientific column: ',field)
    x[[field]] <- value
  }
  # Derive display metrics from the original ORA ratios when older tables lack them.
  if('method'%in%names(x)) {
    for(field in c('NES','Count','geneID','core_enrichment','direction'))if(!field%in%names(x))x[[field]]<-rep(NA,nrow(x))
    ora<-which(x$method=='ORA')
    if(length(ora)&&all(c('BgRatio','GeneRatio')%in%names(x))) {
      ratio<-function(z,part) suppressWarnings(as.numeric(vapply(strsplit(z,'/',fixed=TRUE),function(v)if(length(v)==2)v[part]else NA_character_,character(1))))
      if(!'RichFactor'%in%names(x))x$RichFactor<-rep(NA_real_,nrow(x))
      ix<-ora[!is.finite(x$RichFactor[ora])];x$RichFactor[ix]<-x$Count[ix]/ratio(x$BgRatio[ix],1)
      if(!'FoldEnrichment'%in%names(x))x$FoldEnrichment<-rep(NA_real_,nrow(x))
      ix<-ora[!is.finite(x$FoldEnrichment[ora])];x$FoldEnrichment[ix]<-(ratio(x$GeneRatio[ix],1)/ratio(x$GeneRatio[ix],2))/(ratio(x$BgRatio[ix],1)/ratio(x$BgRatio[ix],2))
    }
  }
  x
}
review_config <- function(cfg, mode, config_path) {
  defaults <- list(fdr_threshold=.05, top_per_direction=if(mode=='atlas')15L else 6L,
    matrix_per_database=4L, terms_per_page=if(mode=='atlas')30L else 12L,
    max_jaccard=.6, gene_examples=4L, max_genes=40L, rows_per_page=30L,
    display_neighbors=3L, style='dot', positive_label='Numerator',negative_label='Denominator')
  for(k in names(defaults))if(is.null(cfg[[k]]))cfg[[k]]<-defaults[[k]]
  for(k in c('population','comparison_id','output_dir'))
    if(!is.character(cfg[[k]])||length(cfg[[k]])!=1L||is.na(cfg[[k]])||!nzchar(cfg[[k]]))stop('Required scalar: ',k)
  if(is.null(cfg$title))cfg$title<-paste(cfg$population,cfg$comparison_id)
  fmt <- figure_format(cfg)
  if(!is.null(cfg$figure_format)) {
    if(length(cfg$figure_format)!=1L||!tolower(cfg$figure_format)%in%c('pdf','png'))stop('figure_format must be pdf or png')
    legacy <- tolower(cfg$figure_format)
    if((!is.null(cfg[['output']])||!is.null(cfg[['plots']])||!is.null(cfg[['enrichment']])||!is.null(cfg[['visualization']]))&&legacy!=fmt)stop('Conflicting figure formats')
    cfg[["output"]]$figure_format <- legacy
  }
  cfg$figure_format <- figure_format(cfg)
  if(mode!='plot' && cfg$figure_format!='pdf')stop('Atlas, dense and semantic module books require PDF; select PDF for this run')
  if(cfg$figure_format=='pdf'&&!nzchar(Sys.which('gs')))stop('Ghostscript (gs) is required to assemble PDF books')
  for(k in c('top_per_direction','matrix_per_database','terms_per_page','gene_examples','max_genes','rows_per_page','display_neighbors')) {
    z<-cfg[[k]]
    if(!is.numeric(z)||length(z)!=1L||!is.finite(z)||z<1||z!=floor(z))stop('Expected positive integer: ',k)
  }
  if(cfg$terms_per_page>30L||cfg$rows_per_page>30L)stop('Page budget cannot exceed 30 rows')
  for(k in c('fdr_threshold','max_jaccard'))if(!is.numeric(cfg[[k]])||length(cfg[[k]])!=1L||!is.finite(cfg[[k]])||cfg[[k]]<0||cfg[[k]]>1)stop('Invalid ',k)
  if(cfg$fdr_threshold<=0)stop('FDR threshold must be positive')
  if(!cfg$style%in%c('dot','lollipop','both'))stop('Invalid style')
  paths<-c(if(mode=='plot')'input_table' else c('enrichment_table','differential_table','mapping_table'),if(mode=='atlas')'status_table')
  out<-normalizePath(cfg$output_dir,mustWork=FALSE)
  for(k in paths){
    cfg[[k]]<-normalizePath(cfg[[k]],mustWork=TRUE)
    if(cfg[[k]]==out||startsWith(cfg[[k]],paste0(out,'/')))stop('Source tables must be outside output_dir')
  }
  cfg$output_dir<-out
  config_path<-normalizePath(config_path,mustWork=TRUE)
  if(startsWith(config_path,paste0(out,'/')))stop('Keep source config outside output_dir')
  if(!is.null(cfg$similarity_cache))cfg$similarity_cache<-normalizePath(cfg$similarity_cache,mustWork=TRUE)
  attr(cfg,'config_path')<-config_path
  input<-if(mode=='plot')cfg$input_table else cfg$enrichment_table
  x<-review_read(input,cfg)
  required<-c('population','comparison_id','database','method','direction','ID','Description','p.adjust','NES','Count','geneID','core_enrichment')
  if(!all(required%in%names(x)))stop('Missing enrichment columns: ',paste(setdiff(required,names(x)),collapse=', '))
  if(anyNA(x[,c('ID','database','method')])||any(!x$method%in%c('GSEA','ORA')))stop('Invalid enrichment identifiers/method')
  if(anyDuplicated(x[,c('database','method','direction','ID')]))stop('Duplicate pathway results; select one task')
  if(any(is.finite(x$p.adjust)&(x$p.adjust<0|x$p.adjust>1)))stop('FDR must lie in [0,1]')
  if(any(x$method=='GSEA'&!is.finite(x$NES)))stop('GSEA requires finite NES')
  if(any(x$method=='ORA'&(!x$direction%in%c('Up','Down')|is.na(x$direction))))stop('ORA requires Up/Down input direction')
  if(mode!='plot'){
    gsea<-which(x$method=='GSEA');ora<-which(x$method=='ORA')
    if(length(gsea)&&(!'setSize'%in%names(x)||any(!is.finite(x$setSize[gsea])|x$setSize[gsea]<=0)))stop('GSEA requires positive finite setSize')
    if(length(ora)&&mode%in%c('atlas','modules')&&(!'FoldEnrichment'%in%names(x)||any(!is.finite(x$FoldEnrichment[ora])|x$FoldEnrichment[ora]<0)))stop('ORA requires finite FoldEnrichment or original GeneRatio/BgRatio')
    de<-review_read(cfg$differential_table,cfg);mp<-review_read(cfg$mapping_table,cfg)
    if(!all(c('gene','log2FoldChange')%in%names(de))||!all(c('SYMBOL','ENTREZID')%in%names(mp)))stop('Expected DE gene/log2FoldChange and SYMBOL/ENTREZID mapping columns')
    if(anyDuplicated(de$gene))stop('Duplicate DE genes; select one task')
    if(!nrow(de))stop('Selected DE task is empty')
    if(any(x$method=='GSEA')&&!'setSize'%in%names(x))stop('GSEA requires setSize')
    if(mode%in%c('atlas','modules')&&any(x$method=='ORA')&&!'FoldEnrichment'%in%names(x))stop('ORA requires FoldEnrichment')
  }
  if(mode=='atlas'){
    status<-review_read(cfg$status_table,cfg)
    if(!all(c('database','method','status')%in%names(status)))stop('Invalid status table')
  }
  if(is.null(cfg$databases))cfg$databases<-unique(x$database)
  if(mode%in%c('atlas','dense')&&(!length(cfg$databases)||anyNA(cfg$databases)||anyDuplicated(cfg$databases)))stop('Declare unique requested databases, including empty ones')
  if(mode=='modules'){
    if(is.null(cfg$ontology)||!cfg$ontology%in%c('BP','MF','CC')||!identical(cfg$database,paste0('GO_',cfg$ontology)))stop('database must match GO ontology BP/MF/CC')
    if(is.null(cfg$semantic_measure))cfg$semantic_measure<-'Wang'
    if(!identical(cfg$semantic_measure,'Wang'))stop('Only Wang semantic similarity is implemented')
    for(k in c('similarity_threshold','random_seed','orgdb_package'))if(is.null(cfg[[k]]))stop('Explicit config required: ',k)
    if(!is.numeric(cfg$similarity_threshold)||length(cfg$similarity_threshold)!=1||!is.finite(cfg$similarity_threshold)||cfg$similarity_threshold<=0||cfg$similarity_threshold>1)stop('Invalid similarity_threshold')
    if(!is.numeric(cfg$random_seed)||length(cfg$random_seed)!=1||!is.finite(cfg$random_seed)||cfg$random_seed<0||cfg$random_seed>2147483647||cfg$random_seed!=floor(cfg$random_seed))stop('Invalid random_seed')
    if(length(cfg$orgdb_package)!=1||!grepl('^org[.].+[.]eg[.]db$',cfg$orgdb_package)||!requireNamespace(cfg$orgdb_package,quietly=TRUE))stop('Requested species OrgDb package unavailable')
    go<-x[x$database==cfg$database,,drop=FALSE]
    ids<-unique(go$ID)
    if(length(ids)){
      ont<-AnnotationDbi::select(GO.db::GO.db,keys=ids,columns='ONTOLOGY',keytype='GOID')
      if(anyNA(ont$ONTOLOGY)||any(ont$ONTOLOGY!=cfg$ontology))stop('GO IDs do not match requested ontology')
    }
  }
  cfg
}
review_merge_details <- function(out, files, fmt) {
  if(fmt!='pdf')return(files)
  families<-unique(sub('_page[0-9]+[.]pdf$','',files))
  answer<-character()
  for(family in families){
    source<-files[sub('_page[0-9]+[.]pdf$','',files)==family]
    dest<-paste0(family,'.pdf')
    result<-system2('gs',c('-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',shQuote(paste0('-sOutputFile=',file.path(out,dest))),vapply(file.path(out,source),shQuote,character(1))))
    if(result!=0)stop('PDF book merge failed')
    unlink(file.path(out,source));answer<-c(answer,dest)
  }
  answer
}
review_provenance <- function(out,cfg,renderer) {
  path<-file.path(out,'_provenance/manifest.json')
  manifest<-if(file.exists(path))jsonlite::fromJSON(path,simplifyVector=FALSE)else list()
  sources<-unlist(cfg[intersect(c('input_table','enrichment_table','differential_table','mapping_table','status_table','similarity_cache'),names(cfg))])
  manifest$source_md5<-tools::md5sum(sources)
  manifest$renderer_md5<-tools::md5sum(c(renderer,file.path(dirname(renderer),'enrichment_review_common.R'),file.path(dirname(renderer),'figure_output.R')))
  manifest$config<-cfg;manifest$config_md5<-unname(tools::md5sum(attr(cfg,'config_path')))
  manifest$statistics_recomputed<-FALSE
  jsonlite::write_json(manifest,path,pretty=TRUE,auto_unbox=TRUE)
  file.copy(attr(cfg,'config_path'),file.path(out,'_provenance/config.json'),overwrite=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'_provenance/session_info.txt'))
}
review_empty <- function(out,cfg,message) {
  p<-ggplot()+theme_void()+annotate('text',x=0,y=0,label=message)+labs(title=cfg$title,subtitle=paste('Display FDR <=',cfg$fdr_threshold))
  ggsave(file.path(out,paste0('empty_review.',cfg$figure_format)),p,width=8,height=3,bg='white')
  write.table(data.frame(status='empty',reason=message),file.path(out,'display_status.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeLines(c('# Enrichment review','','Status: empty','',message,paste0('[Empty review](empty_review.',cfg$figure_format,')'),'Statistics were not recomputed.'),file.path(out,'RESULTS.md'))
  review_provenance(out,cfg,renderer_path)
  quit(save='no',status=0)
}
review_row_book <- function(out,stem,rows,budget,build) {
  pages<-split(sort(rows,decreasing=TRUE),ceiling(seq_along(rows)/budget));files<-character()
  for(i in seq_along(pages)){
    file<-paste0(stem,'_page',i,'.pdf')
    ggsave(file.path(out,file),build(pages[[i]]),width=14,height=max(4,2.5+.3*length(pages[[i]])),bg='white',limitsize=FALSE)
    files<-c(files,file)
  }
  review_merge_details(out,files,'pdf')
}
