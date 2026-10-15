# Exercise installed table-only rendering, including input identity and empty results.
args<-commandArgs(TRUE);scripts<-if(length(args))args[1]else'toolkit/R'
suppressPackageStartupMessages(library(jsonlite))
root<-tempfile('enrichment_review_',tmpdir='/tmp');dir.create(root)
x<-expand.grid(method=c('GSEA','ORA'),direction=c('Up','Down'),j=1:3,stringsAsFactors=FALSE)
x$population<-'EC';x$comparison_id<-'treated_vs_control';x$database<-'GO_BP'
ids<-c('GO:0006954','GO:0006955','GO:0009615');x$ID<-ids[x$j];x$Description<-c('Inflammatory response','Immune response','Response to virus')[x$j]
x$NES<-ifelse(x$method=='GSEA',ifelse(x$direction=='Up',2,-2),NA_real_)
# Each GSEA pathway has one direction only.
x<-x[!(x$method=='GSEA'&((x$j==3&x$direction=='Up')|(x$j!=3&x$direction=='Down'))),]
x$p.adjust<-ifelse(x$j==3,.8,.01);x$setSize<-20;x$Count<-2;x$RichFactor<-.1;x$FoldEnrichment<-2;x$GeneRatio<-'2/10';x$BgRatio<-'20/200';x$geneID<-c('A/B');x$core_enrichment<-'1/2';x$task_id<-'table1_EC'
write.table(x,file.path(root,'enrichment.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
de<-data.frame(gene=c('A','B'),log2FoldChange=c(2,-1),population='EC',comparison_id='treated_vs_control',task_id='table1_EC')
write.table(de,file.path(root,'de.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(SYMBOL=c('A','B'),ENTREZID=c('1','2')),file.path(root,'mapping.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
st<-expand.grid(database=c('GO_BP','KEGG','REACTOME'),method=c('ORA','GSEA'));st$status<-ifelse(st$database=='REACTOME','failed','completed')
write.table(st,file.path(root,'status.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
base<-list(population='EC',comparison_id='treated_vs_control',input_table=file.path(root,'enrichment.tsv'),enrichment_table=file.path(root,'enrichment.tsv'),differential_table=file.path(root,'de.tsv'),mapping_table=file.path(root,'mapping.tsv'),status_table=file.path(root,'status.tsv'),databases=c('GO_BP','KEGG','REACTOME'),title='Treatment versus control',positive_label='Treated',negative_label='Control',top_per_direction=2,terms_per_page=2,rows_per_page=2,max_genes=10,max_jaccard=.6)
source_hash<-tools::md5sum(list.files(root,full.names=TRUE))
run<-function(mode,extra=list(),success=TRUE){
 cfg<-modifyList(base,extra);cfg$output_dir<-file.path(root,paste0(mode,'_',length(list.files(root))))
 config<-paste0(cfg$output_dir,'.json');write_json(cfg,config,auto_unbox=TRUE,pretty=TRUE)
 log<-paste0(cfg$output_dir,'.log')
 status<-system2(file.path(R.home('bin'),'Rscript'),c(file.path(scripts,paste0('enrichment_review_',mode,'.R')),config),stdout=log,stderr=log)
 if(success&&status!=0)stop(paste(readLines(log),collapse='\n'))
 if(!success){stopifnot(status!=0);return(invisible(cfg$output_dir))}
 stopifnot(file.exists(file.path(cfg$output_dir,'RESULTS.md')),file.exists(file.path(cfg$output_dir,'_provenance/manifest.json')))
 # Every result link resolves within output.
 lines<-paste(readLines(file.path(cfg$output_dir,'RESULTS.md')),collapse='\n');links<-regmatches(lines,gregexpr('\\]\\([^)]*\\)',lines))[[1]];links<-substring(links,3,nchar(links)-1);stopifnot(all(file.exists(file.path(cfg$output_dir,links))))
 cfg$output_dir
}
a<-run('atlas');sm<-read.delim(file.path(a,'enrichment_display_summary.tsv'));pg<-read.delim(file.path(a,'page_summary.tsv'))
stopifnot(sum(sm$n_returned)==nrow(x),sum(sm$n_significant)==sum(x$p.adjust<=.05),sum(sm$n_selected)==sum(pg$n_page),all(c('n_page_up','n_page_down')%in%names(pg)),length(list.files(a,pattern='pdf$'))==2,!length(list.files(a,pattern='png$')))
p<-run('plot');stopifnot(file.exists(file.path(p,'evidence_matrix.pdf')),!any(grepl('_page',list.files(p,pattern='pdf$'))))
pp<-run('plot',list(figure_format='png'));stopifnot(length(list.files(pp,pattern='png$'))>0,!length(list.files(pp,pattern='pdf$')))
d<-run('dense');stopifnot(length(list.files(d,pattern='pdf$'))==2)
m<-run('modules',list(database='GO_BP',ontology='BP',orgdb_package='org.Mm.eg.db',semantic_measure='Wang',similarity_threshold=.5,random_seed=14,display_neighbors=3));counts<-read.delim(file.path(m,'group_counts.tsv'));terms<-read.delim(file.path(m,'term_module_membership.tsv'));stopifnot(sum(counts$n_significant)==sum(x$p.adjust<=.05),nrow(terms)==sum(x$p.adjust<=.05))
run('modules',list(database='GO_BP',ontology='MF',orgdb_package='org.Mm.eg.db',similarity_threshold=.5,random_seed=14),FALSE)
run('atlas',list(figure_format='png'),FALSE)
run('plot',list(terms_per_page=1.5),FALSE)
# Duplicate input tasks cannot be collapsed silently.
x2<-rbind(x,transform(x,task_id='table2_EC'));write.table(x2,file.path(root,'ambiguous.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
run('plot',list(input_table=file.path(root,'ambiguous.tsv')),FALSE)
run('plot',list(input_table=file.path(root,'ambiguous.tsv'),task_id='table1_EC'))
run('atlas',list(fdr_threshold=.001));run('dense',list(fdr_threshold=.001));run('modules',list(database='GO_BP',ontology='BP',orgdb_package='org.Mm.eg.db',similarity_threshold=.5,random_seed=14,fdr_threshold=.001))
# Older ORA exports: derive display-only fold enrichment from original ratios.
ora<-x[x$method=='ORA',];ora$FoldEnrichment<-NULL;ora$RichFactor<-NULL
write.table(ora,file.path(root,'ora_only.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
o<-run('atlas',list(enrichment_table=file.path(root,'ora_only.tsv')))
ov<-read.delim(file.path(o,'displayed_terms.tsv'));stopifnot(all(ov$FoldEnrichment==2),all(ov$RichFactor==.1))
# Source-table disambiguation uses input_table as a data column, not the config file path.
x3<-rbind(transform(x,input_table='source_A'),transform(x,input_table='source_B'))
write.table(x3,file.path(root,'sources.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
run('plot',list(input_table=file.path(root,'sources.tsv'),source_table='source_A'))
# Cache binding and source hashes are preserved on reuse.
run('modules',list(database='GO_BP',ontology='BP',orgdb_package='org.Mm.eg.db',semantic_measure='Wang',similarity_threshold=.5,random_seed=14,similarity_cache=file.path(m,'_provenance/semantic_cache.rds')))
stopifnot(identical(unname(source_hash),unname(tools::md5sum(names(source_hash)))))
cat('PASS: installed enrichment review modes, counts, formats, task scope, empty results and source immutability\nArtifacts:',root,'\n')
