#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(jsonlite);library(GOSemSim);library(igraph)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);cfg<-fromJSON(a[1]);renderer_path<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));source(file.path(dirname(renderer_path),'enrichment_review_common.R'));cfg<-review_config(cfg,'modules',a[1]);out<-cfg$output_dir
stopifnot(cfg$semantic_measure=='Wang',cfg$similarity_threshold>0,cfg$similarity_threshold<=1,cfg$rows_per_page>0,cfg$rows_per_page<=30)
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a fresh output directory')
dir.create(file.path(out,'_provenance'),recursive=TRUE,showWarnings=FALSE)
read<-function(p)review_read(p,cfg)
write_tsv<-function(d,name){for(k in names(d))if(is.character(d[[k]]))d[[k]]<-gsub('[\r\n\t]',' ',d[[k]]);write.table(d,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
x<-read(cfg$enrichment_table);x<-x[x$population==cfg$population&x$comparison_id==cfg$comparison_id&x$database==cfg$database&is.finite(x$p.adjust)&x$p.adjust<=cfg$fdr_threshold,,drop=FALSE]
if(!nrow(x))review_empty(out,cfg,'No significant GO terms');stopifnot(!anyDuplicated(x[,c('method','direction','ID')]))
x$group<-paste(x$method,ifelse(x$method=='GSEA',ifelse(x$NES>=0,'Up','Down'),x$direction),sep='_')
ids<-sort(unique(x$ID));cat('SEMANTIC_SIMILARITY_START',length(ids),'unique GO terms\n');flush.console()
if(!is.null(cfg$similarity_cache)){
 cache<-readRDS(cfg$similarity_cache)
 stopifnot(is.list(cache),identical(cache$ontology,cfg$ontology),identical(cache$orgdb_package,cfg$orgdb_package),identical(cache$orgdb_version,as.character(packageVersion(cfg$orgdb_package))),identical(cache$measure,cfg$semantic_measure),all(ids%in%rownames(cache$similarity)),all(ids%in%colnames(cache$similarity)),identical(cache$go_db_version,as.character(packageVersion('GO.db'))),identical(cache$gosemsim_version,as.character(packageVersion('GOSemSim'))))
 sim<-cache$similarity[ids,ids,drop=FALSE]
}else{
 sem<-GOSemSim::godata(OrgDb=cfg$orgdb_package,ont=cfg$ontology,computeIC=FALSE)
 sim<-GOSemSim::mgoSim(ids,ids,semData=sem,measure='Wang',combine=NULL)
}
if(any(!is.finite(sim)))stop('Unresolved semantic similarity; inspect GO coverage before clustering')
if(!isSymmetric(sim,tol=1e-8)||any(sim<0|sim>1))stop('Invalid semantic similarity matrix');diag(sim)<-1;saveRDS(sim,file.path(out,'_provenance/semantic_similarity.rds'))
saveRDS(list(similarity=sim,orgdb_package=cfg$orgdb_package,orgdb_version=as.character(packageVersion(cfg$orgdb_package)),ontology=cfg$ontology,measure=cfg$semantic_measure,go_db_version=as.character(packageVersion('GO.db')),gosemsim_version=as.character(packageVersion('GOSemSim'))),file.path(out,'_provenance/semantic_cache.rds'))
cat('SEMANTIC_SIMILARITY_COMPLETE\n');flush.console()
de<-read(cfg$differential_table);de<-de[de$population==cfg$population&de$comparison_id==cfg$comparison_id,,drop=FALSE];mp<-read(cfg$mapping_table);mp<-merge(mp,de,by.x='SYMBOL',by.y='gene');mp<-mp[order(-abs(mp$log2FoldChange),mp$SYMBOL),];mp<-mp[!duplicated(mp$ENTREZID),]
sym<-setNames(if('original_gene'%in%names(mp))mp$original_gene else mp$SYMBOL,mp$ENTREZID);sym2<-setNames(if('original_gene'%in%names(de))de$original_gene else de$gene,de$gene)
splitgenes<-function(z)if(is.na(z)||!nzchar(z))character()else unique(strsplit(z,'/',fixed=TRUE)[[1]])
clean<-function(s,w=65)paste(strwrap(gsub('_',' ',s),w),collapse='\n')
allterms<-list();allsums<-list();alledges<-list();netplots<-list();counts<-list();gid<-0
set.seed(cfg$random_seed)
for(group in c('GSEA_Up','GSEA_Down','ORA_Up','ORA_Down')){
 d<-x[x$group==group,,drop=FALSE];d<-d[order(d$ID),];n<-nrow(d);if(!n)next
 s<-sim[d$ID,d$ID,drop=FALSE];adj<-s;adj[adj<cfg$similarity_threshold]<-0;diag(adj)<-0
 g<-igraph::graph_from_adjacency_matrix(adj,mode='undirected',weighted=TRUE,diag=FALSE)
 set.seed(cfg$random_seed)
 if(ecount(g)>0){mem<-membership(cluster_louvain(g,weights=E(g)$weight,resolution=1))}else mem<-seq_len(n)
 # Canonical module IDs: decreasing size, then best term FDR, then GO ID.
 groups<-split(seq_len(n),as.integer(mem));ord<-order(-lengths(groups),vapply(groups,function(i)min(d$p.adjust[i]),numeric(1)),vapply(groups,function(i)min(d$ID[i]),character(1)))
 d$module<-NA_character_;summaries<-list()
 for(j in seq_along(ord)){
  ii<-groups[[ord[j]]];module<-sprintf('M%02d',j);d$module[ii]<-module
  central<-if(length(ii)==1)1 else rowMeans(s[ii,ii,drop=FALSE]);rep_i<-ii[order(-central,d$p.adjust[ii],d$ID[ii])[1]]
  field<-if(startsWith(group,'GSEA'))'core_enrichment'else'geneID';freq<-sort(table(unlist(lapply(d[[field]][ii],splitgenes))),decreasing=TRUE)
  labels<-sym[names(freq)];bad<-is.na(labels);labels[bad]<-sym2[names(freq)[bad]];labels[is.na(labels)]<-names(freq)[is.na(labels)]
  shared<-freq>=2;genes<-if(any(shared))paste(head(paste0(labels[shared],'(',as.integer(freq[shared]),')'),5),collapse=', ')else''
  summaries[[j]]<-data.frame(group=group,module=module,n_terms=length(ii),representative_ID=d$ID[rep_i],representative=d$Description[rep_i],best_FDR=min(d$p.adjust[ii]),median_NES=if(startsWith(group,'GSEA'))median(d$NES[ii])else NA_real_,min_NES=if(startsWith(group,'GSEA'))min(d$NES[ii])else NA_real_,max_NES=if(startsWith(group,'GSEA'))max(d$NES[ii])else NA_real_,shared_genes=if(nzchar(genes))genes else '(no repeated contribution genes)',n_semantic_edges=sum(adj[ii,ii,drop=FALSE]>0)/2)
 }
 su<-do.call(rbind,summaries);d<-d[order(d$module,d$p.adjust,d$ID),];gids<-V(g)$name
 # Network positions use the full semantic graph. Sparsify edges for display only.
 set.seed(cfg$random_seed);xy<-if(n==1)matrix(c(0,0),ncol=2)else layout_with_fr(g,weights=if(ecount(g))E(g)$weight else NULL,niter=1000)
 nodes<-data.frame(ID=gids,x=xy[,1],y=xy[,2]);nodes<-merge(nodes,d,by='ID',sort=FALSE)
 edges<-if(ecount(g))as_data_frame(g,what='edges')else data.frame(from=character(),to=character(),weight=numeric())
 if(nrow(edges)){
  edges$edge_id<-seq_len(nrow(edges));chosen<-integer()
  for(id in gids){ix<-which(edges$from==id|edges$to==id);ix<-ix[order(-edges$weight[ix],edges$from[ix],edges$to[ix])];chosen<-c(chosen,head(ix,cfg$display_neighbors))}
  edges$displayed<-edges$edge_id%in%unique(chosen);ed<-edges[edges$displayed,,drop=FALSE]
  ed$x<-nodes$x[match(ed$from,nodes$ID)];ed$y<-nodes$y[match(ed$from,nodes$ID)];ed$xend<-nodes$x[match(ed$to,nodes$ID)];ed$yend<-nodes$y[match(ed$to,nodes$ID)]
 }else {edges$displayed<-logical();ed<-data.frame(x=numeric(),y=numeric(),xend=numeric(),yend=numeric())}
 edges$group<-rep(group,nrow(edges));gid<-gid+1;alledges[[gid]]<-edges
 colors<-umap_palette(su$module)
 write_tsv(data.frame(module=names(colors),colour=unname(colors)),paste0(tolower(group),'_module_colours.tsv'))
 nodes$score<-pmin(10,-log10(pmax(nodes$p.adjust,1e-300)))
 centers<-do.call(rbind,lapply(su$module,function(m){v<-nodes[nodes$module==m,];data.frame(module=m,x=median(v$x),y=median(v$y),n=nrow(v))}))
 # Label multigene modules only; singletons are retained as nodes with explicit counts.
 centers<-centers[centers$n>=2,,drop=FALSE];centers$representative<-su$representative[match(centers$module,su$module)]
 centers$label<-if(!nrow(centers))character()else paste0(centers$module,' (',centers$n,')\n',vapply(centers$representative,function(z){if(nchar(z)>85)z<-paste0(substr(z,1,82),'...');clean(z,30)},character(1)))
 p<-ggplot()+geom_segment(data=ed,aes(x,y,xend=xend,yend=yend),colour='grey65',linewidth=.25,alpha=.35)+geom_point(data=nodes,aes(x,y,fill=module,size=score),shape=21,colour='white',stroke=.25)+scale_fill_manual(values=colors,guide='none')+scale_size_continuous(range=c(1.5,4),name='-log10 FDR\n(capped at 10)')+guides(size=guide_legend(override.aes=list(fill='grey50',colour='grey50')))+theme_void(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(10,10,10,10))+
  labs(title=paste(cfg$database,group,'semantic modules'),subtitle=paste0(cfg$title,'\nAll ',n,' significant terms | ',nrow(su),' modules (',sum(su$n_terms==1),' singleton modules) | Wang >= ',cfg$similarity_threshold,'; Louvain'),caption=paste0('Every significant term is retained. Labels: module ID (term count) and representative GO term; full names in module overview.\nClustering uses all ',ecount(g),' thresholded semantic edges. Display keeps strongest ',cfg$display_neighbors,' neighbors/term (',nrow(ed),' edges); layout uses full graph.\nSpatial distance is a layout aid, not an effect size. Modules are exploratory and have no new P value.'))
 if(nrow(centers)){if(requireNamespace('ggrepel',quietly=TRUE))p<-p+ggrepel::geom_label_repel(data=centers,aes(x,y,label=label),size=3,seed=cfg$random_seed,max.overlaps=Inf,box.padding=.5,min.segment.length=0)else p<-p+geom_label(data=centers,aes(x,y,label=label),size=3)}
 netplots[[group]]<-p;allterms[[group]]<-d;allsums[[group]]<-su
 counts[[group]]<-data.frame(group=group,n_significant=n,n_modules=nrow(su),n_singleton_modules=sum(su$n_terms==1),n_terms_in_multiterm_modules=sum(su$n_terms[su$n_terms>=2]),n_semantic_edges=ecount(g),n_display_edges=nrow(ed))
 write_tsv(nodes,paste0(tolower(group),'_network_nodes.tsv'))
 cat('GROUP_COMPLETE',group,n,'terms',nrow(su),'modules\n');flush.console()
}
terms<-do.call(rbind,allterms);summary<-do.call(rbind,allsums);count<-do.call(rbind,counts);edges<-do.call(rbind,alledges)
write_tsv(terms,'term_module_membership.tsv');write_tsv(summary,'module_summary.tsv');write_tsv(count,'group_counts.tsv');write_tsv(edges,'semantic_edges.tsv')
pdf(file.path(out,'semantic_networks.pdf'),width=13,height=9,useDingbats=FALSE);for(p in netplots)print(p);dev.off()
# Overview and complete term atlas: every module, every significant term; no Top truncation.
makebook<-function(kind){files<-character()
 for(group in names(allsums)){
  z<-if(kind=='overview')allsums[[group]]else allterms[[group]];np<-ceiling(nrow(z)/cfg$rows_per_page)
  for(pg in seq_len(np)){
   v<-z[seq.int((pg-1)*cfg$rows_per_page+1,min(pg*cfg$rows_per_page,nrow(z))),,drop=FALSE];nr<-nrow(v);v$row<-rev(seq_len(nr))
   if(kind=='overview'){v$label<-vapply(v$representative,clean,character(1));v$number<-as.character(v$n_terms);v$fdr<-v$best_FDR;v$gene_text<-v$shared_genes;v$effect_text<-if(startsWith(group,'GSEA'))sprintf('%.2f',v$median_NES)else'-'}else{
    v$label<-vapply(v$Description,clean,character(1));v$number<-if(startsWith(group,'GSEA'))as.character(vapply(v$core_enrichment,function(s)length(splitgenes(s)),integer(1)))else as.character(v$Count);v$fdr<-v$p.adjust;v$effect_text<-if(startsWith(group,'GSEA'))sprintf('%+.2f',v$NES)else sprintf('%.2f',v$FoldEnrichment)
    v$gene_text<-vapply(if(startsWith(group,'GSEA'))v$core_enrichment else v$geneID,function(g){ids<-splitgenes(g);labels<-sym[ids];bad<-is.na(labels);labels[bad]<-sym2[ids[bad]];labels[is.na(labels)]<-ids[is.na(labels)];paste(head(labels,5),collapse=', ')},character(1))
   }
   p<-ggplot(v,aes(y=row))+geom_rect(aes(xmin=0,xmax=16,ymin=row-.48,ymax=row+.48),fill=ifelse(seq_len(nr)%%2,'#F1F4F6','white'),colour=NA)+geom_text(aes(x=.1,label=module),hjust=0,size=3)+geom_text(aes(x=.95,label=number),size=3.1,fontface='bold')+geom_text(aes(x=1.5,label=label),hjust=0,size=3.0,lineheight=.94)+geom_text(aes(x=8.8,label=effect_text),size=3.1)+geom_text(aes(x=9.75,label=formatC(fdr,format='g',digits=2)),size=3.0)+geom_text(aes(x=10.6,label=gene_text),hjust=0,size=2.8)+scale_x_continuous(limits=c(0,16),expand=c(0,0))+scale_y_continuous(limits=c(0,nr+1),expand=c(0,0))+theme_void()+theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(10,10,10,10))+
    labs(title=paste(cfg$database,group,if(kind=='overview')'module overview'else'all terms by module'),subtitle=paste0(cfg$title,'\n',count$n_significant[count$group==group],' significant terms; ',count$n_modules[count$group==group],' modules; no Top truncation | page ',pg,'/',np,' displays ',nr,'/',nrow(z),if(kind=='overview')' modules'else' terms'),caption=if(kind=='overview')'Representative = semantic medoid, not a manually curated theme. Shared genes: number of module terms containing that contribution gene.\nMedian NES is descriptive, best FDR belongs to one term; modules do not have a new significance test.'else'All significant terms are retained. GSEA core and ORA hit examples use contribution genes, not full GO annotations.\nFE is ORA fold enrichment; NES is GSEA signed enrichment. Module membership derived from GO Wang semantics, not those example genes.')
   h<-data.frame(x=c(.1,.95,1.5,8.8,9.75,10.6),label=c('Module',if(kind=='overview')'Terms'else'Core/hit',if(kind=='overview')'Representative GO term'else'GO term',if(kind=='overview')'Median NES'else if(startsWith(group,'GSEA'))'NES'else'FE',if(kind=='overview')'Best FDR'else'FDR',if(kind=='overview')'Shared contribution genes (term count)'else'Contribution gene examples'))
   p<-p+geom_text(data=h,aes(x=x,y=nr+.6,label=label),hjust=c(0,.5,0,.5,.5,0),inherit.aes=FALSE,size=3,fontface='bold')
   f<-tempfile(fileext='.pdf');ggsave(f,p,width=14,height=max(3.2,2+.29*nr),device=grDevices::pdf,useDingbats=FALSE,bg='white');files<-c(files,f)
  }
 }
 dest<-file.path(out,if(kind=='overview')'module_overview.pdf'else'module_term_atlas.pdf');code<-system2('gs',c('-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',shQuote(paste0('-sOutputFile=',dest)),vapply(files,shQuote,character(1))));if(code!=0)stop('PDF merge failed');unlink(files)
}
makebook('overview');makebook('terms')
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));write_json(list(version='semantic_modules_v1',config=cfg,source_md5=tools::md5sum(c(cfg$enrichment_table,cfg$mapping_table,cfg$differential_table)),renderer_md5=unname(tools::md5sum(script)),n_significant_rows=nrow(terms),group_counts=count,seed=cfg$random_seed,similarity='GO Wang semantics',clustering='weighted Louvain, resolution=1',module_labels='semantic medoid representatives, not curated themes',statistics_recomputed=FALSE),file.path(out,'_provenance/manifest.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'_provenance/session_info.txt'))
writeLines(c(paste('#',cfg$database,'semantic modules'),'','Status: completed','','- [模块总览](module_overview.pdf)','- [全部显著条目网络](semantic_networks.pdf)','- [按模块组织的完整条目图册](module_term_atlas.pdf)','- [四组计数](group_counts.tsv)','- [模块摘要与代表条目](module_summary.tsv)','- [全部条目归属](term_module_membership.tsv)','',paste0(cfg$database,', FDR <= ',cfg$fdr_threshold,'; Wang >= ',cfg$similarity_threshold,'; weighted Louvain, seed=',cfg$random_seed,'. Methods and directions clustered separately; all significant terms retained.'),paste0('Display retains the union of strongest ',cfg$display_neighbors,' neighbors per node. Clustering/layout use the full graph; singletons retained; module IDs are local to each group.'),'Module names are semantic medoid terms, not curated themes. No module P values or stability claims are produced; threshold sensitivity requires a separate explicit comparison.','共享贡献基因用于解释，不参与语义聚类；括号内数字表示该模块多少个条目包含此贡献基因。','所有原始NES、FDR、差异结果不改变。'),file.path(out,'RESULTS.md'))
cat('SEMANTIC_MODULES_COMPLETE',nrow(terms),'terms',nrow(summary),'modules\n')

review_provenance(out,cfg,renderer_path)
