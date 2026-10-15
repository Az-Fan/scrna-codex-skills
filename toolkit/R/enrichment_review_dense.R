#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(jsonlite);library(patchwork)})
args<-commandArgs(TRUE);stopifnot(length(args)==1);cfg<-fromJSON(args[1]);renderer_path<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));source(file.path(dirname(renderer_path),'enrichment_review_common.R'));cfg<-review_config(cfg,'dense',args[1]);out<-cfg$output_dir
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a fresh output directory')
dir.create(file.path(out,'_provenance'),recursive=TRUE,showWarnings=FALSE)
read<-function(p)review_read(p,cfg)
e<-read(cfg$enrichment_table);de<-read(cfg$differential_table);mp<-read(cfg$mapping_table)
e<-e[e$population==cfg$population&e$comparison_id==cfg$comparison_id&e$method=='GSEA'&e$database%in%cfg$databases&is.finite(e$p.adjust)&e$p.adjust<=cfg$fdr_threshold,,drop=FALSE]
de<-de[de$population==cfg$population&de$comparison_id==cfg$comparison_id,,drop=FALSE]
if(!nrow(e))review_empty(out,cfg,'No significant GSEA terms');stopifnot(!anyDuplicated(de$gene))
mp<-merge(mp,de,by.x='SYMBOL',by.y='gene');mp<-mp[order(-abs(mp$log2FoldChange),mp$SYMBOL),];mp<-mp[!duplicated(mp$ENTREZID),]
symbol<-setNames(if('original_gene'%in%names(mp))mp$original_gene else mp$SYMBOL,mp$ENTREZID)
fc<-setNames(mp$log2FoldChange,mp$ENTREZID)
genes<-function(s)if(is.na(s)||!nzchar(s))character()else unique(strsplit(s,'/',fixed=TRUE)[[1]])
e$core_n<-vapply(e$core_enrichment,function(z)length(genes(z)),integer(1));e$dir<-ifelse(e$NES>0,'Up','Down')
sel<-integer();audit<-list();k<-0
for(db in cfg$databases)for(dir in c('Up','Down')){
 idx<-which(e$database==db&e$dir==dir);idx<-idx[order(e$p.adjust[idx],e$ID[idx])];kept<-integer()
 for(i in idx){gi<-genes(e$core_enrichment[i]);js<-vapply(kept,function(j){gj<-genes(e$core_enrichment[j]);if(!length(union(gi,gj)))0 else length(intersect(gi,gj))/length(union(gi,gj))},numeric(1))
 reason<-if(length(js)&&max(js)>cfg$max_jaccard)'redundant'else if(length(kept)>=cfg$top_per_direction)'budget'else'selected'
 k<-k+1;audit[[k]]<-data.frame(database=db,ID=e$ID[i],direction=dir,FDR=e$p.adjust[i],decision=reason)
 if(reason=='selected')kept<-c(kept,i)
 };sel<-c(sel,kept)
}
d<-e[sel,,drop=FALSE];d$row<-rev(seq_len(nrow(d)));n<-nrow(d)
clean<-function(s){s<-gsub('_',' ',sub('^(HALLMARK|REACTOME|KEGG)_','',s));if(s==toupper(s))s<-paste0(toupper(substr(s,1,1)),tolower(substr(s,2,nchar(s))));s<-gsub('nfkb','NF-kB',s,ignore.case=TRUE);s<-gsub('tnfa','TNFa',s,ignore.case=TRUE);paste(strwrap(s,47),collapse='\n')}
d$label<-vapply(d$Description,clean,character(1))
d$top_core<-vapply(seq_len(n),function(i){ids<-genes(d$core_enrichment[i]);ids<-ids[ids%in%names(fc)];ids<-ids[order(-sign(d$NES[i])*fc[ids])];paste(head(symbol[ids],4),collapse=', ')},character(1))
colours<-c(Down='#28648C',Up='#B7453A')
base<-theme_void(base_size=10)+theme(plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(8,8,8,8))
# Fixed column positions; rows carry both statistics and concrete contributing genes.
p<-ggplot(d,aes(y=row))+geom_rect(aes(xmin=0,xmax=15,ymin=row-.48,ymax=row+.48),fill=ifelse(seq_len(n)%%2,'#F2F4F6','white'),colour=NA)+geom_segment(aes(x=0,xend=15,y=row-.5,yend=row-.5),colour='grey88',linewidth=.25)+
 geom_text(aes(x=.08,label=database),hjust=0,size=2.8,colour='#59636B')+
 geom_text(aes(x=1.65,label=label),hjust=0,size=3.1,lineheight=.95)+
 geom_rect(aes(xmin=7.5,xmax=7.5+pmin(3,abs(NES))*.55,ymin=row-.17,ymax=row+.17,fill=dir),show.legend=FALSE)+
 geom_text(aes(x=9.3,label=sprintf('%+.2f',NES),colour=dir),size=3.2,fontface='bold',show.legend=FALSE)+
 geom_text(aes(x=10.15,label=formatC(p.adjust,format='g',digits=2)),size=3.1)+
 geom_text(aes(x=11.0,label=paste0(core_n,'/',setSize)),size=3.1)+
 geom_text(aes(x=11.65,label=top_core),hjust=0,size=3.0)+scale_fill_manual(values=colours)+scale_colour_manual(values=colours)+
 scale_x_continuous(limits=c(0,15),expand=c(0,0))+scale_y_continuous(limits=c(0,n+1),expand=c(0,0))+base+
 labs(title='Pathways and their leading-edge genes',subtitle=cfg$title,caption=paste0('GSEA FDR <= ',cfg$fdr_threshold,'. Red: positive NES; blue: negative NES. Bar lengths clipped at |NES|=3. Bar length: |NES|; Core = leading-edge / tested gene-set size.\nCore examples: up to 4 genes ranked by log2FC in the pathway NES direction; individual genes need not pass DE thresholds.'))
headers<-data.frame(x=c(.08,1.65,7.5,9.3,10.15,11,11.65),label=c('Database','Pathway','|NES|','NES','FDR','Core/set','Leading-edge examples'))
# Repeat column headings on every information-table page.
review_row_book(out,'pathway_information_table',d$row,cfg$rows_per_page,function(rows) p+geom_text(data=headers,aes(x=x,label=label),y=max(rows)+.65,inherit.aes=FALSE,hjust=c(0,0,0,.5,.5,.5,0),fontface='bold',size=3.2)+coord_cartesian(ylim=range(rows)+c(-.6,1.1),expand=FALSE))
# Reserve at least two mapped leading-edge genes for every displayed pathway before shared genes.
cores<-lapply(d$core_enrichment,genes);reserved<-unique(unlist(lapply(seq_len(n),function(i){ids<-intersect(cores[[i]],names(fc));head(ids[order(-sign(d$NES[i])*fc[ids])],2)})))
freq<-table(unlist(cores));candidates<-intersect(names(freq),names(fc));candidates<-candidates[order(-as.numeric(freq[candidates]),-abs(fc[candidates]),candidates)]
if(length(reserved)>cfg$max_genes)stop('max_genes too small for per-pathway coverage; increase it')
ids<-head(unique(c(reserved,candidates)),cfg$max_genes);if(!length(ids))stop('No mapped leading-edge genes for the selected pathways')
# group by logFC direction, then shared membership, so common drivers form visible blocks.
ids<-ids[order(fc[ids]<0,-as.numeric(freq[ids]),-abs(fc[ids]),ids)]
grid<-expand.grid(gene_id=ids,i=seq_len(n),stringsAsFactors=FALSE);grid$row<-d$row[grid$i]
grid$member<-mapply(function(g,i)g%in%cores[[i]],grid$gene_id,grid$i);grid$log2FC<-ifelse(grid$member,fc[grid$gene_id],NA_real_)
grid$x<-match(grid$gene_id,ids);top<-data.frame(x=seq_along(ids),row=1,log2FC=as.numeric(fc[ids]))
scale<-scale_fill_gradient2(low=colours['Down'],mid='#FAFAFA',high=colours['Up'],limits=c(-2,2),oob=scales::squish,na.value='#EEEEEE',name='Gene log2FC\n(clipped at +/-2)')
strip<-ggplot(top,aes(x,row,fill=log2FC))+geom_tile(colour='white',linewidth=.3)+scale+theme_void()+theme(legend.position='none',plot.margin=margin(0,8,0,8))
heat<-ggplot(grid,aes(x,row,fill=log2FC))+geom_tile(colour='white',linewidth=.25)+scale+
 scale_x_continuous(breaks=seq_along(ids),labels=symbol[ids],expand=c(0,0))+
 scale_y_continuous(breaks=d$row,labels=paste0(d$database,' | ',gsub('\n',' ',d$label)),expand=c(0,0))+
 theme_minimal(base_size=10)+theme(panel.grid=element_blank(),axis.title=element_blank(),axis.text.x=element_text(angle=60,hjust=1,vjust=1,colour='#222222',size=9),axis.text.y=element_text(size=9,colour='#222222'),legend.position='bottom',plot.margin=margin(0,8,8,8))
q<-strip/heat+plot_layout(heights=c(.35,n))+plot_annotation(title='Leading-edge genes shared across pathways',subtitle=cfg$title,caption='Coloured cells: gene is in this pathway leading edge; grey: not in leading edge. Top strip: gene log2FC.\nColumns selected for shared membership with up to two available mapped examples per pathway; full gene lists retained. No new statistical test.',theme=theme(plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0)))
review_row_book(out,'pathway_core_gene_matrix',d$row,cfg$rows_per_page,function(rows) strip/(heat+coord_cartesian(ylim=range(rows)+c(-.5,.5),expand=FALSE))+plot_layout(heights=c(.35,length(rows)))+plot_annotation(title='Leading-edge gene membership',subtitle=cfg$title,caption='Colour: member gene log2FC, clipped at +/-2. Grey: not a leading-edge member. No new statistical test.'))
write_tsv<-function(v,name){for(k in names(v))if(is.character(v[[k]]))v[[k]]<-gsub('[\r\n\t]',' ',v[[k]]);write.table(v,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
write_tsv(d,'selected_pathways.tsv');write_tsv(do.call(rbind,audit),'pathway_selection_audit.tsv')
write_tsv(data.frame(gene_id=ids,gene=unname(symbol[ids]),log2FC=unname(fc[ids]),n_pathways=as.integer(freq[ids])),'selected_genes.tsv')
grid$pathway_id<-d$ID[grid$i];write_tsv(grid,'pathway_gene_membership.tsv')
missing<-setdiff(unique(unlist(cores)),names(fc));write_tsv(data.frame(unmapped_core_id=missing),'unmapped_core_ids.tsv')
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
write_json(list(version='dense_v1',config=cfg,input_md5=tools::md5sum(c(cfg$enrichment_table,cfg$differential_table,cfg$mapping_table)),renderer_md5=unname(tools::md5sum(script)),n_pathways=n,n_genes=length(ids),n_unmapped_core_ids=length(missing),statistics_recomputed=FALSE),file.path(out,'_provenance/manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(c('# Dense enrichment preview','','Status: completed','','- [通路信息表图](pathway_information_table.pdf)','- [通路–核心基因矩阵](pathway_core_gene_matrix.pdf)','',paste0('GSEA FDR <= ',cfg$fdr_threshold,'; within-database/direction leading-edge Jaccard suppression for display only; no cross-database merging.'),'信息表中Core/set为leading-edge数/受检基因集大小；例示基因按通路NES方向的log2FC排列。','矩阵彩色格为leading-edge成员，颜色来自基因log2FC；灰色仅表示不在该leading edge，不能解释为未表达。','',paste('Pathways:',n,'; displayed genes:',length(ids),' ; unmapped core IDs:',length(missing)), '', '- [完整入图通路](selected_pathways.tsv)','- [基因及log2FC](selected_genes.tsv)','- [通路成员关系](pathway_gene_membership.tsv)','- [选择记录](pathway_selection_audit.tsv)'),file.path(out,'RESULTS.md'))
cat('DENSE_PREVIEW_COMPLETE',n,'pathways',length(ids),'genes\n')

review_provenance(out,cfg,renderer_path)
