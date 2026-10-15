#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);cfg<-fromJSON(a[1]);renderer_path<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));source(file.path(dirname(renderer_path),'enrichment_review_common.R'));cfg<-review_config(cfg,'atlas',a[1]);out<-cfg$output_dir
stopifnot(cfg$top_per_direction>0,cfg$terms_per_page>0,cfg$terms_per_page<=30,cfg$fdr_threshold>0,cfg$fdr_threshold<=1)
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a fresh output directory')
dir.create(file.path(out,'_provenance'),recursive=TRUE,showWarnings=FALSE)
read<-function(p)review_read(p,cfg)
x<-read(cfg$enrichment_table);status<-read(cfg$status_table);de<-read(cfg$differential_table);mp<-read(cfg$mapping_table)
x<-x[x$population==cfg$population&x$comparison_id==cfg$comparison_id,,drop=FALSE]
de<-de[de$population==cfg$population&de$comparison_id==cfg$comparison_id,,drop=FALSE]
stopifnot(!anyDuplicated(x[,c('database','method','direction','ID')]),!anyDuplicated(de$gene))
mp<-merge(mp,de,by.x='SYMBOL',by.y='gene');mp<-mp[order(-abs(mp$log2FoldChange),mp$SYMBOL),];mp<-mp[!duplicated(mp$ENTREZID),]
sym<-setNames(if('original_gene'%in%names(mp))mp$original_gene else mp$SYMBOL,mp$ENTREZID);fc<-setNames(mp$log2FoldChange,mp$ENTREZID)
# ORA may supply readable symbols whereas GSEA leading-edge IDs remain Entrez.
fc_symbol<-setNames(de$log2FoldChange,de$gene);display_symbol<-setNames(if('original_gene'%in%names(de))de$original_gene else de$gene,de$gene)
splitgenes<-function(s){if(is.na(s)||!nzchar(s))character()else unique(strsplit(s,'/',fixed=TRUE)[[1]])}
x$plot_direction<-ifelse(x$method=='GSEA',ifelse(x$NES>=0,'Up','Down'),x$direction)
x$significant<-is.finite(x$p.adjust)&x$p.adjust<=cfg$fdr_threshold
x$display_rank<-rep(NA_integer_,nrow(x))
selected<-list();audit<-list();summary<-list();kk<-0;ss<-0;dd<-0
for(m in c('GSEA','ORA'))for(db in cfg$databases)for(dir in c('Up','Down')){
 idx<-which(x$method==m&x$database==db&x$plot_direction==dir);eligible<-idx[x$significant[idx]];eligible<-eligible[order(x$p.adjust[eligible],x$ID[eligible])];kept<-integer();field<-if(m=='GSEA')'core_enrichment'else'geneID'
 for(i in eligible){gi<-splitgenes(x[[field]][i]);j<-vapply(kept,function(v){gj<-splitgenes(x[[field]][v]);u<-union(gi,gj);if(!length(u))0 else length(intersect(gi,gj))/length(u)},numeric(1));redundant<-length(j)>0&&max(j)>cfg$max_jaccard
 reason<-if(redundant)'redundant'else if(length(kept)>=cfg$top_per_direction)'top_budget'else'selected'
 kk<-kk+1;audit[[kk]]<-data.frame(database=db,method=m,direction=dir,ID=x$ID[i],FDR=x$p.adjust[i],significant_rank=match(i,eligible),decision=reason,redundant_with=if(redundant)x$ID[kept[which.max(j)]]else NA_character_)
 if(reason=='selected'){kept<-c(kept,i);x$display_rank[i]<-length(kept)}
 }
 if(length(kept)){dd<-dd+1;selected[[dd]]<-x[kept,,drop=FALSE]}
 ss<-ss+1;summary[[ss]]<-data.frame(database=db,method=m,direction=dir,n_returned=length(idx),n_finite_fdr=sum(is.finite(x$p.adjust[idx])),n_significant=length(eligible),top_per_direction=cfg$top_per_direction,n_selected=length(kept),n_not_selected=length(eligible)-length(kept),sort_by='FDR ascending, then pathway ID',jaccard_threshold=cfg$max_jaccard)
}
sel<-if(length(selected))do.call(rbind,selected)else x[FALSE,,drop=FALSE];sm<-do.call(rbind,summary)
clean<-function(s){s<-gsub('_',' ',sub('^(HALLMARK|REACTOME|KEGG)_','',s));if(s==toupper(s))s<-paste0(toupper(substr(s,1,1)),tolower(substr(s,2,nchar(s))));for(z in c('DNA','RNA','ATP','MHC','TNFA','NFKB'))s<-gsub(paste0('\\b',tolower(z),'\\b'),if(z=='NFKB')'NF-kB'else if(z=='TNFA')'TNFa'else z,s,ignore.case=TRUE,perl=TRUE);paste(strwrap(s,57),collapse='\n')}
# Explicitly sort by descending effect for Up and ascending effect for Down.
examples<-function(r){ids<-splitgenes(if(r$method=='GSEA')r$core_enrichment else r$geneID);vals<-fc[ids];labels<-sym[ids];unknown<-is.na(vals);vals[unknown]<-fc_symbol[ids[unknown]];labels[unknown]<-display_symbol[ids[unknown]];ok<-is.finite(vals)&!is.na(labels);ids<-ids[ok];vals<-vals[ok];labels<-labels[ok];ord<-order(if(r$plot_direction=='Up')-vals else vals,ids);paste(head(labels[ord],cfg$gene_examples),collapse=', ')}
if(nrow(sel)){sel$label<-vapply(sel$Description,clean,character(1));sel$example_genes<-vapply(seq_len(nrow(sel)),function(i)examples(sel[i,]),character(1))}
write_tsv<-function(d,name){for(k in names(d))if(is.character(d[[k]]))d[[k]]<-gsub('[\r\n\t]',' ',d[[k]]);write.table(d,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
write_tsv(sm,'enrichment_display_summary.tsv');write_tsv(sel,'displayed_terms.tsv');write_tsv(if(length(audit))do.call(rbind,audit)else data.frame(database=character(),method=character(),direction=character(),ID=character(),FDR=numeric(),decision=character()),'selection_audit.tsv')
colours<-c(Up='#B44337',Down='#28668F');page_audit<-list();pi<-0
for(m in c('GSEA','ORA')){
 page_files<-character()
 emit<-function(p,height){f<-tempfile(paste0('atlas_',m,'_'),fileext='.pdf');ggsave(f,p,width=14,height=height,units='in',device=grDevices::pdf,useDingbats=FALSE,bg='white',limitsize=FALSE);page_files<<-c(page_files,f)}
 for(db in cfg$databases){
  z<-sel[sel$method==m&sel$database==db,,drop=FALSE];s<-sm[sm$method==m&sm$database==db,];total<-sum(s$n_returned);nf<-sum(s$n_finite_fdr);sig<-sum(s$n_significant);count<-nrow(z)
  npage<-max(1,ceiling(count/cfg$terms_per_page));sigup<-sum(s$n_significant[s$direction=='Up']);sigdn<-sum(s$n_significant[s$direction=='Down']);selup<-sum(s$n_selected[s$direction=='Up']);seldn<-sum(s$n_selected[s$direction=='Down'])
  dbstatus<-status[status$database==db&status$method%in%c(m,'ALL'),,drop=FALSE];status_text<-if(nrow(dbstatus))paste(unique(dbstatus$status),collapse=', ')else'no status record'
  for(pg in seq_len(npage)){
   v<-if(count)z[seq.int((pg-1)*cfg$terms_per_page+1,min(pg*cfg$terms_per_page,count)),,drop=FALSE]else z
   nr<-nrow(v);pi<-pi+1;page_audit[[pi]]<-data.frame(method=m,database=db,page=pg,pages=npage,n_page=nr,n_page_up=sum(v$plot_direction=='Up'),n_page_down=sum(v$plot_direction=='Down'),n_selected=count,n_significant=sig,n_returned=total,n_finite_fdr=nf,database_status=status_text)
   header<-paste0('Returned ',total,' (Up ',s$n_returned[s$direction=='Up'],' / Down ',s$n_returned[s$direction=='Down'],') | finite FDR ',nf,' | significant ',sig,' (',if(m=='GSEA')'NES+ 'else'Up ',sigup,' / ',if(m=='GSEA')'NES- 'else'Down ',sigdn,') | FDR <= ',cfg$fdr_threshold)
   subheader<-paste0('Top ',cfg$top_per_direction,' per direction AFTER redundancy filtering; FDR ascending then ID | selected ',count,' (Up ',selup,' / Down ',seldn,') | page ',pg,'/',npage,' shows ',nr,' (Up ',sum(v$plot_direction=='Up'),' / Down ',sum(v$plot_direction=='Down'),')')
   if(!nr){
    p<-ggplot()+theme_void()+annotate('text',x=0,y=.5,label=if(any(!dbstatus$status%in%c('completed','empty'))||!nrow(dbstatus))paste('Unavailable or incomplete:',status_text)else if(sig==0)'No significant pathways to display'else'No selected pathways',size=5)+labs(title=paste(db,m),subtitle=paste(cfg$title,header,subheader,paste('Database status:',status_text),sep='\n'))+theme(plot.title=element_text(face='bold',size=16),plot.subtitle=element_text(size=11,lineheight=1.4),plot.margin=margin(15,15,15,15));emit(p,3.5);next
   }
   v$row<-rev(seq_len(nr));v$direction_label<-if(m=='GSEA')ifelse(v$plot_direction=='Up',cfg$positive_label,cfg$negative_label)else paste(v$plot_direction,'genes')
   v$effect<-if(m=='GSEA')v$NES else v$FoldEnrichment
   v$effect_text<-if(m=='GSEA')sprintf('%+.2f',v$effect)else sprintf('%.2f',v$effect)
   v$ratio<-if(m=='GSEA')paste0(vapply(v$core_enrichment,function(g)length(splitgenes(g)),integer(1)),'/',v$setSize)else v$GeneRatio
   # Each PDF page has a content-sized height; sparse databases do not inherit a blank full sheet.
   v$bar<-pmin(3,abs(v$effect))/3*1.0
   if(m=='ORA'){maxeff<-max(z$FoldEnrichment[is.finite(z$FoldEnrichment)]);v$bar<-abs(v$effect)/maxeff*1.0}
   p<-ggplot(v,aes(y=row))+geom_rect(aes(xmin=0,xmax=15.8,ymin=row-.48,ymax=row+.48),fill=ifelse(seq_len(nr)%%2,'#F1F4F6','white'),colour=NA)+geom_segment(aes(x=0,xend=15.8,y=row-.5,yend=row-.5),colour='grey88',linewidth=.2)+
    geom_text(aes(x=.1,label=direction_label,colour=plot_direction),hjust=0,size=2.9,fontface='bold',show.legend=FALSE)+
    geom_text(aes(x=1.35,label=display_rank),size=2.9,colour='grey35')+
    geom_text(aes(x=1.7,label=label),hjust=0,size=3.05,lineheight=.94)+
    geom_rect(aes(xmin=8.1,xmax=8.1+bar,ymin=row-.14,ymax=row+.14,fill=plot_direction),show.legend=FALSE)+
    geom_text(aes(x=9.45,label=effect_text,colour=plot_direction),fontface='bold',size=3.1,show.legend=FALSE)+
    geom_text(aes(x=10.25,label=formatC(p.adjust,format='g',digits=2)),size=3.0)+
    geom_text(aes(x=11.05,label=ratio),size=3.0)+
    geom_text(aes(x=11.65,label=example_genes),hjust=0,size=2.85)+scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
    scale_x_continuous(limits=c(0,15.8),expand=c(0,0))+scale_y_continuous(limits=c(0,nr+1.1),expand=c(0,0))+theme_void(base_size=11)+
    theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10,lineheight=1.25),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(10,10,10,10))+
    labs(title=paste(db,m),subtitle=paste(cfg$title,header,subheader,paste('Database status:',status_text),sep='\n'),caption=paste0(if(m=='GSEA')'Core/set = leading-edge genes / tested set size; signed NES, red positive and blue negative condition. Bar length = min(|NES|, 3).'else'Hit/input = pathway hits / database-annotated foreground. FE = fold enrichment; bars share one scale within this database. Up/Down describe input DE genes, not pathway activation.','\nExamples: up to ',cfg$gene_examples,' core/hit genes ranked by directional log2FC. Jaccard > ',cfg$max_jaccard,' suppressed for display only. Significant totals BEFORE redundancy filtering.\nReturned/finite-FDR counts refer to supplied result rows after enrichment eligibility rules, not every database term. Complete results and selection audit retained.'))
   h<-data.frame(x=c(.1,1.35,1.7,8.1,9.45,10.25,11.05,11.65),label=c('Direction','Rank','Pathway',if(m=='GSEA')'|NES|'else'FE',if(m=='GSEA')'NES'else'FE','FDR',if(m=='GSEA')'Core/set'else'Hit/input',if(m=='GSEA')'Leading-edge examples'else'Hit gene examples'))
   p<-p+geom_text(data=h,aes(x=x,y=nr+.65,label=label),inherit.aes=FALSE,hjust=c(0,.5,0,0,.5,.5,.5,0),size=3,fontface='bold')
   # Render on a content-sized viewport at the top of the standard page.
   emit(p,min(11,max(3.1,2.3+.275*nr)))
  }
 }
 merged<-system2('gs',c('-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',shQuote(paste0('-sOutputFile=',file.path(out,paste0(tolower(m),'_database_atlas.pdf')))),vapply(page_files,shQuote,character(1))))
 if(merged!=0)stop('PDF atlas merge failed');unlink(page_files)
}
write_tsv(do.call(rbind,page_audit),'page_summary.tsv')
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));write_json(list(version='database_atlas_v1',config=cfg,source_md5=tools::md5sum(c(cfg$enrichment_table,cfg$status_table,cfg$mapping_table,cfg$differential_table)),renderer_md5=unname(tools::md5sum(script)),statistics_recomputed=FALSE),file.path(out,'_provenance/manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(c('# Database enrichment atlases','','Status: completed','','- [GSEA：全部数据库](gsea_database_atlas.pdf)','- [ORA：全部数据库](ora_database_atlas.pdf)','- [各数据库方向统计](enrichment_display_summary.tsv)','- [逐页展示计数](page_summary.tsv)','- [全部入图通路与数值](displayed_terms.tsv)','- [显著通路选取及去冗余记录](selection_audit.tsv)','',paste0('Each requested database is retained, including zero-significant/failure pages. FDR <= ',cfg$fdr_threshold,'; top ',cfg$top_per_direction,' per direction after display-only redundancy filtering. Significant totals precede filtering.'),'GSEA正/负方向按NES；ORA按上/下调输入基因集分别计数。返回条目数是完整富集表包含的结果行，不代表数据库中所有通路数。','多页标注当前页/总页数、本页展示数、该数据库选取总数。条目数、显著数、选取数均可追溯至TSV。','GSEA核心数/集合大小与ORA命中数/前景基因数分开解释；例示基因不一定单独显著。'),file.path(out,'RESULTS.md'))
cat('ATLAS_COMPLETE',nrow(sel),'selected rows',pi,'pages\n')

review_provenance(out,cfg,renderer_path)
