#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
args<-commandArgs(TRUE);stopifnot(length(args)==1)
renderer_path<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
cfg<-fromJSON(args[1]);get<-function(k,d) if(is.null(cfg[[k]]))d else cfg[[k]]
source(file.path(dirname(renderer_path),'enrichment_review_common.R'));cfg<-review_config(cfg,'plot',args[1]);input<-cfg$input_table;out<-cfg$output_dir
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a fresh output directory')
dir.create(file.path(out,'_provenance'),recursive=TRUE,showWarnings=FALSE)
x<-review_read(input,cfg)
required<-c('population','comparison_id','ID','Description','database','method','direction','p.adjust','NES','Count','geneID','core_enrichment')
stopifnot(all(required%in%names(x)))
x<-x[x$population==cfg$population & x$comparison_id==cfg$comparison_id,,drop=FALSE]

if(anyDuplicated(x[,c('database','method','direction','ID')]))stop('Duplicated pathway results')
stopifnot(all(x$method%in%c('ORA','GSEA')))
cutoff<-get('fdr_threshold',.05);top<-as.integer(get('top_per_direction',6));budget<-as.integer(get('matrix_per_database',4));page_n<-as.integer(get('terms_per_page',12));overlap<-get('max_jaccard',.6)
stopifnot(cutoff>0,cutoff<=1,top>0,budget>0,page_n>0,overlap>=0,overlap<=1)
style<-get('style','dot');fmt<-get('figure_format','pdf');stopifnot(style%in%c('dot','lollipop','both'),fmt%in%c('pdf','png'))
x$plot_direction<-ifelse(x$method=='GSEA',ifelse(x$NES>=0,'Up','Down'),x$direction)
stopifnot(all(x$plot_direction%in%c('Up','Down')))
x$column<-ifelse(x$method=='GSEA','GSEA',paste('ORA',x$direction))
x$score<-pmin(6,-log10(pmax(x$p.adjust,1e-300)))
x$significant<-is.finite(x$p.adjust)&x$p.adjust<=cutoff
if(!'RichFactor'%in%names(x))x$RichFactor<-NA_real_
missing<-x$method=='ORA'&!is.finite(x$RichFactor)
if(any(missing)&&'BgRatio'%in%names(x))x$RichFactor[missing]<-x$Count[missing]/as.numeric(sub('/.*','',x$BgRatio[missing]))
if(any(!is.finite(x$RichFactor[x$method=='ORA'])))stop('ORA needs finite RichFactor or Count and BgRatio')
clean<-function(s){s<-sub('^(HALLMARK|REACTOME|KEGG)_','',s);s<-gsub('_',' ',s);s<-vapply(s,function(z)if(z==toupper(z))paste0(toupper(substr(z,1,1)),tolower(substr(z,2,nchar(z))))else z,character(1));for(a in c('TNFA','NFKB','RNA','DNA','MHC','ATP','GTP','TGF','IFN'))s<-gsub(paste0('\\b',tolower(a),'\\b'),if(a=='NFKB')'NF-kB'else if(a=='TNFA')'TNFa'else a,s,perl=TRUE,ignore.case=TRUE);vapply(s,function(z)paste(strwrap(z,43),collapse='\n'),character(1))}
x$label<-clean(x$Description)
write_tsv<-function(d,n){for(k in names(d))if(is.character(d[[k]]))d[[k]]<-gsub('[\r\n\t]',' ',d[[k]]);write.table(d,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
# Select within method/direction; overlap comparisons never mix identifier systems.
audit<-list();selected<-integer();ai<-0
for(g in split(seq_len(nrow(x)),interaction(x$database,x$method,x$plot_direction,drop=TRUE))){
 g<-g[x$significant[g]];g<-g[order(x$p.adjust[g],x$ID[g])];kept<-integer()
 for(i in g){
  field<-if(x$method[i]=='GSEA')'core_enrichment' else 'geneID'
  genes<-function(j){s<-x[[field]][j];if(is.na(s)||!nzchar(s))character() else unique(strsplit(s,'/',fixed=TRUE)[[1]])}
  gi<-genes(i);js<-vapply(kept,function(j){gj<-genes(j);u<-union(gi,gj);if(!length(u))0 else length(intersect(gi,gj))/length(u)},numeric(1))
  red<-length(js)>0&&max(js)>overlap
  reason<-if(red)'redundant' else if(length(kept)>=top)'display_budget' else 'selected'
  ai<-ai+1;audit[[ai]]<-data.frame(database=x$database[i],method=x$method[i],direction=x$plot_direction[i],ID=x$ID[i],p.adjust=x$p.adjust[i],decision=reason,redundant_with=if(red)x$ID[kept[which.max(js)]]else NA_character_)
  if(reason=='selected')kept<-c(kept,i)
 };selected<-c(selected,kept)
}
if(!length(selected))review_empty(out,cfg,'No terms pass the display FDR threshold')
d<-x[selected,,drop=FALSE];write_tsv(do.call(rbind,audit),'selection_audit.tsv');write_tsv(d,'detail_selected_terms.tsv')
if(!is.null(cfg$databases)){x<-x[x$database%in%cfg$databases,,drop=FALSE];d<-d[d$database%in%cfg$databases,,drop=FALSE]}
dbs<-intersect(c('HALLMARK','REACTOME','KEGG','GO_BP','GO_MF','GO_CC'),unique(x$database));dbs<-c(dbs,setdiff(unique(x$database),dbs))
cols<-c(Up='#B74738',Down='#326D9B');title<-get('title',paste(cfg$population,cfg$comparison_id))
th<-theme_bw(base_size=11)+theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(colour='#252525'),strip.background=element_rect(fill='#F2F2F2',colour=NA),strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=9),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(8,10,8,8))
files<-character()
save<-function(p,stem,w,h){name<-paste0(stem,'.',fmt);ggsave(file.path(out,name),p,width=w,height=h,units='in',dpi=300,bg='white',limitsize=FALSE);files<<-c(files,name)}
styles<-if(style=='both')c('dot','lollipop')else style
for(db in dbs)for(m in c('GSEA','ORA')){
 z<-d[d$database==db&d$method==m,,drop=FALSE];if(!nrow(z))next
 z<-z[order(z$plot_direction,z$p.adjust,z$ID),,drop=FALSE]
 # Budget page size per direction so ORA panels keep both directions together.
 z$page<-ave(seq_len(nrow(z)),z$plot_direction,FUN=function(v)ceiling(seq_along(v)/max(1,floor(page_n/2))))
 for(pg in sort(unique(z$page))){v<-z[z$page==pg,,drop=FALSE];v$key<-paste(v$ID,v$plot_direction,sep=':');v$key<-factor(v$key,levels=rev(v$key));labels<-setNames(v$label,as.character(v$key))
  v$effect<-if(m=='GSEA')v$NES else v$RichFactor
  v$evidence<-factor(ifelse(v$p.adjust<=.05,'FDR <= 0.05','Exploratory'),levels=c('FDR <= 0.05','Exploratory'))
  v$direction_label<-if(m=='GSEA')ifelse(v$plot_direction=='Up',get('positive_label','Numerator'),get('negative_label','Denominator'))else paste(v$plot_direction,'input genes')
  for(st in styles){
   p<-ggplot(v,aes(effect,key))
   if(st=='lollipop')p<-p+geom_segment(aes(x=0,xend=effect,yend=key),colour='grey80',linewidth=.5)
   p<-p+geom_point(aes(colour=score,shape=evidence),size=3.5,stroke=1)+scale_colour_gradient(low='#C4D6DF',high='#163E5B',limits=c(0,6),name='-log10 FDR\n(capped at 6)')+scale_shape_manual(values=c(16,1),drop=TRUE)+scale_y_discrete(labels=labels)+th+labs(title=paste(db,m),subtitle=title,x=if(m=='GSEA')'Normalized enrichment score (NES)'else 'Rich factor (overlap / pathway genes in background)',y=NULL,caption=paste0('Display FDR <= ',cutoff,'; up to ',top,' nonredundant terms per direction.'))
   if(m=='GSEA')p<-p+scale_x_continuous(limits=c(-1,1)*max(abs(z$NES))*1.08)+geom_vline(xintercept=0,colour='grey65',linewidth=.4)+guides(shape=guide_legend(title=NULL)) else p<-p+scale_x_continuous(limits=c(0,max(z$RichFactor)*1.08))+facet_wrap(~direction_label,ncol=min(2,length(unique(v$direction_label))),scales='free_y')+guides(shape=guide_legend(title=NULL))
   if(length(unique(v$evidence))==1)p<-p+guides(shape='none')
   if(m=='GSEA')p<-p+labs(x=paste0('NES: ',get('negative_label','Denominator'),' <-- 0 --> ',get('positive_label','Numerator')))
   n<-if(m=='ORA')max(table(v$direction_label))else nrow(v)
   save(p,paste(tolower(db),tolower(m),st,paste0('page',pg),sep='_'),if(m=='ORA'&&length(unique(v$direction_label))>1)11 else 8.8,min(8.5,max(3,1.8+.29*n)))
  }
 }
}
# Round-robin representatives across methods and directions, exact-ID matches only.
mi<-integer()
for(db in dbs){z<-d[d$database==db,,drop=FALSE];qs<-list();for(m in c('GSEA','ORA'))for(dir in c('Up','Down')){a<-z[z$method==m&z$plot_direction==dir,,drop=FALSE];qs[[paste(m,dir)]]<-a$ID[order(a$p.adjust,a$ID)]};ids<-character();depth<-1
 while(length(ids)<budget&&any(lengths(qs)>=depth)){
  for(q in qs)if(length(q)>=depth&&!(q[depth]%in%ids)&&length(ids)<budget)ids<-c(ids,q[depth]);depth<-depth+1
 }
 for(id in ids)mi<-c(mi,which(x$database==db&x$ID==id))
}
mat<-x[unique(mi),,drop=FALSE];write_tsv(mat,'matrix_complete_cells.tsv')
# one pathway per row with complete cells, grey nonsignificant, cross absent.
ids<-unique(paste(mat$database,mat$ID,sep='::'));grid<-expand.grid(key=ids,column=c('GSEA','ORA Up','ORA Down'),stringsAsFactors=FALSE)
lookup<-match(paste(grid$key,grid$column),paste(paste(mat$database,mat$ID,sep='::'),mat$column))
grid$database<-sub('::.*','',grid$key);grid$database<-factor(grid$database,levels=dbs)
grid$column<-factor(grid$column,levels=c('GSEA','ORA Up','ORA Down'));grid$key<-factor(grid$key,levels=rev(ids))
grid$available<-!is.na(lookup);grid$fdr<-mat$p.adjust[lookup];grid$fill<-ifelse(grid$available&is.finite(grid$fdr)&grid$fdr<=cutoff,ifelse(mat$plot_direction[lookup]=='Up',1,-1)*mat$score[lookup],NA_real_)
grid$text<-ifelse(!grid$available,'x',ifelse(mat$method[lookup]=='GSEA',sprintf('%.2f',mat$NES[lookup]),as.character(mat$Count[lookup])))
grid$text_color<-ifelse(is.na(grid$fill),'#727272',ifelse(abs(grid$fill)>3,'white','#222222'))
labels<-setNames(mat$label,paste(mat$database,mat$ID,sep='::'));write_tsv(grid,'matrix_plot_values.tsv')
# Bound matrix to 24 rows per page.
for(pg in seq_len(ceiling(length(ids)/24))){rows<-ids[seq.int((pg-1)*24+1,min(pg*24,length(ids)))];v<-grid[as.character(grid$key)%in%rows,,drop=FALSE]
 p<-ggplot(v,aes(column,key))+geom_tile(aes(fill=fill),colour='white',linewidth=.6,width=.96,height=.96)+geom_text(aes(label=text,colour=text_color),size=3.2)+scale_colour_identity()+scale_fill_gradient2(low=cols['Down'],mid='#FAFAFA',high=cols['Up'],limits=c(-6,6),na.value='#EEEEEE',name='Signed -log10 FDR\n(capped at +/-6)')+scale_y_discrete(labels=labels)+facet_grid(database~.,scales='free_y',space='free_y')+th+theme(panel.grid=element_blank(),axis.text.x=element_text(face='bold'),strip.text.y=element_text(angle=0),legend.key.width=grid::unit(1.4,'cm'))+labs(title='Pathway evidence across methods',subtitle=title,x=NULL,y=NULL,caption=paste0('Red: ',get('positive_label','numerator'),'/Up input; blue: ',get('negative_label','denominator'),'/Down input.\nCell text: GSEA NES; ORA overlap count. Grey: FDR > ',cutoff,' or unavailable FDR; x: no result row.\nExact pathway matches within each database; not cross-database equivalence.'))
 save(p,paste0('evidence_matrix_page',pg),9.5,max(4,2.5+.27*length(rows)))
}
files<-review_merge_details(out,files,fmt)
write_json(list(renderer_version='compact_v1',renderer_md5=unname(tools::md5sum(renderer_path)),input=input,input_md5=unname(tools::md5sum(input)),config=cfg,created=as.character(Sys.time()),n_input_selected=nrow(x),n_detail_terms=nrow(d),n_matrix_rows=length(ids),figures=files),file.path(out,'_provenance/manifest.json'),pretty=TRUE,auto_unbox=TRUE)
file.copy(args[1],file.path(out,'_provenance/config.json'));capture.output(sessionInfo(),file=file.path(out,'_provenance/session_info.txt'))
writeLines(c('# Compact enrichment review','','Status: completed','','Statistics were not recomputed. Matrix colours encode evidence, not effect size.','',paste0('- [',files,'](',files,')'),'','Selection and values: [detail terms](detail_selected_terms.tsv), [selection audit](selection_audit.tsv), [matrix values](matrix_plot_values.tsv).'),file.path(out,'RESULTS.md'))
cat('COMPACT_ENRICHMENT_COMPLETE:',length(files),'figures\n')

review_provenance(out,cfg,renderer_path)
