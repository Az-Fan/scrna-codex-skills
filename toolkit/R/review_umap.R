args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: Rscript review_umap.R CONFIG.json")
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"runtime.R"))
source(file.path(dirname(normalizePath(script)),"umap_review.R"))
config <- read_skill_config(args[[1]])
invisible(cfg_get(config,"project.id",required=TRUE))
obj <- load_scrna_object(cfg_get(config,"input.object",required=TRUE),"auto")
out <- prepare_output(config)
files <- review_cluster_umap(obj,config,out)
session <- technical_path(out,"session_info.txt")
writeLines(capture.output(sessionInfo()),session)
write_run_manifest(config,"06-scrna-preprocess-and-cluster",out,c(files,session),
 notes=c("UMAP review only: existing embedding reused; input object unchanged",
         "Sample representation is descriptive; UMAP overlap does not establish absence of batch effects"))
links <- basename(files)
writeLines(c("# UMAP review results","","Status: completed","",
  "Existing embedding reused. Grey: all cells; blue: titled sample. UMAP overlap alone does not establish absence of batch effects.","",
  paste0("- [",links,"](",links,")")),file.path(out,"RESULTS.md"))
