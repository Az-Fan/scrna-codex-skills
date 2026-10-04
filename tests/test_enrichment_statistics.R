source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")
root <- tempfile("scrna-enrichment-statistics-"); dir.create(root)
# Exercise real identifier mapping and GO computation, including genes filtered
# from FDR eligibility that must remain in the GSEA ranking.
genes <- c("Kdr", "Pecam1", "Cdh5", "Ptprc", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Col1a1", "Col3a1",
  "Dcn", "Lyz2", "Ms4a1", "Nkg7", "Ifng", "Tnf", "Il6", "Vegfa", "Acta2", "Rgs5")
x <- data.frame(gene = genes, log2FoldChange = seq(2, -2, length.out = 20),
  stat = seq(5, -5, length.out = 20), pvalue = .001, padj = c(rep(.01, 15), rep(NA_real_, 5)),
  tested = TRUE, multiple_testing_eligible = c(rep(TRUE, 15), rep(FALSE, 5)),
  filter_reason = c(rep("", 15), rep("independent_filtering", 5)))
config <- list(enrichment = list(species = "mouse", gene_id_type = "SYMBOL", databases = list("GO_BP"),
  min_input_genes = 1, min_gene_set_size = 1, max_gene_set_size = 500))
x <- standardize_de_table(x, config)
plot_enrichment_summary <- function(...) invisible(NULL)
set.seed(1)
invisible(run_enrichment(x, root, "all", list(id = "case_vs_control"), config))
directory <- file.path(root, "enrichment")
background <- utils::read.delim(file.path(directory, "ora_universe_genes.tsv"))
rank <- utils::read.delim(file.path(directory, "gsea_ranked_genes.tsv"))
excluded <- background[!background$selected_for_ora & !is.na(background$mapped_id), ]
stopifnot(nrow(excluded) == 5L, all(excluded$mapped_id %in% rank$mapped_id), nrow(rank) == 20L)
for (direction in c("up", "down")) {
  foreground <- utils::read.delim(file.path(directory, paste0("ora_go_bp_", direction, "_foreground.tsv")))
  effective <- utils::read.delim(file.path(directory, paste0("ora_go_bp_", direction, "_effective_universe.tsv")))
  stopifnot(!any(foreground$mapped_id %in% excluded$mapped_id), !any(effective$mapped_id %in% excluded$mapped_id),
    all(foreground$mapped_id %in% background$mapped_id[background$selected_for_ora]))
}
status <- utils::read.delim(file.path(directory, "enrichment_status.tsv"))
stopifnot(nrow(status) == 3L, all(status$status %in% c("completed", "empty")),
  status$n_input[status$method == "GSEA"] == 20L)
unlink(root, recursive = TRUE)
cat("PASS real GO ORA eligible background and foregrounds; GSEA retains independently filtered ranked genes\n")
