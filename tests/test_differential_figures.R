source("toolkit/R/runtime.R")
source("toolkit/R/figure_style.R")
source("toolkit/R/differential_utils.R")

out <- file.path("test-output", "differential-figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
config <- list(output = list(figure_format = "png"),
               figure_style = list(colors = list(condition = list(control = "#4477AA", pah = "#D55E00"))),
               plots = list(top_genes = 4))
comparison <- list(id = "pah_vs_control", numerator = "pah", denominator = "control")
result <- data.frame(
  gene = paste0("Gene", seq_len(8)), log2FoldChange = c(-2, -1, -.5, -.1, .1, .5, 1, 2),
  padj = c(.001, .01, .2, .8, .8, .2, .01, .001), baseMean = seq(10, 80, 10),
  significance = c("Down", "Down", "NS", "NS", "NS", "NS", "Up", "Up"), stringsAsFactors = FALSE
)
plot_de_results(result, out, "endothelial", comparison, config)
norm <- matrix(seq_len(32), nrow = 8, dimnames = list(result$gene, c("S1", "S2", "S3", "S4")))
coldata <- data.frame(condition = c("control", "control", "pah", "pah"), row.names = colnames(norm))
plot_pseudobulk(norm, coldata, "condition", result, out, config)
all_results <- cbind(result, population = "endothelial", comparison_id = comparison$id)
plot_batch_summary(all_results, data.frame(), out, config)
expected <- c("volcano.png", "MA_plot.png", "pseudobulk_PCA.png", "top_DE_heatmap.png", "DEG_count_summary.png")
stopifnot(all(file.exists(file.path(out, expected))))
status <- utils::read.delim(file.path(out, "_provenance/figure_status.tsv"), check.names = FALSE)
stopifnot(setequal(status$family, sub("[.]png$", "", expected)), all(status$style_version == "paper_v1"))
cat("PASS: differential paper figures and provenance\n")
