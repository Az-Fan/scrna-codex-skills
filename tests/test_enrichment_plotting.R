source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")

fixture <- data.frame(
  Description = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_APICAL_JUNCTION",
    "GO: response_to_stress", paste("GO:", paste(rep("endothelial_cell_migration", 5), collapse = "_"))
  ),
  database = c(rep("HALLMARK", 4), rep("GO_BP", 2)),
  method = c(rep("GSEA", 4), rep("ORA", 2)),
  direction = c(rep("ranked", 4), "Up", "Down"),
  enrichment_direction = c("Up", "Up", "Down", "Down", NA, NA),
  NES = c(2.1, 1.8, -2.0, -1.7, NA, NA),
  p.adjust = c(0.01, 0.02, 0.015, 0.20, 0.01, 0.02),
  setSize = c(40, 20, 35, 25, NA, NA),
  RichFactor = c(NA, NA, NA, NA, 0.4, 0.3),
  Count = c(NA, NA, NA, NA, 12, 9),
  population = "vascular_ec",
  comparison_id = "pah_vs_control",
  core_enrichment = c("1/2/3", "1/2/3/4", "8/9/10", "30/31/32", NA, NA),
  geneID = c(NA, NA, NA, NA, "12/13/14", "20/21/22"),
  stringsAsFactors = FALSE
)

selected <- select_enrichment_plot_terms(fixture, top_n = 1, ora_fdr = 0.05, gsea_fdr = 0.25)
gsea_selected <- selected[selected$method == "GSEA", , drop = FALSE]
stopifnot(nrow(gsea_selected) == 2L, setequal(gsea_selected$plot_direction, c("Up", "Down")))
stopifnot(unname(clean_enrichment_label("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", 80)) ==
          "REACTIVE OXYGEN SPECIES PATHWAY")

out <- tempfile("enrichment_plot_test_")
dir.create(out)
config <- list(enrichment = list(plot_top_terms = 2, plot_label_width = 32,
                               plot_terms_per_page = 20, plot_ora_fdr_threshold = 0.05,
                               plot_gsea_fdr_threshold = 0.25, plot_max_gene_overlap = 0.6,
                               plot_population_label = "vascular endothelial cells",
                               plot_comparison_label = "PAH vs control",
                               plot_positive_label = "PAH-enriched",
                               plot_negative_label = "Control-enriched"))
plot_warnings <- character()
withCallingHandlers(plot_enrichment_summary(fixture, out, config), warning = function(warning) {
  plot_warnings <<- c(plot_warnings, conditionMessage(warning))
  invokeRestart("muffleWarning")
})
expected <- c("enrichment_dotplot_overview.pdf", "enrichment_ora_overview.pdf",
              "enrichment_dotplot_go_bp_ora.pdf", "gsea_nes_hallmark.pdf",
              "enrichment_plot_terms_summary.tsv")
stopifnot(all(file.exists(file.path(out, expected))))
stopifnot(!any(grepl("Removed .* rows", plot_warnings)))
plot_terms <- read.delim(file.path(out, "enrichment_plot_terms_summary.tsv"), check.names = FALSE)
stopifnot(all(c("plot_rank", "plot_direction_label", "evidence_class") %in% names(plot_terms)),
          nrow(plot_terms) == 5L,
          !any(grepl("\n", plot_terms$label, fixed = TRUE)),
          !any(grepl("\r", plot_terms$label, fixed = TRUE)),
          !any(grepl("\t", plot_terms$label, fixed = TRUE)),
          !"HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" %in% plot_terms$Description)
unlink(out, recursive = TRUE)

input <- tempfile("differential_", fileext = ".tsv")
writeLines("gene\tlog2FoldChange\nA\t1", input)
config_path <- tempfile("enrichment_config_", fileext = ".json")
jsonlite::write_json(list(project = list(id = "test"), input = list(differential_table = input)),
                     config_path, auto_unbox = TRUE)
parsed <- read_skill_config(config_path)
stopifnot(identical(attr(parsed, "input_record")$path, input),
          is.finite(attr(parsed, "input_record")$bytes),
          nzchar(attr(parsed, "input_record")$sha256))
unlink(c(input, config_path))
cat("PASS: enrichment plotting selection, labels, and PDF outputs\n")
