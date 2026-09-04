source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")

fixture <- data.frame(
  Description = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_APICAL_JUNCTION",
    "GO: response_to_stress", "GO: endothelial_cell_migration"
  ),
  database = c(rep("HALLMARK", 4), rep("GO_BP", 2)),
  method = c(rep("GSEA", 4), rep("ORA", 2)),
  direction = c(rep("ranked", 4), "Up", "Down"),
  enrichment_direction = c("Up", "Up", "Down", "Down", NA, NA),
  NES = c(2.1, 1.8, -2.0, -1.7, NA, NA),
  p.adjust = c(0.01, 0.02, 0.015, 0.03, 0.01, 0.02),
  setSize = c(40, 20, 35, 25, NA, NA),
  RichFactor = c(NA, NA, NA, NA, 0.4, 0.3),
  Count = c(NA, NA, NA, NA, 12, 9),
  population = "vascular_ec",
  comparison_id = "pah_vs_control",
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
                               plot_gsea_fdr_threshold = 0.25))
plot_enrichment_summary(fixture, out, config)
expected <- c("enrichment_dotplot_overview.pdf", "enrichment_ora_overview.pdf",
              "enrichment_dotplot_go_bp_ora.pdf", "gsea_nes_hallmark.pdf")
stopifnot(all(file.exists(file.path(out, expected))))
unlink(out, recursive = TRUE)
cat("PASS: enrichment plotting selection, labels, and PDF outputs\n")
