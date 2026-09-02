args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript differential_analysis.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
source(file.path(dirname(normalizePath(script_file)), "differential_utils.R"))

config <- read_skill_config(args[[1]])
stage <- cfg_get(config, "analysis.stage", "differential")
if (!stage %in% c("differential", "differential_and_enrichment", "enrichment_only")) stop("Unsupported analysis.stage: ", stage)
if (stage == "enrichment_only") {
  run_enrichment_only_workflow(config)
  quit(save = "no", status = 0L)
}
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required")

obj <- load_scrna_object(cfg_get(config, "input.object", required = TRUE), "auto")
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
condition_col <- cfg_get(config, "metadata.condition", required = TRUE)
population_col <- cfg_get(config, "metadata.population", cfg_get(config, "metadata.cell_type"))
covariates <- as_chr(cfg_get(config, "metadata.covariates", list()))
assert_metadata(obj, unique(c(sample_col, condition_col, population_col, covariates)))
meta <- obj[[]]
validate_sample_mapping(meta, sample_col, condition_col, covariates)

assay <- cfg_get(config, "analysis.assay", "RNA")
if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
obj <- join_assay_layers(obj, assay)
method <- cfg_get(config, "analysis.method", "pseudobulk_deseq2")

comparisons <- normalize_comparisons(config)
populations <- select_populations(meta, population_col, cfg_get(config, "population", list()))
out <- prepare_output(config)
dir.create(file.path(out, "comparisons"), showWarnings = FALSE, recursive = TRUE)
thresholds <- list(
  padj = as.numeric(cfg_get(config, "analysis.padj_threshold", 0.05)),
  lfc = as.numeric(cfg_get(config, "analysis.lfc_threshold", 0.25)),
  min_cells = as.integer(cfg_get(config, "analysis.min_cells_per_sample_population", 10)),
  min_samples = as.integer(cfg_get(config, "analysis.min_samples_per_group", 2)),
  min_total_count = as.integer(cfg_get(config, "analysis.min_total_count", 10)),
  min_count_per_sample = as.integer(cfg_get(config, "analysis.min_count_per_sample", 10)),
  min_samples_expressed = as.integer(cfg_get(config, "analysis.min_samples_expressed", cfg_get(config, "analysis.min_samples_per_group", 2)))
)

audit <- make_design_audit(meta, sample_col, condition_col, population_col, populations, comparisons)
write_tsv(audit, file.path(out, "design_audit.tsv"))
all_results <- list(); status_rows <- list(); enrichment_rows <- list(); task_index <- 0L

for (population in populations) for (comparison in comparisons) {
  task_index <- task_index + 1L
  task_id <- safe_name(paste(population, comparison$id, sep = "__"))
  task_dir <- file.path(out, "comparisons", task_id)
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  task <- tryCatch({
    {
      cells <- population_cells(meta, population_col, population)
      task_meta <- meta[cells, , drop = FALSE]
      task_meta <- task_meta[task_meta[[condition_col]] %in% c(comparison$numerator, comparison$denominator), , drop = FALSE]
      cells <- rownames(task_meta)
      if (!length(cells)) stop("No cells remain for this population and comparison")
      cell_audit <- sample_population_audit(task_meta, sample_col, condition_col, thresholds$min_cells)
      write_tsv(cell_audit, file.path(task_dir, "sample_cell_counts.tsv"))
      keep_samples <- cell_audit$sample[cell_audit$n_cells >= thresholds$min_cells]
      task_meta <- task_meta[task_meta[[sample_col]] %in% keep_samples, , drop = FALSE]
      counts_by_group <- table(unique(task_meta[c(sample_col, condition_col)])[[condition_col]])
      missing_groups <- setdiff(c(comparison$numerator, comparison$denominator), names(counts_by_group))
      if (length(missing_groups) || any(counts_by_group[c(comparison$numerator, comparison$denominator)] < thresholds$min_samples)) {
        stop("Insufficient independent samples after minimum-cell filtering")
      }
      sub <- subset(obj, cells = rownames(task_meta))
      if (method == "pseudobulk_deseq2") {
        result <- run_pseudobulk(sub, task_meta, assay, sample_col, condition_col, covariates, comparison, thresholds, config, task_dir)
      } else if (method %in% c("seurat_wilcox", "seurat_mast", "seurat_lr")) {
        result <- run_cell_level(sub, assay, condition_col, comparison, method, config)
      } else stop("Unsupported analysis.method: ", method)
      result <- annotate_de_result(result, population, comparison, method, thresholds, task_meta, sample_col, condition_col)
      write_de_tables(result, task_dir)
      plot_de_results(result, task_dir, population, comparison, config)
    }
    enr <- NULL
    if (isTRUE(cfg_get(config, "enrichment.enabled", FALSE)) || stage %in% c("differential_and_enrichment", "enrichment_only")) {
      enr <- tryCatch(run_enrichment(result, task_dir, population, comparison, config), error = function(e) {
        writeLines(conditionMessage(e), file.path(task_dir, "ENRICHMENT_ERROR.txt"))
        message("Enrichment for ", task_id, " failed without discarding DE results: ", conditionMessage(e))
        data.frame()
      })
    }
    all_results[[task_id]] <- result
    if (!is.null(enr) && nrow(enr)) enrichment_rows[[task_id]] <- enr
    enrichment_status <- if (is.null(enr)) "not_requested" else summarize_enrichment_status(task_dir)
    data.frame(task_id, population, comparison_id = comparison$id, status = "completed", enrichment_status, message = "", n_genes = nrow(result), stringsAsFactors = FALSE)
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(task_dir, "ERROR.txt"))
    message("Task ", task_id, " failed: ", conditionMessage(e))
    data.frame(task_id, population, comparison_id = comparison$id, status = classify_failure(conditionMessage(e)), enrichment_status = "not_run", message = conditionMessage(e), n_genes = 0L, stringsAsFactors = FALSE)
  })
  status_rows[[task_index]] <- task
}

status <- do.call(rbind, status_rows)
write_tsv(status, file.path(out, "task_status.tsv"))
artifacts <- c(file.path(out, "design_audit.tsv"), file.path(out, "task_status.tsv"))
if (length(all_results)) {
  combined <- do.call(rbind, all_results)
  write_tsv(combined, file.path(out, "all_comparisons.tsv"))
  write_tsv(combined[combined$significance %in% c("Up", "Down"), , drop = FALSE], file.path(out, "significant_all_comparisons.tsv"))
  plot_batch_summary(combined, status, out)
  artifacts <- c(artifacts, file.path(out, "all_comparisons.tsv"), file.path(out, "significant_all_comparisons.tsv"))
}
if (length(enrichment_rows)) {
  enrichment <- rbind_fill(enrichment_rows)
  write_tsv(enrichment, file.path(out, "enrichment_all_comparisons.tsv"))
  artifacts <- c(artifacts, file.path(out, "enrichment_all_comparisons.tsv"))
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
artifacts <- c(artifacts, file.path(out, "sessionInfo.txt"))
active_skill <- Sys.getenv("SCRNA_ACTIVE_SKILL", unset = "11-scrna-run-differential-analysis")
write_run_manifest(config, active_skill, out, artifacts,
                   c(paste0("method=", method), paste0("stage=", stage), "Positive log2 fold change means numerator > denominator", "Input object was not rewritten"))
if (!any(status$status == "completed")) stop("No differential-analysis task completed; inspect task_status.tsv")
