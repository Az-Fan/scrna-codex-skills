args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript differential_analysis.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
source(file.path(dirname(normalizePath(script_file)), "figure_style.R"))
source(file.path(dirname(normalizePath(script_file)), "differential_utils.R"))

config <- read_skill_config(args[[1]])
stage <- cfg_get(config, "analysis.stage", "differential")
if (!stage %in% c("differential", "differential_and_enrichment", "enrichment_only")) stop("Unsupported analysis.stage: ", stage)
if (stage == "enrichment_only") {
  run_enrichment_only_workflow(config)
  quit(save = "no", status = 0L)
}
thresholds <- validate_de_thresholds(config)
if (cfg_get(config, "analysis.method", "pseudobulk_deseq2") == "pseudobulk_deseq2") {
  validate_de_design(config, cfg_get(config, "metadata.condition", required = TRUE), as_chr(cfg_get(config, "metadata.covariates", list())))
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
counts_source_audit <- if (method == "pseudobulk_deseq2") assert_pseudobulk_counts_source(obj, assay, config) else NULL

comparisons <- normalize_comparisons(config)
populations <- select_populations(meta, population_col, cfg_get(config, "population", list()))
task_ids <- unlist(lapply(populations, function(population) vapply(comparisons, function(comparison) safe_name(paste(population, comparison$id, sep = "__")), character(1))), use.names = FALSE)
if (anyDuplicated(task_ids)) stop("Differential task IDs collide after sanitization; use distinct population/comparison identifiers")
out <- prepare_output(config)
if (!is.null(counts_source_audit)) write_tsv(counts_source_audit, file.path(out, "counts_source_audit.tsv"))
dir.create(file.path(out, "comparisons"), showWarnings = FALSE, recursive = TRUE)
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
      replicate_counts <- counts_by_group[match(c(comparison$numerator, comparison$denominator), names(counts_by_group))]
      replicate_counts[is.na(replicate_counts)] <- 0L
      low_replication <- method == "pseudobulk_deseq2" && any(replicate_counts < 3L)
      write_tsv(data.frame(condition = c(comparison$numerator, comparison$denominator), n_samples = as.integer(replicate_counts),
        low_replication_warning = low_replication,
        interpretation = if (low_replication) "exploratory_low_confidence" else if (method == "pseudobulk_deseq2") "replicated_sample_level" else "cell_level_exploratory"), file.path(task_dir, "replication_audit.tsv"))
      missing_groups <- setdiff(c(comparison$numerator, comparison$denominator), names(counts_by_group))
      if (length(missing_groups) || any(counts_by_group[c(comparison$numerator, comparison$denominator)] < thresholds$min_samples)) {
        stop("Insufficient independent samples after minimum-cell filtering")
      }
      sub <- subset(obj, cells = rownames(task_meta))
      if (low_replication) warning(task_id, ": fewer than 3 biological replicates in at least one group; exploratory / low-confidence inference")
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
    data.frame(task_id, population, comparison_id = comparison$id, status = "completed", enrichment_status, message = "", n_genes = nrow(result), low_replication_warning = low_replication, stringsAsFactors = FALSE)
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(task_dir, "ERROR.txt"))
    message("Task ", task_id, " failed: ", conditionMessage(e))
    data.frame(task_id, population, comparison_id = comparison$id, status = classify_failure(conditionMessage(e)), enrichment_status = "not_run", message = conditionMessage(e), n_genes = 0L, low_replication_warning = NA, stringsAsFactors = FALSE)
  })
  status_rows[[task_index]] <- task
}

status <- do.call(rbind, status_rows)
write_tsv(status, file.path(out, "task_status.tsv"))
artifacts <- c(file.path(out, "design_audit.tsv"), file.path(out, "task_status.tsv"))
if (!is.null(counts_source_audit)) artifacts <- c(artifacts, file.path(out, "counts_source_audit.tsv"))
if (length(all_results)) {
  combined <- do.call(rbind, all_results)
  write_tsv(combined, file.path(out, "all_comparisons.tsv"))
  write_tsv(combined[combined$significance %in% c("Up", "Down"), , drop = FALSE], file.path(out, "significant_all_comparisons.tsv"))
  plot_batch_summary(combined, status, out, config)
  artifacts <- c(artifacts, file.path(out, "all_comparisons.tsv"), file.path(out, "significant_all_comparisons.tsv"))
}
if (length(enrichment_rows)) {
  enrichment <- rbind_fill(enrichment_rows)
  write_tsv(enrichment, file.path(out, "enrichment_all_comparisons.tsv"))
  artifacts <- c(artifacts, file.path(out, "enrichment_all_comparisons.tsv"))
}
writeLines(capture.output(sessionInfo()), technical_path(out, "session_info.txt"))
artifacts <- c(artifacts, technical_path(out, "session_info.txt"))
comparison_artifacts <- list.files(file.path(out, "comparisons"), recursive = TRUE, full.names = TRUE)
artifacts <- c(artifacts, comparison_artifacts[file.info(comparison_artifacts)$isdir %in% FALSE])
active_skill <- Sys.getenv("SCRNA_ACTIVE_SKILL", unset = "11-scrna-run-differential-analysis")
write_run_manifest(config, active_skill, out, artifacts,
                   c(paste0("method=", method), paste0("stage=", stage), "Positive log2 fold change means numerator > denominator", "Input object was not rewritten"),
                   exit_status = if (any(status$status == "completed")) 0L else 1L)
if (!any(status$status == "completed")) stop("No differential-analysis task completed; inspect task_status.tsv")
