#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript cell_abundance.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))

as_chr <- function(x) unname(unlist(x %||% list(), use.names = FALSE))
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
rbind_fill <- function(frames) {
  frames <- Filter(function(x) !is.null(x) && is.data.frame(x), frames)
  if (!length(frames)) return(data.frame())
  columns <- unique(unlist(lapply(frames, colnames), use.names = FALSE))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(columns, colnames(frame))
    for (column in missing) frame[[column]] <- NA
    frame[columns]
  })
  result <- do.call(rbind, frames)
  rownames(result) <- NULL
  result
}
require_pkg <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing required package: ", package)
}
classify_failure <- function(message) {
  text <- tolower(message)
  if (grepl("missing required package|executable not found|requires.*package", text)) return("missing_dependency")
  if (grepl("insufficient|at least .* samples|replicate", text)) return("skipped_low_replicates")
  if (grepl("rank deficient|confound|unique.*condition", text)) return("invalid_design")
  "failed"
}
bool_value <- function(x, default = FALSE) if (is.null(x)) default else isTRUE(x)

config <- read_skill_config(args[[1]])
set.seed(as.integer(cfg_get(config, "analysis.random_seed", 1L)))
attr(config, "resolved_random_seed") <- as.integer(cfg_get(config, "analysis.random_seed", 1L))
out <- prepare_output(config)
dir.create(file.path(out, "comparisons"), recursive = TRUE, showWarnings = FALSE)

sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
condition_col <- cfg_get(config, "metadata.condition", required = TRUE)
cell_type_col <- cfg_get(config, "metadata.cell_type", required = TRUE)
covariates <- as_chr(cfg_get(config, "metadata.covariates", list()))
methods <- unique(tolower(as_chr(cfg_get(config, "analysis.methods", required = TRUE))))
fdr <- as.numeric(cfg_get(config, "analysis.fdr", 0.05))
min_samples <- as.integer(cfg_get(config, "analysis.min_samples_per_group", 2L))
min_cells <- as.integer(cfg_get(config, "analysis.min_cells_per_sample", 20L))
denominator_mode <- cfg_get(config, "analysis.denominator.mode", required = TRUE)
denominator_description <- cfg_get(config, "analysis.denominator.description", required = TRUE)
denominator_include <- as_chr(cfg_get(config, "analysis.denominator.include", list()))

if (!is.finite(fdr) || fdr <= 0 || fdr >= 1) stop("analysis.fdr must be between 0 and 1")
if (min_samples < 2L) stop("analysis.min_samples_per_group must be at least 2")
if (min_cells < 1L) stop("analysis.min_cells_per_sample must be positive")

read_counts_table <- function(path) {
  table <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  count_col <- cfg_get(config, "input.count_column", "n_cells")
  required <- unique(c(sample_col, condition_col, cell_type_col, covariates, count_col))
  missing <- setdiff(required, colnames(table))
  if (length(missing)) stop("Counts table is missing columns: ", paste(missing, collapse = ", "))
  if (anyNA(table[required])) stop("Counts table contains missing required values")
  values <- suppressWarnings(as.numeric(table[[count_col]]))
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
    stop("Counts table values must be non-negative integers")
  }
  table$n_cells_internal <- as.integer(round(values))
  table
}

object_path <- cfg_get(config, "input.object")
counts_path <- cfg_get(config, "input.counts_table")
object <- NULL
cell_meta <- NULL
if (!is.null(object_path)) {
  require_pkg("Seurat")
  object <- load_scrna_object(object_path, "auto")
  assert_metadata(object, unique(c(sample_col, condition_col, cell_type_col, covariates)))
  cell_meta <- object[[]]
  required <- unique(c(sample_col, condition_col, cell_type_col, covariates))
  if (anyNA(cell_meta[required]) || any(vapply(cell_meta[required], function(x) any(trimws(as.character(x)) == ""), logical(1)))) {
    stop("Input object contains missing or blank required metadata")
  }
  raw_table <- cell_meta
  raw_table$n_cells_internal <- 1L
} else {
  raw_table <- read_counts_table(counts_path)
}

validate_sample_metadata <- function(table) {
  columns <- unique(c(condition_col, covariates))
  samples <- unique(as.character(table[[sample_col]]))
  rows <- lapply(samples, function(sample_id) {
    sub <- table[as.character(table[[sample_col]]) == sample_id, , drop = FALSE]
    values <- lapply(columns, function(column) unique(as.character(sub[[column]])))
    names(values) <- columns
    bad <- names(values)[vapply(values, length, integer(1)) != 1L]
    if (length(bad)) stop("Sample ", sample_id, " maps to multiple values for: ", paste(bad, collapse = ", "))
    data.frame(sample = sample_id, as.data.frame(values, stringsAsFactors = FALSE), check.names = FALSE)
  })
  result <- do.call(rbind, rows)
  colnames(result)[1] <- sample_col
  rownames(result) <- as.character(result[[sample_col]])
  result
}

sample_meta <- validate_sample_metadata(raw_table)
input_totals <- stats::aggregate(raw_table$n_cells_internal, list(sample = as.character(raw_table[[sample_col]])), sum)
colnames(input_totals) <- c(sample_col, "n_input_cells")

if (denominator_mode == "selected_cell_types") {
  missing_labels <- setdiff(denominator_include, unique(as.character(raw_table[[cell_type_col]])))
  if (length(missing_labels)) stop("Requested denominator labels are absent: ", paste(missing_labels, collapse = ", "))
  analysis_table <- raw_table[as.character(raw_table[[cell_type_col]]) %in% denominator_include, , drop = FALSE]
  if (!is.null(cell_meta)) cell_meta <- cell_meta[rownames(cell_meta) %in% rownames(analysis_table), , drop = FALSE]
} else {
  analysis_table <- raw_table
}
if (!nrow(analysis_table)) stop("No cells remain in the declared denominator")

types <- sort(unique(as.character(analysis_table[[cell_type_col]])))
samples <- as.character(sample_meta[[sample_col]])
long_counts <- stats::aggregate(
  analysis_table$n_cells_internal,
  list(sample = as.character(analysis_table[[sample_col]]), cell_type = as.character(analysis_table[[cell_type_col]])),
  sum
)
colnames(long_counts)[3] <- "n_cells"
complete_grid <- expand.grid(sample = samples, cell_type = types, stringsAsFactors = FALSE)
long_counts <- merge(complete_grid, long_counts, by = c("sample", "cell_type"), all.x = TRUE, sort = FALSE)
long_counts$n_cells[is.na(long_counts$n_cells)] <- 0L
long_counts$n_cells <- as.integer(long_counts$n_cells)
long_counts <- merge(long_counts, sample_meta, by.x = "sample", by.y = sample_col, all.x = TRUE, sort = FALSE)
colnames(long_counts)[colnames(long_counts) == "sample"] <- sample_col

count_matrix <- matrix(0L, nrow = length(types), ncol = length(samples), dimnames = list(types, samples))
for (i in seq_len(nrow(long_counts))) count_matrix[long_counts$cell_type[i], long_counts[[sample_col]][i]] <- long_counts$n_cells[i]
denominator_totals <- colSums(count_matrix)
if (any(denominator_totals == 0)) stop("At least one sample has zero cells in the declared denominator: ", paste(names(denominator_totals)[denominator_totals == 0], collapse = ", "))
proportion_matrix <- sweep(count_matrix, 2, denominator_totals, "/")

counts_output <- data.frame(cell_type = rownames(count_matrix), count_matrix, check.names = FALSE)
props_output <- data.frame(cell_type = rownames(proportion_matrix), proportion_matrix, check.names = FALSE)
write_tsv(counts_output, file.path(out, "sample_cell_counts.tsv"))
write_tsv(props_output, file.path(out, "sample_cell_proportions.tsv"))

denom_totals_df <- data.frame(sample = names(denominator_totals), n_denominator_cells = as.integer(denominator_totals), stringsAsFactors = FALSE)
colnames(denom_totals_df)[1] <- sample_col
design_audit <- merge(sample_meta, input_totals, by = sample_col, all.x = TRUE, sort = FALSE)
design_audit <- merge(design_audit, denom_totals_df, by = sample_col, all.x = TRUE, sort = FALSE)
design_audit$denominator_fraction <- design_audit$n_denominator_cells / design_audit$n_input_cells
design_audit$passes_min_cells <- design_audit$n_denominator_cells >= min_cells
design_audit$denominator_mode <- denominator_mode
design_audit$denominator_description <- denominator_description
write_tsv(design_audit, file.path(out, "design_audit.tsv"))

eligibility <- data.frame(
  cell_type = types,
  total_cells = rowSums(count_matrix),
  samples_present = rowSums(count_matrix > 0),
  samples_at_min_cells = rowSums(count_matrix >= min_cells),
  stringsAsFactors = FALSE
)
write_tsv(eligibility, file.path(out, "cell_type_eligibility.tsv"))

sample_meta_for <- function(comparison) {
  md <- sample_meta[as.character(sample_meta[[condition_col]]) %in% c(comparison$denominator, comparison$numerator), , drop = FALSE]
  md$contrast_group <- factor(
    ifelse(as.character(md[[condition_col]]) == comparison$numerator, "numerator", "denominator"),
    levels = c("denominator", "numerator")
  )
  for (column in covariates) if (is.character(md[[column]])) md[[column]] <- factor(md[[column]])
  counts_by_group <- table(md$contrast_group)
  if (any(counts_by_group < min_samples)) stop("Insufficient independent samples: ", paste(names(counts_by_group), counts_by_group, collapse = "; "))
  terms <- c(covariates, "contrast_group")
  design <- stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), md)
  if (qr(design)$rank < ncol(design)) stop("Rank deficient or confounded sample-level design")
  rownames(md) <- as.character(md[[sample_col]])
  md
}

standard_result <- function(method, comparison, feature_id, effect, effect_scale, p_value = NA_real_, adjusted_p_value = NA_real_,
                            ci_lower = NA_real_, ci_upper = NA_real_, posterior_probability = NA_real_, evidence_type = "",
                            credible = FALSE, significant = FALSE, reference = "", annotation = feature_id, note = "") {
  data.frame(
    method = method, comparison_id = comparison$id, feature_type = "cell_type", feature_id = as.character(feature_id),
    annotation = as.character(annotation), numerator = comparison$numerator, denominator = comparison$denominator,
    effect = as.numeric(effect), effect_scale = effect_scale, ci_lower = as.numeric(ci_lower), ci_upper = as.numeric(ci_upper),
    p_value = as.numeric(p_value), adjusted_p_value = as.numeric(adjusted_p_value), posterior_probability = as.numeric(posterior_probability),
    evidence_type = evidence_type, credible = as.logical(credible), significant = as.logical(significant), tested = TRUE,
    reference = as.character(reference), status = "completed", note = note, stringsAsFactors = FALSE
  )
}

run_propeller <- function(comparison, task_dir) {
  require_pkg("speckle"); require_pkg("limma")
  md <- sample_meta_for(comparison)
  cm <- count_matrix[, rownames(md), drop = FALSE]
  transform_name <- cfg_get(config, "method_options.propeller.transform", "logit")
  prop_list <- speckle::convertDataToList(cm, data.type = "counts", transform = transform_name)
  design_formula <- stats::as.formula(paste("~ 0 +", paste(c("contrast_group", covariates), collapse = " + ")))
  design <- stats::model.matrix(design_formula, md)
  colnames(design) <- make.names(colnames(design))
  numerator_column <- make.names("contrast_groupnumerator")
  denominator_column <- make.names("contrast_groupdenominator")
  if (!all(c(numerator_column, denominator_column) %in% colnames(design))) stop("Could not construct propeller contrast")
  contrast <- limma::makeContrasts(contrasts = paste0(numerator_column, "-", denominator_column), levels = design)
  result <- speckle::propeller.ttest(
    prop.list = prop_list, design = design, contrasts = contrast,
    robust = bool_value(cfg_get(config, "method_options.propeller.robust", TRUE), TRUE),
    trend = bool_value(cfg_get(config, "method_options.propeller.trend", FALSE), FALSE), sort = FALSE
  )
  raw <- data.frame(cell_type = rownames(result), result, check.names = FALSE, row.names = NULL)
  write_tsv(raw, file.path(task_dir, "propeller_results.tsv"))
  p <- result[[grep("^P\\.Value$|P.Value", colnames(result), value = TRUE)[1]]]
  adj <- result[[grep("^FDR$|adj", colnames(result), ignore.case = TRUE, value = TRUE)[1]]]
  mean_num <- rowMeans(proportion_matrix[, rownames(md)[md$contrast_group == "numerator"], drop = FALSE])
  mean_den <- rowMeans(proportion_matrix[, rownames(md)[md$contrast_group == "denominator"], drop = FALSE])
  effect <- log2((mean_num + 0.5 / mean(denominator_totals[rownames(md)])) / (mean_den + 0.5 / mean(denominator_totals[rownames(md)])))
  out_result <- standard_result("propeller", comparison, rownames(result), effect, "log2_mean_proportion_ratio",
    p, adj, evidence_type = "moderated_t_test_bh_fdr", significant = !is.na(adj) & adj <= fdr,
    note = paste0("Proportions transformed with ", transform_name, "; effect is descriptive log2 ratio of group mean proportions."))
  out_result$mean_proportion_numerator <- mean_num[out_result$feature_id]
  out_result$mean_proportion_denominator <- mean_den[out_result$feature_id]
  out_result
}

run_sccomp <- function(comparison, task_dir) {
  require_pkg("sccomp"); require_pkg("cmdstanr")
  env_prefix <- Sys.getenv("CONDA_PREFIX")
  if (!nzchar(env_prefix)) env_prefix <- normalizePath(file.path(R.home(), "../.."), mustWork = TRUE)
  cmdstan_candidate <- file.path(env_prefix, "bin", "cmdstan")
  current_cmdstan <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(current_cmdstan) && dir.exists(cmdstan_candidate)) cmdstanr::set_cmdstan_path(cmdstan_candidate)
  if (!nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
    stop("cmdstan is required to proceed; expected the registered environment at ", cmdstan_candidate)
  }
  md <- sample_meta_for(comparison)
  long <- long_counts[as.character(long_counts[[sample_col]]) %in% rownames(md), c(sample_col, "cell_type", "n_cells", covariates), drop = FALSE]
  long$contrast_group <- md[as.character(long[[sample_col]]), "contrast_group"]
  colnames(long)[colnames(long) == sample_col] <- "sample_internal"
  long$sample_internal <- as.character(long$sample_internal)
  long$n_cells <- as.integer(long$n_cells)
  terms <- c(covariates, "contrast_group")
  composition_formula <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  cache <- cfg_get(config, "method_options.sccomp.cache_dir", file.path(env_prefix, "var", "cache", "sccomp"))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  draws_dir <- file.path(task_dir, "draws")
  dir.create(draws_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- sccomp::sccomp_estimate(
    long, formula_composition = composition_formula, formula_variability = ~1,
    sample = "sample_internal", cell_group = "cell_type", abundance = "n_cells",
    cores = as.integer(cfg_get(config, "method_options.sccomp.cores", 1L)),
    inference_method = cfg_get(config, "method_options.sccomp.inference_method", "pathfinder"),
    percent_false_positive = 100 * fdr,
    output_directory = draws_dir, cache_stan_model = cache,
    max_sampling_iterations = as.integer(cfg_get(config, "method_options.sccomp.max_sampling_iterations", 20000L)),
    mcmc_seed = as.integer(cfg_get(config, "analysis.random_seed", 1L)), verbose = TRUE
  )
  result <- sccomp::sccomp_test(
    fit, contrasts = "contrast_groupnumerator", percent_false_positive = 100 * fdr,
    test_composition_above_logit_fold_change = as.numeric(cfg_get(config, "method_options.sccomp.minimum_logit_effect", 0.1))
  )
  raw <- as.data.frame(result)
  write_tsv(raw, file.path(task_dir, "sccomp_results.tsv"))
  if (bool_value(cfg_get(config, "method_options.sccomp.save_fit", FALSE), FALSE)) saveRDS(fit, file.path(task_dir, "sccomp_fit.rds"))
  cell_group_column <- if ("cell_group" %in% colnames(raw)) "cell_group" else if ("cell_type" %in% colnames(raw)) "cell_type" else NA_character_
  if (is.na(cell_group_column)) stop("sccomp result lacks a cell-group identifier column")
  standard_result("sccomp", comparison, raw[[cell_group_column]], raw$c_effect, "logit_composition_effect",
    adjusted_p_value = raw$c_FDR, ci_lower = raw$c_lower, ci_upper = raw$c_upper,
    posterior_probability = 1 - raw$c_pH0, evidence_type = "posterior_fdr",
    credible = !is.na(raw$c_FDR) & raw$c_FDR <= fdr, significant = !is.na(raw$c_FDR) & raw$c_FDR <= fdr,
    note = "sccomp robust multi-beta-binomial composition model; adjusted_p_value is posterior FDR, not a BH p-value.")
}

read_similarity_matrix <- function(path) {
  value <- utils::read.delim(path, check.names = FALSE, row.names = 1)
  result <- matrix(as.numeric(as.matrix(value)), nrow = nrow(value), dimnames = dimnames(value))
  if (nrow(result) != ncol(result)) stop("DCATS similarity matrix must be square")
  if (is.null(rownames(result)) || is.null(colnames(result)) || anyDuplicated(rownames(result)) || anyDuplicated(colnames(result))) {
    stop("DCATS similarity matrix requires unique row and column labels")
  }
  if (!setequal(rownames(result), types) || !setequal(colnames(result), types)) {
    stop("DCATS similarity matrix row and column labels must exactly match analysed cell types")
  }
  if (any(!is.finite(result))) stop("DCATS similarity matrix contains non-finite values")
  result[types, types, drop = FALSE]
}

run_dcats <- function(comparison, task_dir) {
  require_pkg("DCATS")
  md <- sample_meta_for(comparison)
  cm <- t(count_matrix[, rownames(md), drop = FALSE])
  cov_design <- if (length(covariates)) stats::model.matrix(stats::as.formula(paste("~", paste(covariates, collapse = " + "))), md)[, -1, drop = FALSE] else matrix(numeric(), nrow = nrow(md), ncol = 0)
  design <- cbind(condition_effect = as.integer(md$contrast_group == "numerator"), cov_design)
  rownames(design) <- rownames(md)
  colnames(design) <- make.names(colnames(design))
  similarity_file <- cfg_get(config, "method_options.dcats.similarity_matrix")
  similarity <- if (is.null(similarity_file)) NULL else read_similarity_matrix(similarity_file)
  reference_labels <- as_chr(cfg_get(config, "method_options.dcats.reference_cell_types", list()))
  reference <- if (length(reference_labels)) match(reference_labels, colnames(cm)) else NULL
  if (length(reference) && anyNA(reference)) stop("DCATS reference cell types are absent")
  result <- DCATS::dcats_GLM(
    cm, design, similarity_mat = similarity,
    base_model = if (ncol(design) > 1L) "FULL" else "NULL", reference = reference
  )
  raw <- data.frame(
    cell_type = rownames(result$ceoffs), effect = result$ceoffs[, "condition_effect"],
    standard_error = result$coeffs_err[, "condition_effect"], likelihood_ratio = result$LR_vals[, "condition_effect"],
    p_value = result$LRT_pvals[, "condition_effect"], fdr = result$fdr[, "condition_effect"], stringsAsFactors = FALSE
  )
  write_tsv(raw, file.path(task_dir, "dcats_results.tsv"))
  standard_result("dcats", comparison, raw$cell_type, raw$effect, "beta_binomial_log_odds",
    raw$p_value, raw$fdr, raw$effect - 1.96 * raw$standard_error, raw$effect + 1.96 * raw$standard_error,
    evidence_type = "likelihood_ratio_bh_fdr", significant = !is.na(raw$fdr) & raw$fdr <= fdr,
    reference = if (length(reference_labels)) paste(reference_labels, collapse = ";") else "total_denominator",
    note = if (is.null(similarity)) "DCATS without misclassification correction." else "DCATS with configured similarity-matrix correction.")
}

find_sccoda_python <- function() {
  configured <- cfg_get(config, "runtime.sccoda_python")
  if (!is.null(configured) && file.exists(configured)) return(normalizePath(configured))
  pixi_root <- path.expand(cfg_get(config, "runtime.pixi_root", "~/projects/scrna_envs"))
  candidate <- file.path(pixi_root, "07-cell-abundance", ".pixi", "envs", "sccoda", "bin", "python")
  if (file.exists(candidate)) return(normalizePath(candidate))
  stop("scCODA executable not found; set runtime.sccoda_python")
}

sccoda_candidates <- function(md) {
  props <- proportion_matrix[, rownames(md), drop = FALSE]
  presence <- rowMeans(count_matrix[, rownames(md), drop = FALSE] > 0)
  dispersion <- apply(props, 1, function(x) stats::var(x) / max(mean(x), .Machine$double.eps))
  eligible <- names(sort(dispersion[presence >= 1 - as.numeric(cfg_get(config, "method_options.sccoda.automatic_reference_absence_threshold", 0.05))]))
  head(eligible, as.integer(cfg_get(config, "method_options.sccoda.sensitivity_reference_count", 2L)))
}

run_sccoda <- function(comparison, task_dir) {
  md <- sample_meta_for(comparison)
  python <- find_sccoda_python()
  engine_candidates <- c(
    file.path(dirname(normalizePath(script_file)), "cell_abundance_sccoda.py"),
    file.path(dirname(dirname(normalizePath(script_file))), "python", "cell_abundance_sccoda.py")
  )
  engine <- engine_candidates[file.exists(engine_candidates)][1]
  if (!file.exists(engine)) stop("Missing scCODA engine: ", engine)
  counts_file <- file.path(task_dir, "sccoda_input_counts.tsv")
  metadata_file <- file.path(task_dir, "sccoda_input_metadata.tsv")
  cm <- t(count_matrix[, rownames(md), drop = FALSE])
  write_tsv(data.frame(sample = rownames(cm), cm, check.names = FALSE), counts_file)
  meta_write <- md
  meta_write$sample <- rownames(meta_write)
  write_tsv(meta_write[c("sample", condition_col, covariates)], metadata_file)
  refs <- as_chr(cfg_get(config, "method_options.sccoda.reference_cell_types", list("automatic")))
  if ("automatic" %in% refs && bool_value(cfg_get(config, "method_options.sccoda.reference_sensitivity", TRUE), TRUE)) refs <- unique(c(refs, sccoda_candidates(md)))
  statuses <- list(); results <- list()
  for (i in seq_along(refs)) {
    reference <- refs[i]
    ref_dir <- file.path(task_dir, paste0("reference_", safe_name(reference)))
    dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    argv <- c(
      engine, "--counts", counts_file, "--metadata", metadata_file, "--condition", condition_col,
      "--numerator", comparison$numerator, "--denominator", comparison$denominator,
      "--comparison-id", comparison$id, "--reference", reference,
      "--covariates", paste(covariates, collapse = ","), "--fdr", as.character(fdr),
      "--num-samples", as.character(cfg_get(config, "method_options.sccoda.num_samples", 10000L)),
      "--num-warmup", as.character(cfg_get(config, "method_options.sccoda.num_warmup", 1000L)),
      "--seed", as.character(cfg_get(config, "analysis.random_seed", 1L)),
      "--absence-threshold", as.character(cfg_get(config, "method_options.sccoda.automatic_reference_absence_threshold", 0.05)),
      "--output-dir", ref_dir
    )
    if (bool_value(cfg_get(config, "method_options.sccoda.save_posterior", FALSE), FALSE)) argv <- c(argv, "--save-posterior")
    log_file <- file.path(ref_dir, "run.log")
    status <- system2(python, args = shQuote(argv), stdout = log_file, stderr = log_file)
    output_file <- file.path(ref_dir, "all_results.tsv")
    if (identical(status, 0L) && file.exists(output_file)) {
      result <- utils::read.delim(output_file, check.names = FALSE, stringsAsFactors = FALSE)
      for (column in intersect(c("credible", "significant", "tested"), colnames(result))) {
        result[[column]] <- tolower(as.character(result[[column]])) == "true"
      }
      result$reference_run <- reference
      results[[reference]] <- result
      statuses[[i]] <- data.frame(reference, status = "completed", message = "", stringsAsFactors = FALSE)
    } else {
      message <- paste(readLines(log_file, warn = FALSE), collapse = " | ")
      statuses[[i]] <- data.frame(reference, status = "failed", message = substr(message, 1, 2000), stringsAsFactors = FALSE)
    }
  }
  write_tsv(do.call(rbind, statuses), file.path(task_dir, "reference_status.tsv"))
  if (!length(results)) stop("All scCODA reference runs failed; inspect reference_status.tsv")
  combined <- rbind_fill(results)
  combined
}

run_milo <- function(comparison, task_dir) {
  require_pkg("miloR"); require_pkg("SingleCellExperiment"); require_pkg("S4Vectors"); require_pkg("Matrix"); require_pkg("Seurat")
  if (is.null(object) || is.null(cell_meta)) stop("Milo requires a Seurat object with cell-level metadata")
  md <- sample_meta_for(comparison)
  cells <- rownames(cell_meta)[as.character(cell_meta[[sample_col]]) %in% rownames(md)]
  reduction <- cfg_get(config, "method_options.milo.reduction", required = TRUE)
  if (!reduction %in% names(object@reductions)) stop("Milo reduction not found: ", reduction)
  embedding <- Seurat::Embeddings(object, reduction = reduction)[cells, , drop = FALSE]
  d <- min(as.integer(cfg_get(config, "method_options.milo.d", min(30L, ncol(embedding)))), ncol(embedding))
  k <- as.integer(cfg_get(config, "method_options.milo.k", 21L))
  if (nrow(embedding) <= k) stop("Milo k must be smaller than the number of selected cells")
  coldata <- data.frame(
    sample_internal = as.character(cell_meta[cells, sample_col, drop = TRUE]),
    cell_type = as.character(cell_meta[cells, cell_type_col, drop = TRUE]), row.names = cells, stringsAsFactors = FALSE
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(0, nrow = 1, ncol = length(cells), sparse = TRUE)),
    colData = S4Vectors::DataFrame(coldata)
  )
  colnames(sce) <- cells
  SingleCellExperiment::reducedDim(sce, "PCA") <- embedding[, seq_len(d), drop = FALSE]
  milo <- miloR::Milo(sce)
  milo <- miloR::buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  milo <- miloR::makeNhoods(milo, prop = as.numeric(cfg_get(config, "method_options.milo.prop", 0.1)), k = k, d = d, refined = TRUE, reduced_dims = "PCA")
  milo <- miloR::countCells(milo, samples = "sample_internal", meta.data = coldata)
  milo <- miloR::calcNhoodDistance(milo, d = d, reduced.dim = "PCA")
  formula <- stats::as.formula(paste("~ 0 +", paste(c("contrast_group", covariates), collapse = " + ")))
  result <- miloR::testNhoods(
    milo, design = formula, design.df = md,
    model.contrasts = "contrast_groupnumerator - contrast_groupdenominator",
    fdr.weighting = cfg_get(config, "method_options.milo.fdr_weighting", "graph-overlap"),
    reduced.dim = "PCA", robust = TRUE
  )
  annotated <- miloR::annotateNhoods(milo, result, coldata_col = "cell_type")
  annotated$nhood_id <- rownames(annotated)
  write_tsv(annotated, file.path(task_dir, "milo_neighborhood_results.tsv"))
  if (bool_value(cfg_get(config, "method_options.milo.save_object", TRUE), TRUE)) saveRDS(milo, file.path(task_dir, "milo_object.rds"))
  adj_col <- if ("SpatialFDR" %in% colnames(annotated)) "SpatialFDR" else "FDR"
  label_col <- if ("cell_type" %in% colnames(annotated)) "cell_type" else grep("cell_type", colnames(annotated), value = TRUE)[1]
  standardized <- standard_result("milo", comparison, annotated$nhood_id, annotated$logFC, "log2_neighborhood_count_fold_change",
    annotated$PValue, annotated[[adj_col]], evidence_type = paste0("negative_binomial_glm_", adj_col),
    significant = !is.na(annotated[[adj_col]]) & annotated[[adj_col]] <= fdr,
    annotation = annotated[[label_col]], note = paste0("KNN neighborhoods from explicit reduction '", reduction, "'; k=", k, ", d=", d, "."))
  standardized$feature_type <- "neighborhood"
  try({
    p <- miloR::plotDAbeeswarm(annotated, group.by = label_col, alpha = fdr)
    ggplot2::ggsave(file.path(task_dir, "milo_da_beeswarm.pdf"), p, width = 9, height = 6, units = "in")
  }, silent = TRUE)
  try({
    milo <- miloR::buildNhoodGraph(milo)
    p <- miloR::plotNhoodGraphDA(milo, annotated, alpha = fdr)
    ggplot2::ggsave(file.path(task_dir, "milo_da_graph.pdf"), p, width = 8, height = 7, units = "in")
  }, silent = TRUE)
  standardized
}

plot_base_abundance <- function() {
  require_pkg("ggplot2")
  plot_long <- long_counts
  plot_long$proportion <- plot_long$n_cells / denominator_totals[as.character(plot_long[[sample_col]])]
  sample_order <- as.character(sample_meta[[sample_col]])[order(as.character(sample_meta[[condition_col]]), as.character(sample_meta[[sample_col]]))]
  plot_long[[sample_col]] <- factor(plot_long[[sample_col]], levels = sample_order)
  pdf(file.path(out, "sample_composition.pdf"), width = 11, height = 7, onefile = TRUE)
  for (chunk in split(sample_order, ceiling(seq_along(sample_order) / 24))) {
    p <- ggplot2::ggplot(plot_long[plot_long[[sample_col]] %in% chunk, ], ggplot2::aes(x = .data[[sample_col]], y = .data$proportion, fill = .data$cell_type)) +
      ggplot2::geom_col(width = 0.85) + ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(title = "Sample-level cell composition", subtitle = denominator_description, x = NULL, y = "Relative abundance", fill = cell_type_col) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    print(p)
  }
  dev.off()
  pdf(file.path(out, "cell_type_proportions_by_condition.pdf"), width = 10, height = 7, onefile = TRUE)
  for (chunk in split(types, ceiling(seq_along(types) / 12))) {
    p <- ggplot2::ggplot(plot_long[plot_long$cell_type %in% chunk, ], ggplot2::aes(x = .data[[condition_col]], y = .data$proportion, colour = .data[[condition_col]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, colour = "grey55", width = 0.5) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1), size = 2) +
      ggplot2::geom_text(ggplot2::aes(label = .data[[sample_col]]), position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1), vjust = -0.7, size = 2.2, show.legend = FALSE) +
      ggplot2::facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
      ggplot2::labs(title = "Per-sample relative abundance", subtitle = "Each point is one biological sample", x = NULL, y = "Proportion", colour = condition_col) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
    print(p)
  }
  dev.off()
  heat <- plot_long
  heat$logit_proportion <- qlogis((heat$n_cells + 0.5) / (denominator_totals[as.character(heat[[sample_col]])] + 1))
  pdf(file.path(out, "sample_proportion_heatmap.pdf"), width = 11, height = 7, onefile = TRUE)
  for (chunk in split(types, ceiling(seq_along(types) / 25))) {
    p <- ggplot2::ggplot(heat[heat$cell_type %in% chunk, ], ggplot2::aes(x = .data[[sample_col]], y = .data$cell_type, fill = .data$logit_proportion)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.15) + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
      ggplot2::labs(title = "Sample × cell-type composition", subtitle = "Logit proportion with 0.5-cell continuity correction", x = NULL, y = NULL, fill = "logit(prop)") +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    print(p)
  }
  dev.off()
}

plot_method_effects <- function(result, task_dir, title) {
  require_pkg("ggplot2")
  if (!nrow(result)) return(invisible(NULL))
  result$display <- ifelse(result$feature_type == "neighborhood", paste0(result$annotation, " / ", result$feature_id), result$feature_id)
  ordering <- order(ifelse(is.na(result$adjusted_p_value), 1, result$adjusted_p_value), -abs(result$effect))
  result <- result[ordering, , drop = FALSE]
  pdf(file.path(task_dir, "effect_summary.pdf"), width = 9, height = 7, onefile = TRUE)
  for (indices in split(seq_len(nrow(result)), ceiling(seq_len(nrow(result)) / 30))) {
    d <- result[indices, , drop = FALSE]
    d$display <- factor(d$display, levels = rev(d$display))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = display, colour = significant)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, orientation = "y", na.rm = TRUE) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
      ggplot2::labs(title = title, subtitle = unique(d$effect_scale)[1], x = "Method-specific effect", y = NULL, colour = paste0("Evidence <= ", fdr)) +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "top")
    print(p)
  }
  dev.off()
}

plot_base_abundance()

comparisons <- cfg_get(config, "comparisons", required = TRUE)
all_results <- list(); status_rows <- list(); task_index <- 0L
for (comparison in comparisons) {
  comparison$id <- as.character(comparison$id)
  comparison$numerator <- as.character(comparison$numerator)
  comparison$denominator <- as.character(comparison$denominator)
  for (method in methods) {
    task_index <- task_index + 1L
    task_dir <- file.path(out, "comparisons", safe_name(comparison$id), method)
    dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
    started <- Sys.time()
    task <- tryCatch({
      result <- switch(method,
        propeller = run_propeller(comparison, task_dir),
        sccomp = run_sccomp(comparison, task_dir),
        sccoda = run_sccoda(comparison, task_dir),
        dcats = run_dcats(comparison, task_dir),
        milo = run_milo(comparison, task_dir),
        stop("Unsupported method: ", method)
      )
      write_tsv(result, file.path(task_dir, "all_results.tsv"))
      write_tsv(result[result$significant %in% TRUE, , drop = FALSE], file.path(task_dir, "significant_results.tsv"))
      plot_method_effects(result, task_dir, paste(method, comparison$id))
      key <- paste(comparison$id, method, sep = "__")
      all_results[[key]] <- result
      data.frame(comparison_id = comparison$id, method, status = "completed", n_results = nrow(result), n_significant = sum(result$significant %in% TRUE, na.rm = TRUE), message = "", duration_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), stringsAsFactors = FALSE)
    }, error = function(e) {
      message <- conditionMessage(e)
      writeLines(message, file.path(task_dir, "ERROR.txt"))
      data.frame(comparison_id = comparison$id, method, status = classify_failure(message), n_results = 0L, n_significant = 0L, message, duration_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), stringsAsFactors = FALSE)
    })
    status_rows[[task_index]] <- task
  }
}

status <- do.call(rbind, status_rows)
write_tsv(status, file.path(out, "task_status.tsv"))
if (length(all_results)) {
  combined <- rbind_fill(all_results)
  write_tsv(combined, file.path(out, "all_method_results.tsv"))
  write_tsv(combined[combined$significant %in% TRUE, , drop = FALSE], file.path(out, "significant_method_results.tsv"))
  comparable <- combined[combined$feature_type == "cell_type", , drop = FALSE]
  # scCODA sensitivity runs are repeated fits of one method, not independent methods.
  # Use its automatic-reference fit for cross-method concordance when present.
  if ("reference_run" %in% colnames(comparable)) {
    is_sccoda_extra <- comparable$method == "sccoda" & !is.na(comparable$reference_run) & comparable$reference_run != "automatic"
    automatic_keys <- unique(paste(
      comparable$comparison_id[comparable$method == "sccoda" & comparable$reference_run == "automatic"],
      comparable$feature_id[comparable$method == "sccoda" & comparable$reference_run == "automatic"], sep = "\r"
    ))
    row_keys <- paste(comparable$comparison_id, comparable$feature_id, sep = "\r")
    comparable <- comparable[!(is_sccoda_extra & row_keys %in% automatic_keys), , drop = FALSE]
  }
  if (nrow(comparable)) {
    keys <- unique(paste(comparable$comparison_id, comparable$feature_id, sep = "\r"))
    concordance <- do.call(rbind, lapply(keys, function(key) {
      parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
      d <- comparable[comparable$comparison_id == parts[1] & comparable$feature_id == parts[2], , drop = FALSE]
      # One vote per method and feature. This prevents any sensitivity fit from
      # being counted as another source of statistical support.
      d <- d[!duplicated(d$method), , drop = FALSE]
      sig <- d[d$significant %in% TRUE, , drop = FALSE]
      sig_signs <- unique(sign(sig$effect[is.finite(sig$effect) & sig$effect != 0]))
      n_significant_methods <- length(unique(sig$method))
      classification <- if (length(sig_signs) > 1L) "discordant" else if (n_significant_methods >= 2L) "supported" else if (n_significant_methods == 1L) "partial" else "not_supported"
      data.frame(comparison_id = parts[1], cell_type = parts[2], n_methods_tested = length(unique(d$method)), n_methods_significant = n_significant_methods, significant_methods = paste(unique(sig$method), collapse = ";"), significant_directions = paste(ifelse(sig$effect > 0, "up", "down"), collapse = ";"), classification, note = "Method-specific effects are not numerically interchangeable; scCODA sensitivity references do not add method votes.", stringsAsFactors = FALSE)
    }))
    write_tsv(concordance, file.path(out, "method_concordance.tsv"))
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
artifacts <- list.files(out, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[file.info(artifacts)$isdir %in% FALSE]
write_run_manifest(config, "13-scrna-test-cell-abundance", out, artifacts,
  c(paste0("methods=", paste(methods, collapse = ",")), paste0("denominator=", denominator_description),
    "Biological samples, not cells, are the inferential replicates.",
    "Observed proportions are relative compositions, not absolute tissue cell abundance.",
    "Input object was not rewritten."))
if (!any(status$status == "completed")) stop("No cell-abundance task completed; inspect task_status.tsv")
