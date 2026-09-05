args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript apply_qc_filter.R config.json")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "runtime.R"))
config <- read_skill_config(args[[1]])
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required")

approval_status <- tolower(as.character(cfg_get(config, "approval.status", required = TRUE)))
if (!identical(approval_status, "approved")) stop("approval.status must be exactly 'approved'")

object_path <- cfg_get(config, "input.object", required = TRUE)
decision_path <- normalizePath(cfg_get(config, "input.decision_table", required = TRUE), mustWork = TRUE)
assay <- cfg_get(config, "input.assay", "RNA")
sample_col <- cfg_get(config, "metadata.sample", required = TRUE)
condition_col <- cfg_get(config, "metadata.condition")
cell_id_col <- cfg_get(config, "decision.cell_id_column", "cell_id")
include_cols <- unlist(cfg_get(config, "decision.include_all_true", required = TRUE), use.names = FALSE)
exclude_cols <- unlist(cfg_get(config, "decision.exclude_any_true", list()), use.names = FALSE)
reason_col <- cfg_get(config, "decision.reason_column")
carry_cols <- unlist(cfg_get(config, "decision.carry_columns", list()), use.names = FALSE)
expected_retained <- as.integer(cfg_get(config, "decision.expected_retained_cells", required = TRUE))
if (!length(include_cols)) stop("decision.include_all_true must contain at least one column")
if (is.na(expected_retained) || expected_retained < 1L) stop("decision.expected_retained_cells must be positive")

obj <- load_scrna_object(object_path, "auto")
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object")
if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
assert_metadata(obj, c(sample_col, condition_col %||% character()))

decisions <- utils::read.delim(decision_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- unique(c(cell_id_col, include_cols, exclude_cols, reason_col %||% character(), carry_cols))
missing <- setdiff(required, names(decisions))
if (length(missing)) stop("Decision table is missing columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(decisions[[cell_id_col]])) stop("Decision table contains duplicated cell IDs")
if (!setequal(colnames(obj), decisions[[cell_id_col]])) stop("Decision table and object cell sets differ")
decisions <- decisions[match(colnames(obj), decisions[[cell_id_col]]), , drop = FALSE]

as_flag <- function(x, name) {
  if (is.logical(x)) {
    if (anyNA(x)) stop("Decision column contains missing values: ", name)
    return(x)
  }
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y)); out[y %in% c("true", "1", "yes")] <- TRUE; out[y %in% c("false", "0", "no")] <- FALSE
  if (anyNA(out)) stop("Decision column is not unambiguous boolean data: ", name)
  out
}

include_matrix <- vapply(include_cols, function(column) as_flag(decisions[[column]], column), logical(nrow(decisions)))
if (is.null(dim(include_matrix))) include_matrix <- matrix(include_matrix, ncol = 1L)
keep <- rowSums(include_matrix) == ncol(include_matrix)
if (length(exclude_cols)) {
  exclude_matrix <- vapply(exclude_cols, function(column) as_flag(decisions[[column]], column), logical(nrow(decisions)))
  if (is.null(dim(exclude_matrix))) exclude_matrix <- matrix(exclude_matrix, ncol = 1L)
  keep <- keep & rowSums(exclude_matrix) == 0L
}
if (sum(keep) != expected_retained) stop("Retained-cell count differs from approved expectation: actual=", sum(keep), ", expected=", expected_retained)

reason <- if (!is.null(reason_col)) as.character(decisions[[reason_col]]) else ifelse(keep, "retained", "excluded_by_approved_decision")
reason[keep] <- "retained"
decision_output <- decisions
names(decision_output)[names(decision_output) == cell_id_col] <- "cell_id"
decision_output$retained <- keep
decision_output$filter_decision_reason <- reason
decision_output[[sample_col]] <- as.character(obj[[]][[sample_col]])
if (!is.null(condition_col)) decision_output[[condition_col]] <- as.character(obj[[]][[condition_col]])

meta <- obj[[]][keep, , drop = FALSE]
if (length(carry_cols)) meta[, carry_cols] <- decisions[keep, carry_cols, drop = FALSE]
meta$qc_filter_status <- "approved"
meta$qc_filter_source <- decision_path
counts <- get_raw_counts(obj, assay = assay)[, keep, drop = FALSE]
filtered <- Seurat::CreateSeuratObject(counts = counts, assay = assay, project = cfg_get(config, "project.id", required = TRUE), meta.data = meta, min.cells = 0, min.features = 0)
if (ncol(filtered) != expected_retained || nrow(filtered) != nrow(obj)) stop("Filtered object dimensions are invalid")
if (length(filtered@reductions) || length(filtered@graphs)) stop("Filtered handoff unexpectedly inherited reductions or graphs")

out <- prepare_output(config)
object_name <- cfg_get(config, "output.object_name", "filtered_object.qs")
if (!grepl("\\.(qs|rds)$", object_name, ignore.case = TRUE)) stop("output.object_name must end in .qs or .rds")
object_out <- file.path(out, object_name)
if (file.exists(object_out)) stop("Refusing to overwrite filtered object: ", object_out)
extension <- if (grepl("\\.qs$", object_out, ignore.case = TRUE)) ".qs" else ".rds"
partial <- sub(paste0("\\", extension, "$"), paste0(".partial", extension), object_out, ignore.case = TRUE)
on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
save_scrna_object(filtered, partial)
check <- load_scrna_object(partial, if (identical(extension, ".qs")) "qs" else "rds")
if (!identical(dim(check), dim(filtered)) || !identical(colnames(check), colnames(filtered))) stop("Filtered object reread validation failed")
if (!file.rename(partial, object_out)) stop("Atomic filtered-object finalization failed")

decision_file <- gzfile(file.path(out, "cell_filter_decisions.tsv.gz"), "wt")
utils::write.table(decision_output, decision_file, sep = "\t", quote = FALSE, row.names = FALSE); close(decision_file)
sample <- as.character(obj[[]][[sample_col]])
condition <- if (!is.null(condition_col)) as.character(obj[[]][[condition_col]]) else rep(NA_character_, ncol(obj))
groups <- if (!is.null(condition_col)) list(sample = sample, condition = condition) else list(sample = sample)
retention <- aggregate(cbind(n_input = rep(1L, length(keep)), n_retained = as.integer(keep)), groups, sum)
retention$retention_fraction <- retention$n_retained / retention$n_input
sample_file <- file.path(out, "filter_summary_by_sample.tsv")
utils::write.table(retention, sample_file, sep = "\t", quote = FALSE, row.names = FALSE)
condition_file <- NULL
if (!is.null(condition_col)) {
  condition_summary <- aggregate(cbind(n_input = rep(1L, length(keep)), n_retained = as.integer(keep)), list(condition = condition), sum)
  condition_summary$retention_fraction <- condition_summary$n_retained / condition_summary$n_input
  condition_file <- file.path(out, "filter_summary_by_condition.tsv")
  utils::write.table(condition_summary, condition_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
reason_counts <- as.data.frame(table(filter_decision_reason = decision_output$filter_decision_reason), stringsAsFactors = FALSE)
names(reason_counts)[names(reason_counts) == "Freq"] <- "n_cells"
reason_file <- file.path(out, "filter_decision_counts.tsv")
utils::write.table(reason_counts, reason_file, sep = "\t", quote = FALSE, row.names = FALSE)
approval_file <- file.path(out, "approved_filter_record.tsv")
approval <- data.frame(status = approval_status, approved_at = cfg_get(config, "approval.approved_at", ""), note = cfg_get(config, "approval.note", ""), expected_retained_cells = expected_retained, decision_table = decision_path, decision_table_sha256 = .scrna_sha256(decision_path), stringsAsFactors = FALSE)
utils::write.table(approval, approval_file, sep = "\t", quote = FALSE, row.names = FALSE)
session_file <- technical_path(out, "session_info.txt"); writeLines(capture.output(utils::sessionInfo()), session_file)
artifacts <- c(object_out, file.path(out, "cell_filter_decisions.tsv.gz"), sample_file, condition_file %||% character(), reason_file, approval_file, session_file)
write_run_manifest(config, "04-scrna-apply-qc-filter", out, artifacts, c("Scientific decisions were read from the approved decision table", paste0("retained_cells=", expected_retained), "Output object contains raw counts and retained metadata only"))
