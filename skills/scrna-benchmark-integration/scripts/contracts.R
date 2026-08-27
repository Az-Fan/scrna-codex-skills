assert_metadata_contract <- function(meta, sample_col, condition_col, batch_col = NULL) {
  required <- c(sample_col, condition_col, batch_col)
  required <- required[!is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(meta))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  if (anyNA(meta[[sample_col]]) || any(meta[[sample_col]] == "")) stop("Sample identifiers contain missing values")
  invisible(TRUE)
}

sample_cell_counts <- function(meta, sample_col, group_col) {
  as.data.frame(table(sample = meta[[sample_col]], group = meta[[group_col]]), stringsAsFactors = FALSE)
}

assert_replicates <- function(meta, sample_col, condition_col, minimum = 2L) {
  design <- unique(meta[c(sample_col, condition_col)])
  counts <- table(design[[condition_col]])
  if (any(counts < minimum)) warning("At least one condition has fewer than ", minimum, " independent samples")
  invisible(counts)
}
