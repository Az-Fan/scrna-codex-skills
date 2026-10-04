source("toolkit/R/runtime.R")
root <- tempfile("de-inference-modes-"); dir.create(root)
set.seed(31)
counts <- matrix(rnbinom(60 * 40, mu = 10, size = 3), nrow = 60,
  dimnames = list(paste0("gene", 1:60), paste0("cell", 1:40)))
counts[1:8, 21:40] <- counts[1:8, 21:40] + 15L
object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
object$sample <- rep(c("sample_A", "sample_B"), each = 20)
object$condition <- rep(c("A", "B"), each = 20)
object <- Seurat::NormalizeData(object, verbose = FALSE)
input <- file.path(root, "object.rds"); saveRDS(object, input)
config <- list(project = list(id = "inference_modes"), input = list(object = input),
  metadata = list(sample = "sample", condition = "condition"),
  comparison = list(id = "B_vs_A", numerator = "B", denominator = "A"),
  analysis = list(method = "seurat_wilcox", min_cells_per_sample_population = 3L,
    counts_source = list(kind = "raw_umi"), min_samples_per_group = 2L))
run_case <- function(label, config, succeeds) {
  config$output_dir <- file.path(root, label)
  path <- file.path(root, paste0(label, ".json"))
  jsonlite::write_json(config, path, auto_unbox = TRUE)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c("toolkit/R/differential_analysis.R", shQuote(path)), stdout = TRUE, stderr = TRUE))
  exit <- attr(output, "status") %||% 0L
  if ((exit == 0L) != succeeds) stop(paste(output, collapse = "\n"))
  status <- read.delim(file.path(config$output_dir, "task_status.tsv"))
  if (!succeeds) {
    stopifnot(status$status == "skipped_low_replicates", !file.exists(file.path(config$output_dir, "all_comparisons.tsv")))
  } else {
    result <- read.delim(file.path(config$output_dir, "all_comparisons.tsv"))
    stopifnot(nrow(result) > 0L, all(result$inference_level == "cell_level_exploratory"),
      all(result$n_samples_numerator == 1L), all(result$n_samples_denominator == 1L),
      all(is.na(result$wald_tested)))
  }
}
run_case("exploratory", config, TRUE)
formal <- config; formal$analysis$method <- "pseudobulk_deseq2"
run_case("formal", formal, FALSE)
single <- subset(object, cells = colnames(object)[object$condition == "A"])
single_input <- file.path(root, "single.rds"); saveRDS(single, single_input)
config$input$object <- single_input
run_case("missing_condition", config, FALSE)
unlink(root, recursive = TRUE)
cat("PASS one-sample exploratory DE remains explicitly cell-level; formal DE and missing-condition cases remain blocked\n")
