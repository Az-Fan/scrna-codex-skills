source("toolkit/R/runtime.R")
source("toolkit/R/differential_utils.R")
root <- tempfile("scrna-runtime-contracts-")
dir.create(root)
directory <- file.path(root, "tenx"); dir.create(directory)
writeLines("counts", file.path(directory, "matrix.mtx"))
writeLines("gene", file.path(directory, "features.tsv"))
first <- .scrna_sha256(directory)
copy <- file.path(root, "copy"); dir.create(copy)
file.copy(list.files(directory, full.names = TRUE), copy)
stopifnot(identical(first, .scrna_sha256(copy)))
writeLines("changed counts", file.path(copy, "matrix.mtx"))
stopifnot(!identical(first, .scrna_sha256(copy)))
config_path <- file.path(root, "config.json")
jsonlite::write_json(list(input = list(path = directory)), config_path, auto_unbox = TRUE)
record <- attr(read_skill_config(config_path), "input_record")
stopifnot(record$type == "directory", length(record$members) == 2L, identical(record$sha256, first))

table <- data.frame(gene = c("Kdr", "Pecam1"), log2FoldChange = c(1, 2), padj = c(.01, .02),
  population = "EC", comparison_id = "case_vs_control", numerator = "case", denominator = "control")
paths <- file.path(root, c("a.tsv", "b.tsv"))
write_tsv(table, paths[1]); write_tsv(table, paths[2])
tasks <- make_table_tasks(list(input = list(differential_tables = list(list(path = paths[1], id = "a"), list(path = paths[2], id = "b")))))
jsonlite::write_json(list(input = list(differential_tables = list(list(path = paths[1], id = "a"), list(path = paths[2], id = "b")))), config_path, auto_unbox = TRUE)
multi_record <- attr(read_skill_config(config_path), "input_record")
stopifnot(multi_record$type == "tables", length(multi_record$tables) == 2L,
  identical(multi_record$tables[[1]]$sha256, .scrna_sha256(paths[1])))
stopifnot(length(tasks) == 2L, !anyDuplicated(vapply(tasks, function(x) x$id, "")),
  all(vapply(tasks, function(x) x$comparison$numerator == "case" && x$comparison$denominator == "control", logical(1))))
# Exercise the real multi-task writer with controlled nonempty computation
# outputs, independently of whether the tiny genes have significant GO terms.
table$log2FoldChange <- -table$log2FoldChange
write_tsv(table, paths[2])
workflow_config <- list(project = list(id = "contract_test"), output_dir = file.path(root, "enrichment"),
  input = list(differential_tables = list(list(path = paths[1], id = "a"), list(path = paths[2], id = "b"))))
jsonlite::write_json(workflow_config, config_path, auto_unbox = TRUE)
run_enrichment <- function(result, task_dir, population, comparison, config) {
  data.frame(term = "controlled_output", effect = result$log2FoldChange[1], population = population,
    numerator = comparison$numerator, denominator = comparison$denominator)
}
summarize_enrichment_status <- function(task_dir) "completed"
run_enrichment_only_workflow(read_skill_config(config_path))
combined <- utils::read.delim(file.path(workflow_config$output_dir, "enrichment_all_comparisons.tsv"))
stopifnot(nrow(combined) == 2L, length(unique(combined$task_id)) == 2L,
  setequal(combined$effect, c(-1, 1)), setequal(combined$input_table, paths))
failure <- tryCatch(normalize_comparisons(list(comparisons = list(
  list(id = "T cell", numerator = "case", denominator = "control"),
  list(id = "T/cell", numerator = "other", denominator = "control")))), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("unique after sanitization", failure))
cat("PASS content-addressed directories, independent table tasks, direction handoff and comparison collision rejection\n")
unlink(root, recursive = TRUE)
