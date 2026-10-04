source("toolkit/R/runtime.R")
root <- tempfile("executor-preflight-"); dir.create(root)
run_rejected <- function(script, config, expected) {
  path <- file.path(root, "config.json")
  jsonlite::write_json(config, path, auto_unbox = TRUE)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(path)), stdout = TRUE, stderr = TRUE))
  stopifnot(!is.null(attr(output, "status")), attr(output, "status") != 0L,
            any(grepl(expected, output, fixed = TRUE)), !dir.exists(config$output_dir))
}
base <- list(input = list(object = file.path(root, "missing.rds")),
             metadata = list(condition = "condition", covariates = list("sex")),
             output_dir = file.path(root, "results"))
config <- base; config$analysis <- list(min_samples_per_group = 1)
run_rejected("toolkit/R/differential_analysis.R", config, "analysis.min_samples_per_group")
config <- base; config$analysis <- list(design = "~ (sex + condition)^2")
run_rejected("toolkit/R/differential_analysis.R", config, "Unsupported analysis.design")

# Fake object bytes suffice: rejected approval must be checked before loading it.
object <- file.path(root, "prepared.qs"); writeLines("object fixture", object)
table <- file.path(root, "decisions.tsv"); writeLines(c("cluster\tbroad\tfine\tdecision", "0\tEndothelial\tEndothelial\tconfirmed"), table)
record <- file.path(root, "record.json")
jsonlite::write_json(list(schema_version = 1L, kind = "annotation_review", review_run_id = "run-1",
  cluster_column = "cluster", apply_object = list(sha256 = .scrna_sha256(object))), record, auto_unbox = TRUE)
config <- list(workflow = list(action = "apply_confirmed"),
  input = list(object = object, decisions = table, review_record = record),
  metadata = list(cluster = "cluster"), output_dir = file.path(root, "results"),
  approval = list(status = "approved", source = "human", review_run_id = "run-1",
    decision_sha256 = .scrna_sha256(table), review_record_sha256 = .scrna_sha256(record)))
bad <- config; bad$approval$source <- "agent"
run_rejected("toolkit/R/annotate_cells.R", bad, "approval.source")
bad <- config; bad$approval$review_run_id <- "another-run"
run_rejected("toolkit/R/annotate_cells.R", bad, "approval.review_run_id")
writeLines("changed table", table)
run_rejected("toolkit/R/annotate_cells.R", config, "decisions SHA256")
config$approval$decision_sha256 <- .scrna_sha256(table)
writeLines("changed object", object)
run_rejected("toolkit/R/annotate_cells.R", config, "input.object SHA256")
writeLines("changed review", record)
run_rejected("toolkit/R/annotate_cells.R", config, "review record SHA256")
unlink(root, recursive = TRUE)
cat("PASS direct R executor rejects unsupported design, replication bypass and stale annotation approval before object loading or scientific output\n")
