source("toolkit/R/runtime.R")
args <- commandArgs(trailingOnly = TRUE)
fixture <- if (length(args)) args[[1]] else "tests/fixtures/tiny_scrna.rds"
obj <- readRDS(fixture)
expected <- get_raw_counts(obj)
split_obj <- obj
split_obj[["RNA"]] <- split(split_obj[["RNA"]], f = split_obj$sample_label)
original <- split_obj
observed <- get_raw_counts(split_obj)
stopifnot(identical(dimnames(observed), dimnames(expected)), all(observed == expected), identical(split_obj, original))
selected <- subset(split_obj, cells = colnames(split_obj)[split_obj$cell_type == "Endothelial"])
exported <- get_raw_counts(selected)
stopifnot(ncol(exported) == 40L, identical(colnames(exported), colnames(selected)), all(exported == expected[, colnames(selected)]))
bad <- obj
bad[["partial"]] <- SeuratObject::CreateAssay5Object(counts = expected[, 1:20])
failure <- tryCatch(get_raw_counts(bad, "partial"), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("do not match all object cells", failure))
fractional <- obj
counts <- expected
counts@x[[1]] <- counts@x[[1]] + 0.5
SeuratObject::LayerData(fractional, assay = "RNA", layer = "counts") <- counts
failure <- tryCatch(get_raw_counts(fractional), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("non-negative integers", failure))
cat("PASS complete multi-layer counts, subset alignment, unchanged source, partial-assay and noninteger rejection\n")
