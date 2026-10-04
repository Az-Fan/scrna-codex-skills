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
legacy <- obj
legacy[["RNA_v4"]] <- SeuratObject::CreateAssayObject(counts = expected)
stopifnot(all(get_raw_counts(legacy, "RNA_v4") == expected))
source("toolkit/R/differential_utils.R")
stopifnot(identical(join_assay_layers(legacy, "RNA_v4"), legacy))
declared <- list(analysis = list(counts_source = list(kind = "raw_umi")))
stopifnot(assert_pseudobulk_counts_source(legacy, "RNA_v4", declared)$counts_kind == "raw_umi")
failure <- tryCatch(assert_pseudobulk_counts_source(obj, "RNA", list()), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("counts_source.kind", failure))
corrected <- obj
corrected[["integrated"]] <- SeuratObject::CreateAssay5Object(counts = round(log1p(expected)))
failure <- tryCatch(assert_pseudobulk_counts_source(corrected, "integrated", declared), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("cannot use", failure))
reference_path <- tempfile(fileext = ".rds"); saveRDS(obj, reference_path)
declared$analysis$counts_source$reference_object <- reference_path
stopifnot(assert_pseudobulk_counts_source(selected, "RNA", declared)$reference_sha256 == .scrna_sha256(reference_path))
Sys.setenv(SCRNA_TEST_COUNTS_REFERENCE = reference_path)
declared$analysis$counts_source$reference_object <- "${SCRNA_TEST_COUNTS_REFERENCE}"
stopifnot(assert_pseudobulk_counts_source(obj, "RNA", declared)$reference_object == normalizePath(reference_path))
Sys.unsetenv("SCRNA_TEST_COUNTS_REFERENCE")
declared$analysis$counts_source$reference_object <- reference_path
altered <- obj
changed_counts <- expected; changed_counts@x[[1]] <- changed_counts@x[[1]] + 1
SeuratObject::LayerData(altered, assay = "RNA", layer = "counts") <- changed_counts
failure <- tryCatch(assert_pseudobulk_counts_source(altered, "RNA", declared), error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("differ from", failure))
unlink(reference_path)
cat("PASS v4/v5 counts, split-layer/subset alignment, raw-source declaration, transformed-assay rejection and reference mismatch\n")
