args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "tests/fixtures/tiny_scrna.rds"
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
set.seed(20260827)
genes <- c("mt-Nd1", "mt-Co1", "Kdr", "Pecam1", "Cdh5", "Col1a1", "Col3a1", "Dcn", "Actb", "Gapdh", paste0("Gene", 1:10))
samples <- rep(c("ctrl_1", "ctrl_2", "stz_1", "stz_2"), each = 20)
condition <- ifelse(grepl("^ctrl", samples), "control", "stz")
batch <- rep(c("batch_1", "batch_2"), each = 20, times = 2)
cell_type <- rep(rep(c("Endothelial", "Fibroblast"), each = 10), 4)
cells <- sprintf("cell_%03d", seq_along(samples))
mu <- matrix(1.5, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
mu[c("Kdr", "Pecam1", "Cdh5"), cell_type == "Endothelial"] <- 7
mu[c("Col1a1", "Col3a1", "Dcn"), cell_type == "Fibroblast"] <- 7
mu[c("Gene1", "Gene2"), condition == "stz"] <- 5
counts <- matrix(stats::rpois(length(mu), lambda = as.vector(mu)), nrow = nrow(mu), dimnames = dimnames(mu))
obj <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), project = "tiny_scrna")
obj$sample_label <- samples; obj$condition <- condition; obj$batch_id <- batch; obj$cell_type <- cell_type
obj$seurat_clusters <- ifelse(cell_type == "Endothelial", "0", "1")
obj$annotation_marker_pred <- cell_type
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-")
embedding <- cbind(UMAP_1 = rep(c(-1, 1), each = 10, times = 4) + stats::rnorm(length(cells), 0, .2), UMAP_2 = rep(seq(-1, 1, length.out = 20), 4) + stats::rnorm(length(cells), 0, .2))
rownames(embedding) <- cells
obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = embedding, key = "UMAP_", assay = "RNA")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output)
cat(normalizePath(output), "\n")
