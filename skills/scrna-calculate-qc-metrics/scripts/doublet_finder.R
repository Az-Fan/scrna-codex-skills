cpm <- function(counts) {
  sf <- Matrix::colSums(counts) / 1e6
  if (is.matrix(counts)) return(t(t(counts) / sf))
  if (inherits(counts, "dgCMatrix")) {
    sep <- counts@p
    sep <- sep[-1] - sep[-length(sep)]
    j <- S4Vectors::Rle(seq_along(sep), sep)
    counts@x <- counts@x / sf[as.integer(j)]
  } else if (inherits(counts, "dgTMatrix")) {
    counts@x <- counts@x / sf[counts@j + 1]
  } else stop("Unsupported matrix class: ", class(counts)[1])
  counts
}

logCPM <- function(counts) {
  x <- cpm(counts)
  if (is.matrix(x)) x <- log2(x + 1) else x@x <- log2(x@x + 1)
  x
}

rd_PCA <- function(norm.dat, select.genes = rownames(norm.dat), select.cells = colnames(norm.dat), sampled.cells = select.cells, max.pca = 10, th = 2) {
  pca <- stats::prcomp(t(as.matrix(norm.dat[select.genes, sampled.cells])), tol = 0.01)
  v <- summary(pca)$importance[2, ]
  selected <- head(which((v - mean(v)) / sd(v) > th), max.pca)
  if (!length(selected)) return(NULL)
  if (length(sampled.cells) < length(select.cells)) {
    rotation <- pca$rotation[, selected, drop = FALSE]
    reduced <- as.matrix(Matrix::t(norm.dat[rownames(rotation), select.cells, drop = FALSE]) %*% rotation)
  } else reduced <- pca$x[, selected, drop = FALSE]
  list(rd.dat = reduced, pca = pca)
}

sc06SyntheticDoubletScore <- function(data, select.genes, proportion.artificial = 0.20, k = NULL) {
  if (is.null(k)) k <- round(max(10, min(100, ncol(data) * 0.01)))
  real.cells <- colnames(data)
  n_doublets <- round(length(real.cells) / (1 - proportion.artificial) - length(real.cells))
  doublets <- data[, sample(real.cells, n_doublets, TRUE)] + data[, sample(real.cells, n_doublets, TRUE)]
  colnames(doublets) <- paste0("X", seq_len(n_doublets))
  combined <- cbind(data, doublets)
  normalized <- logCPM(combined)
  if (ncol(data) > 10000) {
    sampled <- sample(seq_len(ncol(data)), min(ncol(data), 10000))
    reduced <- rd_PCA(normalized, select.genes, colnames(combined), sampled.cells = sampled, th = 0, max.pca = 50)
  } else reduced <- rd_PCA(normalized, select.genes, colnames(combined), th = 0, max.pca = 50)
  if (is.null(reduced)) stop("PCA did not retain any component for doublet scoring")
  reduced <- reduced$rd.dat
  neighbours <- RANN::nn2(reduced, k = k)
  num <- ncol(data)
  artificial_dist <- RANN::nn2(reduced[seq_len(num), , drop = FALSE], reduced[(num + 1):nrow(reduced), , drop = FALSE], k = 10)[[2]]
  threshold <- mean(artificial_dist) + 1.64 * stats::sd(as.vector(artificial_dist))
  freq <- neighbours[[1]] > num & neighbours[[2]] < threshold
  score <- pmax(rowMeans(freq), rowMeans(freq[, seq_len(ceiling(k / 2)), drop = FALSE]))
  names(score) <- rownames(reduced)
  list(score[real.cells], score[colnames(doublets)])
}
