as_chr <- function(x) if (is.null(x)) character() else as.character(unlist(x, use.names = FALSE))
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
rbind_fill <- function(xs) {
  if (!length(xs)) return(data.frame())
  columns <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) { for (nm in setdiff(columns, names(x))) x[[nm]] <- NA; x[, columns, drop = FALSE] })
  do.call(rbind, xs)
}

first_column <- function(x, explicit, candidates, required = FALSE) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (!explicit %in% names(x)) stop("Configured table column not found: ", explicit)
    return(explicit)
  }
  hit <- candidates[candidates %in% names(x)][1]
  if (is.na(hit) && required) stop("Could not detect required column; tried: ", paste(candidates, collapse = ", "))
  if (is.na(hit)) NULL else hit
}

read_de_table <- function(path, sheet = NULL) {
  path <- normalizePath(path, mustWork = TRUE); ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for Excel differential tables")
    return(as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1), check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (ext == "csv") return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

standardize_de_table <- function(x, config) {
  cols <- cfg_get(config, "enrichment.table_columns", list())
  gene_col <- first_column(x, cols$gene, c("gene", "gene_symbol", "symbol", "SYMBOL", "Gene", "genes"), TRUE)
  lfc_col <- first_column(x, cols$log2fc, c("log2FoldChange", "avg_log2FC", "avg_logFC", "logFC", "log2FC", "avg_diff"))
  p_col <- first_column(x, cols$pvalue, c("pvalue", "p_val", "PValue", "P.Value", "p_value", "pval"))
  padj_col <- first_column(x, cols$padj, c("padj", "p_val_adj", "FDR", "adj.P.Val", "p_adjust", "qvalue"))
  stat_col <- first_column(x, cols$stat, c("stat", "WaldStatistic", "t", "F", "score", "signed_stat"))
  sig_col <- first_column(x, cols$significance, c("significance", "threshold"))
  if (is.null(lfc_col) && is.null(stat_col)) stop("Differential table needs a log-fold-change or signed statistic column for ranked enrichment")
  out <- x; out$gene <- as.character(x[[gene_col]])
  out$log2FoldChange <- if (is.null(lfc_col)) NA_real_ else suppressWarnings(as.numeric(x[[lfc_col]]))
  out$pvalue <- if (is.null(p_col)) NA_real_ else suppressWarnings(as.numeric(x[[p_col]]))
  out$padj <- if (is.null(padj_col)) NA_real_ else suppressWarnings(as.numeric(x[[padj_col]]))
  out$stat <- if (is.null(stat_col)) out$log2FoldChange else suppressWarnings(as.numeric(x[[stat_col]]))
  if (anyDuplicated(out$gene)) warning("Differential table contains duplicate gene identifiers; GSEA will retain the statistic with greatest absolute magnitude")
  thresholds <- list(padj = as.numeric(cfg_get(config, "analysis.padj_threshold", 0.05)), lfc = as.numeric(cfg_get(config, "analysis.lfc_threshold", 0.25)))
  if (is.null(sig_col)) {
    out$significance <- ifelse(!is.na(out$padj) & out$padj <= thresholds$padj & out$log2FoldChange >= thresholds$lfc, "Up",
                               ifelse(!is.na(out$padj) & out$padj <= thresholds$padj & out$log2FoldChange <= -thresholds$lfc, "Down", "NS"))
  } else {
    raw <- tolower(trimws(as.character(x[[sig_col]])))
    out$significance <- ifelse(raw %in% c("up", "upregulated", "up-regulated", "significant up"), "Up",
                               ifelse(raw %in% c("down", "downregulated", "down-regulated", "significant down"), "Down", "NS"))
  }
  out$tested <- if (!is.null(p_col)) !is.na(out$pvalue) else if (!is.null(padj_col)) !is.na(out$padj) else is.finite(out$stat)
  out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
  attr(out, "column_mapping") <- data.frame(standard = c("gene", "log2FoldChange", "pvalue", "padj", "stat", "significance"), source = c(gene_col, lfc_col %||% "", p_col %||% "", padj_col %||% "", stat_col %||% if (is.null(lfc_col)) "" else lfc_col, sig_col %||% "derived"), stringsAsFactors = FALSE)
  out
}

make_table_tasks <- function(config) {
  specs <- cfg_get(config, "input.differential_tables")
  if (is.null(specs)) {
    path <- cfg_get(config, "input.differential_table", cfg_get(config, "enrichment.input_results", required = TRUE))
    specs <- list(list(path = path, id = "input_table"))
  }
  if (!is.list(specs) || !length(specs)) stop("input.differential_tables must be a non-empty array")
  tasks <- list()
  for (i in seq_along(specs)) {
    spec <- specs[[i]]; x <- standardize_de_table(read_de_table(spec$path, spec$sheet), config)
    pop_col <- first_column(x, spec$population_column %||% cfg_get(config, "enrichment.table_columns.population"), c("population", "celltype", "cell_type", "cluster"))
    cmp_col <- first_column(x, spec$comparison_column %||% cfg_get(config, "enrichment.table_columns.comparison"), c("comparison_id", "comparison", "contrast"))
    keys <- data.frame(population = if (is.null(pop_col)) spec$population %||% "all_genes" else as.character(x[[pop_col]]), comparison_id = if (is.null(cmp_col)) spec$comparison_id %||% spec$id %||% paste0("table_", i) else as.character(x[[cmp_col]]), stringsAsFactors = FALSE)
    for (key in split(seq_len(nrow(x)), interaction(keys$population, keys$comparison_id, drop = TRUE))) {
      pop <- keys$population[key[1]]; cmp <- keys$comparison_id[key[1]]; task_id <- safe_name(paste(pop, cmp, sep = "__"))
      tasks[[paste0(task_id, "__", length(tasks) + 1L)]] <- list(id = task_id, population = pop, comparison = list(id = cmp, numerator = spec$numerator %||% "numerator", denominator = spec$denominator %||% "denominator"), result = x[key, , drop = FALSE], mapping = attr(x, "column_mapping"), source = spec$path)
    }
  }
  tasks
}

run_enrichment_only_workflow <- function(config) {
  out <- prepare_output(config); dir.create(file.path(out, "comparisons"), recursive = TRUE, showWarnings = FALSE)
  tasks <- make_table_tasks(config); statuses <- list(); enriched <- list()
  for (i in seq_along(tasks)) {
    task <- tasks[[i]]; task_dir <- file.path(out, "comparisons", task$id); dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(task$mapping, file.path(task_dir, "input_column_mapping.tsv")); write_tsv(task$result, file.path(task_dir, "standardized_input_table.tsv"))
    ans <- tryCatch(run_enrichment(task$result, task_dir, task$population, task$comparison, config), error = function(e) e)
    if (inherits(ans, "error")) {
      writeLines(conditionMessage(ans), file.path(task_dir, "ENRICHMENT_ERROR.txt")); statuses[[i]] <- data.frame(task_id = task$id, population = task$population, comparison_id = task$comparison$id, status = "failed", message = conditionMessage(ans), n_genes = nrow(task$result), source = task$source, stringsAsFactors = FALSE)
    } else {
      if (nrow(ans)) enriched[[task$id]] <- ans
      statuses[[i]] <- data.frame(task_id = task$id, population = task$population, comparison_id = task$comparison$id, status = if (nrow(ans)) "completed" else "empty", message = "", n_genes = nrow(task$result), source = task$source, stringsAsFactors = FALSE)
    }
  }
  status <- do.call(rbind, statuses); write_tsv(status, file.path(out, "task_status.tsv"))
  artifacts <- file.path(out, "task_status.tsv")
  if (length(enriched)) { combined <- rbind_fill(enriched); write_tsv(combined, file.path(out, "enrichment_all_comparisons.tsv")); artifacts <- c(artifacts, file.path(out, "enrichment_all_comparisons.tsv")) }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt")); artifacts <- c(artifacts, file.path(out, "sessionInfo.txt"))
  write_run_manifest(config, "10-scrna-run-differential-analysis", out, artifacts, c("stage=enrichment_only", "Differential input tables were not modified"))
  if (!any(status$status %in% c("completed", "empty"))) stop("No table enrichment task completed; inspect task_status.tsv")
}

normalize_comparisons <- function(config) {
  xs <- cfg_get(config, "comparisons")
  if (is.null(xs)) xs <- list(cfg_get(config, "comparison", required = TRUE))
  if (!is.list(xs) || !length(xs)) stop("comparisons must be a non-empty JSON array")
  lapply(seq_along(xs), function(i) {
    x <- xs[[i]]; numerator <- x$numerator; denominator <- x$denominator
    if (is.null(numerator) || is.null(denominator) || numerator == denominator) stop("Each comparison needs different numerator and denominator")
    list(id = x$id %||% safe_name(paste(numerator, "vs", denominator)), numerator = as.character(numerator), denominator = as.character(denominator))
  })
}

validate_sample_mapping <- function(meta, sample_col, condition_col, covariates) {
  design <- unique(meta[c(sample_col, condition_col, covariates)])
  n_condition <- tapply(as.character(design[[condition_col]]), design[[sample_col]], function(x) length(unique(x)))
  if (any(n_condition > 1L)) stop("A sample maps to multiple conditions: ", paste(names(n_condition)[n_condition > 1L], collapse = ", "))
  invisible(TRUE)
}

join_assay_layers <- function(obj, assay) {
  if (requireNamespace("SeuratObject", quietly = TRUE) && utils::packageVersion("SeuratObject") >= "5.0.0") {
    layers <- SeuratObject::Layers(obj[[assay]])
    if (sum(grepl("^counts(\\.|$)", layers)) > 1L || sum(grepl("^data(\\.|$)", layers)) > 1L) obj <- Seurat::JoinLayers(obj, assay = assay)
  }
  obj
}

select_populations <- function(meta, population_col, spec) {
  if (is.null(population_col) || !nzchar(population_col)) return("all_cells")
  available <- sort(unique(as.character(meta[[population_col]][!is.na(meta[[population_col]])])))
  mode <- spec$mode %||% "all"
  include <- as_chr(spec$include); exclude <- as_chr(spec$exclude)
  selected <- if (mode == "all") available else intersect(include, available)
  selected <- setdiff(selected, exclude)
  if (!length(selected)) stop("Population selection matched zero metadata values")
  selected
}

population_cells <- function(meta, population_col, population) {
  if (population == "all_cells" || is.null(population_col)) return(rownames(meta))
  rownames(meta)[!is.na(meta[[population_col]]) & as.character(meta[[population_col]]) == population]
}

make_design_audit <- function(meta, sample_col, condition_col, population_col, populations, comparisons) {
  rows <- lapply(populations, function(pop) {
    x <- meta[population_cells(meta, population_col, pop), , drop = FALSE]
    z <- as.data.frame(table(sample = x[[sample_col]], condition = x[[condition_col]]), stringsAsFactors = FALSE)
    z <- z[z$Freq > 0, , drop = FALSE]; names(z)[3] <- "n_cells"; z$population <- pop; z
  })
  out <- do.call(rbind, rows); out[, c("population", "sample", "condition", "n_cells")]
}

sample_population_audit <- function(meta, sample_col, condition_col, min_cells) {
  x <- as.data.frame(table(sample = meta[[sample_col]], condition = meta[[condition_col]]), stringsAsFactors = FALSE)
  names(x)[3] <- "n_cells"; x <- x[x$n_cells > 0, , drop = FALSE]; x$passes_min_cells <- x$n_cells >= min_cells; x
}

build_coldata <- function(meta, sample_col, condition_col, covariates, samples, denominator) {
  fields <- unique(c(sample_col, condition_col, covariates)); x <- unique(meta[fields])
  if (anyDuplicated(x[[sample_col]])) stop("Sample-level covariates are not constant within sample")
  x <- x[match(samples, x[[sample_col]]), , drop = FALSE]; rownames(x) <- x[[sample_col]]
  x[[condition_col]] <- stats::relevel(factor(x[[condition_col]]), denominator); x
}

run_pseudobulk <- function(obj, meta, assay, sample_col, condition_col, covariates, comparison, thresholds, config, task_dir) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Package 'DESeq2' is required for pseudobulk")
  counts <- get_raw_counts(obj, assay)
  groups <- factor(meta[colnames(counts), sample_col])
  mm <- Matrix::sparse.model.matrix(~ 0 + groups)
  pb <- counts %*% mm; colnames(pb) <- sub("^groups", "", colnames(pb))
  keep <- Matrix::rowSums(pb) >= thresholds$min_total_count & Matrix::rowSums(pb >= thresholds$min_count_per_sample) >= thresholds$min_samples_expressed
  pb <- pb[keep, , drop = FALSE]
  if (!nrow(pb)) stop("No genes passed minimum total count")
  coldata <- build_coldata(meta, sample_col, condition_col, covariates, colnames(pb), comparison$denominator)
  design_text <- cfg_get(config, "analysis.design")
  if (is.null(design_text)) design_text <- paste("~", paste(c(covariates, condition_col), collapse = " + "))
  design <- stats::as.formula(design_text); mm_design <- stats::model.matrix(design, coldata)
  if (qr(mm_design)$rank < ncol(mm_design)) stop("Design matrix is rank deficient")
  dds <- DESeq2::DESeqDataSetFromMatrix(as.matrix(pb), coldata, design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c(condition_col, comparison$numerator, comparison$denominator), independentFiltering = TRUE)
  if (isTRUE(cfg_get(config, "analysis.lfc_shrink", TRUE))) {
    if (!requireNamespace("apeglm", quietly = TRUE)) warning("lfc_shrink was requested but apeglm is unavailable; retaining unshrunk DESeq2 effect sizes") else {
      coef_name <- grep(paste0("^", condition_col, "_", comparison$numerator, "_vs_", comparison$denominator, "$"), DESeq2::resultsNames(dds), value = TRUE)
      if (length(coef_name) == 1L) res <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm") else warning("Could not identify one DESeq2 coefficient for apeglm shrinkage; retaining unshrunk effect sizes")
    }
  }
  out <- as.data.frame(res); out$gene <- rownames(out)
  norm <- DESeq2::counts(dds, normalized = TRUE)
  gr <- coldata[[condition_col]]
  out$mean_numerator <- rowMeans(norm[, gr == comparison$numerator, drop = FALSE])
  out$mean_denominator <- rowMeans(norm[, gr == comparison$denominator, drop = FALSE])
  write_tsv(data.frame(sample = colnames(pb), condition = as.character(gr), coldata, check.names = FALSE), file.path(task_dir, "sample_design.tsv"))
  saveRDS(list(counts = pb, normalized_counts = norm, coldata = coldata, design = design_text), file.path(task_dir, "pseudobulk_data.rds"))
  plot_pseudobulk(norm, coldata, condition_col, out, task_dir, config)
  out
}

run_cell_level <- function(obj, assay, condition_col, comparison, method, config) {
  Seurat::Idents(obj) <- factor(obj[[]][[condition_col]])
  test <- c(seurat_wilcox = "wilcox", seurat_mast = "MAST", seurat_lr = "LR")[[method]]
  x <- Seurat::FindMarkers(obj, assay = assay, ident.1 = comparison$numerator, ident.2 = comparison$denominator, test.use = test,
                           logfc.threshold = as.numeric(cfg_get(config, "analysis.test_logfc_threshold", 0)),
                           min.pct = as.numeric(cfg_get(config, "analysis.min_pct", 0.01)), return.thresh = 1, verbose = FALSE)
  x$gene <- rownames(x); fc <- intersect(c("avg_log2FC", "avg_logFC", "avg_diff"), names(x))[1]
  names(x)[names(x) == fc] <- "log2FoldChange"; names(x)[names(x) == "p_val"] <- "pvalue"; names(x)[names(x) == "p_val_adj"] <- "padj"; x
}

annotate_de_result <- function(x, population, comparison, method, thresholds, meta, sample_col, condition_col) {
  if (!"log2FoldChange" %in% names(x)) stop("Differential result has no recognized log2 fold-change column")
  if (!"pvalue" %in% names(x)) x$pvalue <- NA_real_; if (!"padj" %in% names(x)) x$padj <- NA_real_
  x$population <- population; x$comparison_id <- comparison$id; x$numerator <- comparison$numerator; x$denominator <- comparison$denominator
  x$direction <- ifelse(is.na(x$log2FoldChange), "Not_tested", ifelse(x$log2FoldChange > 0, "Up", ifelse(x$log2FoldChange < 0, "Down", "Stable")))
  x$significance <- ifelse(is.na(x$padj), "Not_tested", ifelse(x$padj <= thresholds$padj & x$log2FoldChange >= thresholds$lfc, "Up", ifelse(x$padj <= thresholds$padj & x$log2FoldChange <= -thresholds$lfc, "Down", "NS")))
  x$tested <- !is.na(x$pvalue); x$filter_reason <- ifelse(x$tested, "", "independent_filtering_or_unavailable")
  x$method <- method; x$inference_level <- if (method == "pseudobulk_deseq2") "sample_level_formal" else "cell_level_exploratory"
  tab <- table(meta[[condition_col]]); smp <- unique(meta[c(sample_col, condition_col)]); stab <- table(smp[[condition_col]])
  x$n_cells_numerator <- unname(tab[comparison$numerator]); x$n_cells_denominator <- unname(tab[comparison$denominator])
  x$n_samples_numerator <- unname(stab[comparison$numerator]); x$n_samples_denominator <- unname(stab[comparison$denominator])
  front <- c("gene", "population", "comparison_id", "numerator", "denominator", "log2FoldChange", "pvalue", "padj", "direction", "significance", "tested", "filter_reason", "method", "inference_level")
  x[, c(front, setdiff(names(x), front)), drop = FALSE]
}

write_de_tables <- function(x, out) {
  write_tsv(x, file.path(out, "all_genes.tsv")); sig <- x[x$significance %in% c("Up", "Down"), , drop = FALSE]
  write_tsv(sig, file.path(out, "significant_genes.tsv")); write_tsv(x[x$significance == "Up", , drop = FALSE], file.path(out, "upregulated_genes.tsv")); write_tsv(x[x$significance == "Down", , drop = FALSE], file.path(out, "downregulated_genes.tsv"))
}

plot_de_results <- function(x, out, population, comparison, config) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  x$plot_p <- -log10(pmax(x$padj, .Machine$double.xmin)); x$plot_p[!is.finite(x$plot_p)] <- 0
  p <- ggplot2::ggplot(x, ggplot2::aes(log2FoldChange, plot_p, color = significance)) + ggplot2::geom_point(alpha = .65, size = .9) +
    ggplot2::scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70", Not_tested = "grey90"), drop = FALSE) +
    ggplot2::labs(title = paste(population, comparison$id), x = paste0("log2FC (", comparison$numerator, " vs ", comparison$denominator, ")"), y = "-log10 adjusted P") + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(out, "volcano.pdf"), p, width = 7, height = 6)
  if ("baseMean" %in% names(x)) {
    ma <- ggplot2::ggplot(x, ggplot2::aes(log10(baseMean + 1), log2FoldChange, color = significance)) + ggplot2::geom_point(alpha = .6, size = .8) + ggplot2::scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70", Not_tested = "grey90"), drop = FALSE) + ggplot2::theme_bw() + ggplot2::labs(x = "log10(baseMean + 1)", y = "log2FC")
    ggplot2::ggsave(file.path(out, "MA_plot.pdf"), ma, width = 7, height = 6)
  }
}

plot_pseudobulk <- function(norm, coldata, condition_col, result, out, config) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  vst <- log2(norm + 1); pc <- stats::prcomp(t(vst), scale. = FALSE); pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  d <- data.frame(sample = rownames(pc$x), PC1 = pc$x[,1], PC2 = pc$x[,2], condition = coldata[rownames(pc$x), condition_col])
  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = condition, label = sample)) + ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.7, size = 3) + ggplot2::theme_bw() + ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
  ggplot2::ggsave(file.path(out, "pseudobulk_PCA.pdf"), p, width = 7, height = 6)
  top_n <- as.integer(cfg_get(config, "plots.top_genes", 30)); ord <- order(result$padj, -abs(result$log2FoldChange), na.last = NA); genes <- head(result$gene[ord], top_n)
  if (length(genes) >= 2L) {
    z <- t(scale(t(vst[genes, , drop = FALSE]))); grDevices::pdf(file.path(out, "top_DE_heatmap.pdf"), width = 8, height = max(5, length(genes) * .16 + 2)); stats::heatmap(z, Colv = NA, scale = "none", margins = c(8, 8)); grDevices::dev.off()
  }
}

plot_batch_summary <- function(x, status, out) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  sig <- x[x$significance %in% c("Up", "Down"), , drop = FALSE]
  if (!nrow(sig)) return(invisible(NULL))
  z <- aggregate(gene ~ population + comparison_id + significance, sig, length); names(z)[4] <- "n_genes"
  if (nrow(z)) { p <- ggplot2::ggplot(z, ggplot2::aes(population, n_genes, fill = significance)) + ggplot2::geom_col(position = "dodge") + ggplot2::facet_wrap(~comparison_id, scales = "free_x") + ggplot2::scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4")) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)); ggplot2::ggsave(file.path(out, "DEG_count_summary.pdf"), p, width = max(8, length(unique(z$population)) * .45 + 4), height = 6) }
}

run_enrichment <- function(result, out, population, comparison, config) {
  enr_dir <- file.path(out, "enrichment"); dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Package 'clusterProfiler' is required when enrichment is enabled")
  species <- tolower(cfg_get(config, "enrichment.species", required = TRUE)); id_type <- cfg_get(config, "enrichment.gene_id_type", "SYMBOL")
  human <- species %in% c("human", "homo_sapiens", "homo sapiens")
  mouse <- species %in% c("mouse", "mus_musculus", "mus musculus")
  if (!human && !mouse) stop("Built-in comprehensive enrichment currently supports human or mouse")
  org_pkg <- if (human) "org.Hs.eg.db" else "org.Mm.eg.db"
  msig_species <- if (human) "Homo sapiens" else "Mus musculus"
  if (!requireNamespace(org_pkg, quietly = TRUE)) stop("Required annotation package is unavailable: ", org_pkg)
  orgdb <- get(org_pkg, envir = asNamespace(org_pkg)); universe <- result$gene[result$tested]
  mapping <- clusterProfiler::bitr(unique(result$gene), fromType = id_type, toType = "ENTREZID", OrgDb = orgdb)
  mapping$ambiguous_input_id <- duplicated(mapping[[id_type]]) | duplicated(mapping[[id_type]], fromLast = TRUE)
  write_tsv(mapping, file.path(enr_dir, "gene_id_mapping.tsv")); mapping_one <- mapping[!duplicated(mapping[[id_type]]), , drop = FALSE]
  map <- setNames(mapping_one$ENTREZID, mapping_one[[id_type]])
  uni <- unique(stats::na.omit(unname(map[universe]))); min_input <- as.integer(cfg_get(config, "enrichment.min_input_genes", 5))
  requested <- toupper(as_chr(cfg_get(config, "enrichment.databases", c("GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME", "HALLMARK"))))
  min_gs <- as.integer(cfg_get(config, "enrichment.min_gene_set_size", 10)); max_gs <- as.integer(cfg_get(config, "enrichment.max_gene_set_size", 500))
  rows <- list(); statuses <- list(); k <- 0L
  add_status <- function(database, method, direction, status, n_input, n_terms, message = "") {
    k <<- k + 1L; statuses[[k]] <<- data.frame(database, method, direction, status, n_input, n_terms, message, stringsAsFactors = FALSE)
  }
  add_result <- function(tab, database, method, direction) {
    if (!nrow(tab)) return(invisible(NULL)); tab$database <- database; tab$method <- method; tab$direction <- direction
    tab$analysis <- paste(method, database, sep = "_"); rows[[paste(database, method, direction, length(rows), sep = "_")]] <<- tab
    write_tsv(tab, file.path(enr_dir, paste0(tolower(method), "_", tolower(database), "_", tolower(direction), "_full.tsv")))
  }
  rank0 <- result$stat; if (is.null(rank0)) rank0 <- result$log2FoldChange; names(rank0) <- result$gene
  rank0 <- rank0[is.finite(rank0) & names(rank0) %in% names(map)]
  collapsed <- tapply(rank0, unname(map[names(rank0)]), function(v) v[which.max(abs(v))])
  rank <- as.numeric(collapsed); names(rank) <- names(collapsed); rank <- sort(rank, decreasing = TRUE)

  run_one <- function(database, ont = NULL, term2gene = NULL) {
    for (direction in c("Up", "Down")) {
      genes <- unique(stats::na.omit(unname(map[result$gene[result$significance == direction]])))
      if (length(genes) < min_input) { add_status(database, "ORA", direction, "skipped_too_few_genes", length(genes), 0L); next }
      ans <- tryCatch({
        x <- if (!is.null(ont)) clusterProfiler::enrichGO(genes, OrgDb = orgdb, keyType = "ENTREZID", ont = ont, universe = uni, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE) else clusterProfiler::enricher(genes, universe = uni, TERM2GENE = term2gene, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs)
        tab <- as.data.frame(x); add_result(tab, database, "ORA", direction); add_status(database, "ORA", direction, if (nrow(tab)) "completed" else "empty", length(genes), nrow(tab))
      }, error = function(e) add_status(database, "ORA", direction, "failed", length(genes), 0L, conditionMessage(e)))
    }
    if (length(rank) < 10L) return(add_status(database, "GSEA", "ranked", "skipped_too_few_genes", length(rank), 0L))
    tryCatch({
      x <- if (!is.null(ont)) clusterProfiler::gseGO(rank, OrgDb = orgdb, keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs, eps = 0, verbose = FALSE) else clusterProfiler::GSEA(rank, TERM2GENE = term2gene, pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs, eps = 0, verbose = FALSE)
      tab <- as.data.frame(x); if (nrow(tab)) { tab$enrichment_direction <- ifelse(tab$NES >= 0, "Up", "Down"); add_result(tab, database, "GSEA", "ranked") }
      add_status(database, "GSEA", "ranked", if (nrow(tab)) "completed" else "empty", length(rank), nrow(tab))
    }, error = function(e) add_status(database, "GSEA", "ranked", "failed", length(rank), 0L, conditionMessage(e)))
  }

  go_defs <- list(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
  for (database in intersect(names(go_defs), requested)) run_one(database, ont = go_defs[[database]])
  msig_defs <- list(KEGG = c("C2", "CP:KEGG"), REACTOME = c("C2", "CP:REACTOME"), HALLMARK = c("H", NA_character_))
  if (length(intersect(names(msig_defs), requested)) && !requireNamespace("msigdbr", quietly = TRUE)) {
    for (database in intersect(names(msig_defs), requested)) add_status(database, "ALL", "all", "missing_dependency", 0L, 0L, "Package 'msigdbr' is unavailable")
  } else for (database in intersect(names(msig_defs), requested)) {
    def <- msig_defs[[database]]
    msig <- tryCatch(msigdbr::msigdbr(species = msig_species, category = def[1], subcategory = if (is.na(def[2])) NULL else def[2]), error = function(e) e)
    if (inherits(msig, "error")) { add_status(database, "ALL", "all", "failed", 0L, 0L, conditionMessage(msig)); next }
    term2gene <- unique(msig[c("gs_name", "entrez_gene")]); names(term2gene) <- c("term", "gene"); term2gene$gene <- as.character(term2gene$gene)
    if (!nrow(term2gene)) add_status(database, "ALL", "all", "empty_gene_set", 0L, 0L) else run_one(database, term2gene = term2gene)
  }
  status <- if (length(statuses)) do.call(rbind, statuses) else data.frame(); write_tsv(status, file.path(enr_dir, "enrichment_status.tsv"))
  if (!length(rows)) return(data.frame())
  ans <- rbind_fill(rows); ans$population <- population; ans$comparison_id <- comparison$id
  plot_enrichment_summary(ans, enr_dir, config); ans
}

plot_enrichment_summary <- function(x, out, config) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(x)) return(invisible(NULL))
  padj_col <- intersect(c("p.adjust", "pvalue"), names(x))[1]; if (is.na(padj_col)) return(invisible(NULL))
  x$plot_score <- -log10(pmax(x[[padj_col]], .Machine$double.xmin)); x$label <- x$Description
  split_key <- interaction(x$database, x$method, x$direction, drop = TRUE)
  top_n <- as.integer(cfg_get(config, "enrichment.plot_top_terms", 10)); keep <- unlist(lapply(split(seq_len(nrow(x)), split_key), function(i) head(i[order(x[[padj_col]][i])], top_n)))
  d <- x[unique(keep), , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(plot_score, stats::reorder(label, plot_score), color = database)) + ggplot2::geom_point(size = 2) + ggplot2::facet_grid(method + direction ~ database, scales = "free_y", space = "free_y") + ggplot2::theme_bw() + ggplot2::labs(x = "-log10 adjusted P", y = NULL)
  ggplot2::ggsave(file.path(out, "enrichment_dotplot.pdf"), p, width = 14, height = max(8, nrow(d) * .12 + 3), limitsize = FALSE)
  gsea <- d[d$method == "GSEA" & "NES" %in% names(d), , drop = FALSE]
  if (nrow(gsea)) { q <- ggplot2::ggplot(gsea, ggplot2::aes(NES, stats::reorder(label, NES), color = database)) + ggplot2::geom_point(size = 2) + ggplot2::facet_wrap(~database, scales = "free_y") + ggplot2::geom_vline(xintercept = 0, linetype = 2) + ggplot2::theme_bw() + ggplot2::labs(y = NULL); ggplot2::ggsave(file.path(out, "gsea_NES_dotplot.pdf"), q, width = 12, height = max(7, nrow(gsea) * .12 + 3), limitsize = FALSE) }
}

classify_failure <- function(x) {
  if (grepl("Insufficient independent samples", x)) "skipped_low_replicates" else if (grepl("rank deficient|multiple conditions|not constant", x)) "invalid_design" else if (grepl("required|unavailable", x, ignore.case = TRUE)) "missing_dependency" else "failed"
}
