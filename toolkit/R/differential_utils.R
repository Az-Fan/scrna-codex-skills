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

initialize_enrichment_seed <- function(config) {
  seed <- suppressWarnings(as.integer(cfg_get(config, "random_seed", 1L)))
  if (length(seed) != 1L || is.na(seed)) stop("random_seed must be one integer")
  attr(config, "resolved_random_seed") <- seed
  set.seed(seed)
  config
}

run_enrichment_only_workflow <- function(config) {
  config <- initialize_enrichment_seed(config)
  out <- prepare_output(config); dir.create(file.path(out, "comparisons"), recursive = TRUE, showWarnings = FALSE)
  tasks <- make_table_tasks(config); statuses <- list(); enriched <- list()
  for (i in seq_along(tasks)) {
    task <- tasks[[i]]; task_dir <- file.path(out, "comparisons", task$id); dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
    stale_error <- file.path(task_dir, "ENRICHMENT_ERROR.txt")
    if (file.exists(stale_error)) unlink(stale_error)
    write_tsv(task$mapping, file.path(task_dir, "input_column_mapping.tsv")); write_tsv(task$result, file.path(task_dir, "standardized_input_table.tsv"))
    ans <- tryCatch(run_enrichment(task$result, task_dir, task$population, task$comparison, config), error = function(e) e)
    if (inherits(ans, "error")) {
      writeLines(conditionMessage(ans), file.path(task_dir, "ENRICHMENT_ERROR.txt")); statuses[[i]] <- data.frame(task_id = task$id, population = task$population, comparison_id = task$comparison$id, status = "failed", message = conditionMessage(ans), n_genes = nrow(task$result), source = task$source, stringsAsFactors = FALSE)
    } else {
      if (nrow(ans)) enriched[[task$id]] <- ans
      enrichment_state <- summarize_enrichment_status(task_dir)
      statuses[[i]] <- data.frame(task_id = task$id, population = task$population, comparison_id = task$comparison$id, status = enrichment_state, message = "", n_genes = nrow(task$result), source = task$source, stringsAsFactors = FALSE)
    }
  }
  status <- do.call(rbind, statuses); write_tsv(status, file.path(out, "task_status.tsv"))
  artifacts <- file.path(out, "task_status.tsv")
  if (length(enriched)) { combined <- rbind_fill(enriched); write_tsv(combined, file.path(out, "enrichment_all_comparisons.tsv")); artifacts <- c(artifacts, file.path(out, "enrichment_all_comparisons.tsv")) }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt")); artifacts <- c(artifacts, file.path(out, "sessionInfo.txt"))
  comparison_artifacts <- list.files(file.path(out, "comparisons"), recursive = TRUE, full.names = TRUE)
  artifacts <- c(artifacts, comparison_artifacts[file.info(comparison_artifacts)$isdir %in% FALSE])
  active_skill <- Sys.getenv("SCRNA_ACTIVE_SKILL", unset = "11-scrna-run-differential-analysis")
  write_run_manifest(config, active_skill, out, artifacts, c("stage=enrichment_only", "Differential input tables were not modified"))
  if (!any(status$status %in% c("completed", "partial", "empty"))) stop("No table enrichment task completed; inspect task_status.tsv")
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
  raw_res <- DESeq2::results(dds, contrast = c(condition_col, comparison$numerator, comparison$denominator), independentFiltering = TRUE)
  res <- raw_res
  shrink_requested <- isTRUE(cfg_get(config, "analysis.lfc_shrink", TRUE))
  allow_unshrunk <- isTRUE(cfg_get(config, "analysis.allow_unshrunk_lfc", FALSE))
  shrink_applied <- FALSE
  shrink_status <- if (shrink_requested) "requested" else "not_requested"
  shrink_reason <- ""
  if (shrink_requested) {
    if (!requireNamespace("apeglm", quietly = TRUE)) {
      shrink_status <- "missing_dependency"
      shrink_reason <- "Package 'apeglm' is unavailable"
      if (!allow_unshrunk) {
        stop("lfc_shrink=true requires package 'apeglm'; install it or explicitly set analysis.allow_unshrunk_lfc=true")
      }
      warning(shrink_reason, "; explicit analysis.allow_unshrunk_lfc=true permits unshrunk DESeq2 effect sizes")
    } else {
      coef_name <- grep(paste0("^", condition_col, "_", comparison$numerator, "_vs_", comparison$denominator, "$"), DESeq2::resultsNames(dds), value = TRUE)
      if (length(coef_name) != 1L) {
        shrink_status <- "coefficient_not_identified"
        shrink_reason <- "Could not identify exactly one DESeq2 coefficient for apeglm shrinkage"
        if (!allow_unshrunk) stop(shrink_reason, "; set analysis.allow_unshrunk_lfc=true only if the fallback is intentional")
        warning(shrink_reason, "; explicit analysis.allow_unshrunk_lfc=true permits unshrunk DESeq2 effect sizes")
      } else {
        shrunken <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
        res$log2FoldChange <- shrunken$log2FoldChange
        res$lfcSE <- shrunken$lfcSE
        shrink_applied <- TRUE
        shrink_status <- "applied"
      }
    }
  }
  shrink_audit <- data.frame(
    requested = shrink_requested,
    applied = shrink_applied,
    method = if (shrink_applied) "apeglm" else "none",
    inferential_statistics_source = "unshrunk_DESeq2_Wald",
    status = shrink_status,
    reason = shrink_reason,
    allow_unshrunk_lfc = allow_unshrunk,
    stringsAsFactors = FALSE
  )
  write_tsv(shrink_audit, file.path(task_dir, "effect_size_audit.tsv"))
  out <- as.data.frame(res); out$gene <- rownames(out)
  out$lfc_shrink_requested <- shrink_requested
  out$lfc_shrink_applied <- shrink_applied
  out$lfc_shrink_method <- if (shrink_applied) "apeglm" else "none"
  out$lfc_shrink_status <- shrink_status
  out$inferential_statistics_source <- "unshrunk_DESeq2_Wald"
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
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  }
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
  msig_defs <- list(
    KEGG = list(collection = "C2", subcollections = c("CP:KEGG_MEDICUS", "CP:KEGG_LEGACY", "CP:KEGG")),
    REACTOME = list(collection = "C2", subcollections = "CP:REACTOME"),
    HALLMARK = list(collection = "H", subcollections = character())
  )
  if (length(intersect(names(msig_defs), requested)) && !requireNamespace("msigdbr", quietly = TRUE)) {
    for (database in intersect(names(msig_defs), requested)) add_status(database, "ALL", "all", "missing_dependency", 0L, 0L, "Package 'msigdbr' is unavailable")
  } else for (database in intersect(names(msig_defs), requested)) {
    def <- msig_defs[[database]]
    available <- tryCatch(msigdbr::msigdbr_collections(), error = function(e) NULL)
    subcollections <- def$subcollections
    if (length(subcollections) && !is.null(available)) {
      subcollections <- intersect(subcollections, available$gs_subcollection[available$gs_collection == def$collection])
    }
    load_one <- function(subcollection = NULL) {
      args <- list(species = msig_species, collection = def$collection)
      if (!is.null(subcollection) && nzchar(subcollection)) args$subcollection <- subcollection
      do.call(msigdbr::msigdbr, args)
    }
    msig <- tryCatch({
      if (length(def$subcollections) && !length(subcollections)) {
        stop("No supported MSigDB subcollection is available for ", database)
      }
      if (length(subcollections)) do.call(rbind, lapply(subcollections, load_one)) else load_one()
    }, error = function(e) e)
    if (inherits(msig, "error")) { add_status(database, "ALL", "all", "failed", 0L, 0L, conditionMessage(msig)); next }
    gene_column <- intersect(c("ncbi_gene", "entrez_gene"), names(msig))[1]
    if (is.na(gene_column)) {
      add_status(database, "ALL", "all", "failed", 0L, 0L, "MSigDB result lacks an NCBI/Entrez gene column")
      next
    }
    term2gene <- unique(msig[c("gs_name", gene_column)]); names(term2gene) <- c("term", "gene"); term2gene$gene <- as.character(term2gene$gene)
    if (!nrow(term2gene)) add_status(database, "ALL", "all", "empty_gene_set", 0L, 0L) else run_one(database, term2gene = term2gene)
  }
  status <- if (length(statuses)) do.call(rbind, statuses) else data.frame(); write_tsv(status, file.path(enr_dir, "enrichment_status.tsv"))
  if (!length(rows)) return(data.frame())
  ans <- rbind_fill(rows); ans$population <- population; ans$comparison_id <- comparison$id
  plot_enrichment_summary(ans, enr_dir, config); ans
}

summarize_enrichment_status <- function(task_dir) {
  path <- file.path(task_dir, "enrichment", "enrichment_status.tsv")
  if (!file.exists(path)) return("failed")
  status <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(status) || !"status" %in% names(status)) return("empty")
  bad <- status$status %in% c("failed", "missing_dependency", "empty_gene_set")
  usable <- status$status %in% c("completed", "empty", "skipped_too_few_genes")
  if (any(bad) && any(usable)) "partial" else if (any(bad)) "failed" else if (any(status$status == "completed")) "completed" else "empty"
}

clean_enrichment_label <- function(value, width = 55) {
  value <- as.character(value)
  value <- gsub("^(HALLMARK|REACTOME)_", "", value)
  value <- gsub("^KEGG_(MEDICUS_REFERENCE_)?", "", value)
  value <- gsub("_+", " ", value)
  value <- trimws(gsub("[[:space:]]+", " ", value))
  value <- vapply(value, function(item) {
    if (grepl("[A-Z]", item) && !grepl("[a-z]", item)) {
      item <- tolower(item)
      item <- paste0(toupper(substr(item, 1, 1)), substr(item, 2, nchar(item)))
    }
    acronyms <- c(
      dna = "DNA", rna = "RNA", mrna = "mRNA", trna = "tRNA", rrna = "rRNA",
      mirna = "miRNA", snrnp = "snRNP", mhc = "MHC", gpcr = "GPCR", rqc = "RQC",
      ros = "ROS", ecm = "ECM", tnf = "TNF", tnfa = "TNFA", nfkb = "NFKB",
      tgfb = "TGFB", vegf = "VEGF", egfr = "EGFR", pi3k = "PI3K", akt = "AKT",
      mtorc1 = "MTORC1", p53 = "P53", cd28 = "CD28", ccr5 = "CCR5", gnb = "GNB",
      plcb = "PLCB", pkc = "PKC", eif = "EIF", eifs = "EIFS", eif2ak4 = "EIF2AK4",
      gcn2 = "GCN2", htt = "HTT", creb = "CREB", atr = "ATR", prc2 = "PRC2",
      cams = "CAMs", cxcr4 = "CXCR4", itpr = "ITPR", g = "G",
      i = "I", ii = "II", iii = "III", iv = "IV"
    )
    for (key in names(acronyms)) {
      item <- gsub(paste0("\\b", key, "\\b"), acronyms[[key]], item,
                   ignore.case = TRUE, perl = TRUE)
    }
    item
  }, character(1), USE.NAMES = FALSE)
  vapply(value, function(item) paste(strwrap(item, width = width), collapse = "\n"), character(1))
}

enrichment_term_members <- function(x, i) {
  column <- if (x$method[[i]] == "GSEA" && "core_enrichment" %in% names(x)) "core_enrichment" else "geneID"
  if (!column %in% names(x) || is.na(x[[column]][[i]]) || !nzchar(x[[column]][[i]])) return(character())
  unique(strsplit(as.character(x[[column]][[i]]), "[/;,]")[[1]])
}

enrichment_gene_overlap <- function(a, b) {
  if (!length(a) || !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

select_enrichment_plot_terms <- function(x, top_n, ora_fdr = 0.05, gsea_fdr = 0.25,
                                         max_gene_overlap = 0.6) {
  padj_col <- intersect(c("p.adjust", "pvalue"), names(x))[1]
  if (!nrow(x) || is.na(padj_col)) return(x[FALSE, , drop = FALSE])
  x$plot_fdr <- suppressWarnings(as.numeric(x[[padj_col]]))
  x$plot_direction <- as.character(x$direction)
  is_gsea <- x$method == "GSEA"
  if ("enrichment_direction" %in% names(x)) {
    use <- is_gsea & !is.na(x$enrichment_direction) & nzchar(x$enrichment_direction)
    x$plot_direction[use] <- as.character(x$enrichment_direction[use])
  }
  x$plot_threshold <- ifelse(is_gsea, gsea_fdr, ora_fdr)
  x <- x[is.finite(x$plot_fdr) & x$plot_fdr <= x$plot_threshold & x$plot_direction %in% c("Up", "Down"), , drop = FALSE]
  if (!nrow(x)) return(x)
  split_key <- interaction(x$database, x$method, x$plot_direction, drop = TRUE)
  keep <- unlist(lapply(split(seq_len(nrow(x)), split_key), function(i) {
    secondary <- if (all(x$method[i] == "GSEA") && "NES" %in% names(x)) -abs(x$NES[i]) else seq_along(i)
    candidates <- i[order(x$plot_fdr[i], secondary, na.last = TRUE)]
    selected <- integer()
    for (candidate in candidates) {
      members <- enrichment_term_members(x, candidate)
      duplicate_label <- length(selected) && x$Description[[candidate]] %in% x$Description[selected]
      redundant <- length(selected) && length(members) && any(vapply(selected, function(previous) {
        enrichment_gene_overlap(members, enrichment_term_members(x, previous)) > max_gene_overlap
      }, logical(1)))
      if (!duplicate_label && !redundant) selected <- c(selected, candidate)
      if (length(selected) >= top_n) break
    }
    selected
  }), use.names = FALSE)
  selected <- x[unique(keep), , drop = FALSE]
  selected$plot_rank <- ave(seq_len(nrow(selected)),
                            interaction(selected$database, selected$method, selected$plot_direction, drop = TRUE),
                            FUN = seq_along)
  selected
}

enrichment_plot_context <- function(x, config) {
  population <- unique(stats::na.omit(as.character(x$population)))
  comparison <- unique(stats::na.omit(as.character(x$comparison_id)))
  list(
    population = cfg_get(config, "enrichment.plot_population_label",
                         if (length(population)) gsub("_+", " ", population[[1]]) else ""),
    comparison = cfg_get(config, "enrichment.plot_comparison_label",
                         if (length(comparison)) gsub("_+", " ", comparison[[1]]) else ""),
    positive = cfg_get(config, "enrichment.plot_positive_label", "Positive NES"),
    negative = cfg_get(config, "enrichment.plot_negative_label", "Negative NES")
  )
}

enrichment_detail_height <- function(n_terms) {
  min(10, max(3.5, 0.42 * as.integer(n_terms) + 2.6))
}

format_enrichment_fdr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.finite(x), "",
         ifelse(x < 0.001, "q<0.001", paste0("q=", formatC(x, format = "f", digits = 3))))
}

map_enrichment_plot_values <- function(x) {
  if (!nrow(x)) return(x)
  is_gsea <- x$method == "GSEA"
  x$plot_fdr_label <- format_enrichment_fdr(x$plot_fdr)
  x$plot_x_value <- NA_real_
  x$plot_size_value <- NA_real_
  x$plot_x_metric <- ifelse(is_gsea, "NES", if ("RichFactor" %in% names(x)) "signed RichFactor" else "signed -log10 FDR")
  x$plot_size_metric <- ifelse(is_gsea, "setSize", "Count")
  if (any(is_gsea)) {
    x$plot_x_value[is_gsea] <- suppressWarnings(as.numeric(x$NES[is_gsea]))
    x$plot_size_value[is_gsea] <- suppressWarnings(as.numeric(x$setSize[is_gsea]))
  }
  is_ora <- !is_gsea
  if (any(is_ora)) {
    magnitude <- if ("RichFactor" %in% names(x)) {
      suppressWarnings(as.numeric(x$RichFactor[is_ora]))
    } else {
      -log10(pmax(x$plot_fdr[is_ora], .Machine$double.xmin))
    }
    x$plot_x_value[is_ora] <- ifelse(x$plot_direction[is_ora] == "Up", magnitude, -magnitude)
    x$plot_size_value[is_ora] <- suppressWarnings(as.numeric(x$Count[is_ora]))
  }
  x
}

enrichment_plot_formats <- function(config) {
  requested <- tolower(as.character(unlist(cfg_get(config, "enrichment.plot_format", "png"))))
  if (identical(requested, "both")) requested <- c("png", "pdf")
  if (!length(requested) || any(!requested %in% c("png", "pdf"))) {
    stop("enrichment.plot_format must be 'png', 'pdf', or 'both'")
  }
  unique(requested)
}

save_enrichment_plot <- function(stem, plot, width, height, formats, dpi = 300) {
  for (format in formats) {
    ggplot2::ggsave(paste0(stem, ".", format), plot, width = width, height = height,
                    units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  }
  invisible(NULL)
}

plot_enrichment_summary <- function(x, out, config) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(x)) return(invisible(NULL))
  padj_col <- intersect(c("p.adjust", "pvalue"), names(x))[1]; if (is.na(padj_col)) return(invisible(NULL))
  label_width <- as.integer(cfg_get(config, "enrichment.plot_label_width", 55))
  top_n <- as.integer(cfg_get(config, "enrichment.plot_top_terms", 10))
  ora_fdr <- as.numeric(cfg_get(config, "enrichment.plot_ora_fdr_threshold", 0.05))
  gsea_fdr <- as.numeric(cfg_get(config, "enrichment.plot_gsea_fdr_threshold", 0.25))
  max_gene_overlap <- as.numeric(cfg_get(config, "enrichment.plot_max_gene_overlap", 0.6))
  terms_per_page <- as.integer(cfg_get(config, "enrichment.plot_terms_per_page", 20))
  plot_formats <- enrichment_plot_formats(config)
  plot_dpi <- as.integer(cfg_get(config, "enrichment.plot_dpi", 300))
  if (is.na(top_n) || top_n < 1L) top_n <- 10L
  if (is.na(label_width) || label_width < 15L) label_width <- 55L
  if (is.na(max_gene_overlap) || max_gene_overlap < 0 || max_gene_overlap > 1) max_gene_overlap <- 0.6
  if (is.na(terms_per_page) || terms_per_page < 4L) terms_per_page <- 20L
  if (is.na(plot_dpi) || plot_dpi < 72L) plot_dpi <- 300L

  old_plots <- list.files(out, pattern = "^(enrichment_|gsea_).*[.](pdf|png)$", full.names = TRUE)
  if (length(old_plots)) unlink(old_plots)
  d <- select_enrichment_plot_terms(x, top_n, ora_fdr, gsea_fdr, max_gene_overlap)
  if (!nrow(d)) return(invisible(NULL))
  d$label <- clean_enrichment_label(d$Description, .Machine$integer.max)
  context <- enrichment_plot_context(x, config)
  d$plot_direction_label <- ifelse(d$plot_direction == "Up", context$positive, context$negative)
  d$plot_direction_label <- factor(d$plot_direction_label, levels = c(context$negative, context$positive))
  evidence_levels <- c("FDR <= 0.05", paste0("FDR 0.05-", format(gsea_fdr)))
  d$evidence_class <- ifelse(d$plot_fdr <= 0.05, evidence_levels[[1]], evidence_levels[[2]])
  d$evidence_class <- factor(d$evidence_class, levels = evidence_levels)
  d <- map_enrichment_plot_values(d)
  write_tsv(d, file.path(out, "enrichment_plot_terms_summary.tsv"))
  d$label <- clean_enrichment_label(d$Description, label_width)
  title_suffix <- if (nzchar(context$population)) paste0(" | ", context$population) else ""
  comparison_prefix <- if (nzchar(context$comparison)) paste0(context$comparison, "; ") else ""
  fdr_fill_scale <- ggplot2::scale_fill_gradient(low = "#17365D", high = "#D7E6F2", name = "FDR")
  direction_scale <- ggplot2::scale_color_manual(values = setNames(c("#B35806", "#2166AC"), c(context$negative, context$positive)), name = NULL)
  evidence_scale <- ggplot2::scale_shape_manual(values = stats::setNames(c(16, 1), evidence_levels), name = NULL, drop = FALSE)
  common_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = .35),
                   axis.line.x = ggplot2::element_line(color = "grey35", linewidth = .35),
                   axis.text.x = ggplot2::element_text(color = "#303030", size = 8.5),
                   axis.text.y = ggplot2::element_text(color = "#303030", size = 8),
                   plot.title = ggplot2::element_text(face = "bold", color = "#202020", size = 12),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   strip.text = ggplot2::element_text(face = "bold", color = "#303030", size = 9),
                   legend.text = ggplot2::element_text(size = 8.5), legend.title = ggplot2::element_text(size = 9),
                   legend.position = "top")

  overview_n <- min(3L, top_n)
  overview <- select_enrichment_plot_terms(x, overview_n, ora_fdr, gsea_fdr, max_gene_overlap)
  overview$label <- clean_enrichment_label(overview$Description, min(label_width, 34L))
  overview$plot_direction_label <- ifelse(overview$plot_direction == "Up", context$positive, context$negative)
  overview$plot_direction_label <- factor(overview$plot_direction_label, levels = c(context$negative, context$positive))
  overview$evidence_class <- ifelse(overview$plot_fdr <= 0.05, evidence_levels[[1]], evidence_levels[[2]])
  overview$evidence_class <- factor(overview$evidence_class, levels = evidence_levels)
  overview <- map_enrichment_plot_values(overview)
  gsea_overview <- overview[overview$method == "GSEA" & is.finite(overview$NES), , drop = FALSE]
  if (nrow(gsea_overview)) {
    gsea_overview$fdr_hjust <- ifelse(gsea_overview$NES >= 0, -0.15, 1.15)
    p <- ggplot2::ggplot(gsea_overview, ggplot2::aes(y = stats::reorder(label, NES))) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = stats::reorder(label, NES), color = plot_direction_label), linewidth = .65, alpha = .7) +
      ggplot2::geom_point(ggplot2::aes(x = NES, size = setSize, color = plot_direction_label, shape = evidence_class), stroke = 1.1) +
      ggplot2::geom_text(ggplot2::aes(x = NES, label = plot_fdr_label, hjust = fdr_hjust), size = 2.35, color = "#303030") +
      ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = .45) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.22, .22))) +
      ggplot2::facet_wrap(~database, scales = "free_y", ncol = 3) + direction_scale + evidence_scale + common_theme +
      ggplot2::labs(title = paste0("GSEA overview", title_suffix),
                    subtitle = paste0(comparison_prefix, "top ", overview_n, " non-redundant terms per NES direction; q labels show FDR; displayed FDR <= ", format(gsea_fdr)),
                    x = "Normalized enrichment score (NES)", y = NULL, size = "Gene-set size", shape = NULL) +
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), size = ggplot2::guide_legend(order = 2),
                      shape = ggplot2::guide_legend(order = 3))
    save_enrichment_plot(file.path(out, "enrichment_dotplot_overview"), p, 13.5, 9, plot_formats, plot_dpi)
  }

  ora <- d[d$method == "ORA", , drop = FALSE]
  ora_overview <- overview[overview$method == "ORA", , drop = FALSE]
  if (nrow(ora_overview)) {
    ora_overview$label_panel <- paste(ora_overview$label, ora_overview$database, sep = "|||DB|||")
    p <- ggplot2::ggplot(ora_overview, ggplot2::aes(y = stats::reorder(label_panel, plot_x_value))) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = plot_x_value, yend = stats::reorder(label_panel, plot_x_value), color = plot_direction_label), linewidth = .7, alpha = .75) +
      ggplot2::geom_point(ggplot2::aes(x = plot_x_value, size = plot_size_value, color = plot_direction_label, fill = plot_fdr), shape = 21, stroke = .9) +
      ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = .45) +
      ggplot2::facet_wrap(~database, scales = "free_y", ncol = 3) +
      ggplot2::scale_y_discrete(labels = function(value) sub("[|][|][|]DB[|][|][|].*$", "", value)) +
      direction_scale + fdr_fill_scale + common_theme +
      ggplot2::labs(title = paste0("ORA overview", title_suffix),
                    subtitle = paste0(comparison_prefix, "top ", overview_n, " non-redundant terms per direction; displayed FDR <= ", format(ora_fdr)),
                    x = if ("RichFactor" %in% names(ora)) "Directional rich factor" else "Directional -log10 FDR", y = NULL, size = "Gene count") +
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), fill = ggplot2::guide_colorbar(order = 2),
                      size = ggplot2::guide_legend(order = 3))
    save_enrichment_plot(file.path(out, "enrichment_ora_overview"), p, 13.5, 7.8, plot_formats, plot_dpi)
  }
  if (nrow(ora)) for (idx in split(seq_len(nrow(ora)), ora$database)) {
    full <- ora[idx, , drop = FALSE]
    full <- full[order(factor(full$plot_direction, levels = c("Up", "Down")), full$plot_fdr), , drop = FALSE]
    database <- as.character(full$database[[1]])
    pages <- split(seq_len(nrow(full)), ceiling(seq_len(nrow(full)) / terms_per_page))
    for (page in seq_along(pages)) {
      z <- full[pages[[page]], , drop = FALSE]
      z$plot_direction_label <- factor(ifelse(z$plot_direction == "Up", context$positive, context$negative), levels = c(context$negative, context$positive))
      q <- ggplot2::ggplot(z, ggplot2::aes(y = stats::reorder(label, plot_x_value))) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = plot_x_value, yend = stats::reorder(label, plot_x_value), color = plot_direction_label), linewidth = .75, alpha = .75) +
        ggplot2::geom_point(ggplot2::aes(x = plot_x_value, size = plot_size_value, color = plot_direction_label, fill = plot_fdr), shape = 21, stroke = .9) +
        ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = .45) +
        direction_scale + fdr_fill_scale + common_theme +
        ggplot2::labs(title = paste0(database, " ORA", title_suffix),
                      subtitle = paste0(comparison_prefix, "top ", top_n, " non-redundant terms per direction; FDR <= ", format(ora_fdr), "; page ", page, "/", length(pages)),
                      x = if ("RichFactor" %in% names(z)) "Directional rich factor" else "Directional -log10 FDR", y = NULL, size = "Gene count") +
        ggplot2::guides(color = ggplot2::guide_legend(order = 1), fill = ggplot2::guide_colorbar(order = 2),
                        size = ggplot2::guide_legend(order = 3))
      suffix <- if (length(pages) > 1L) paste0("_page", page) else ""
      save_enrichment_plot(file.path(out, paste0("enrichment_dotplot_", tolower(safe_name(database)), "_ora", suffix)), q,
                           9, enrichment_detail_height(max(table(z$plot_direction))), plot_formats, plot_dpi)
    }
  }

  gsea <- d[d$method == "GSEA" & is.finite(d$NES), , drop = FALSE]
  if (nrow(gsea)) for (idx in split(seq_len(nrow(gsea)), gsea$database)) {
    full <- gsea[idx, , drop = FALSE]
    full <- full[order(factor(full$plot_direction, levels = c("Up", "Down")), full$plot_fdr), , drop = FALSE]
    database <- as.character(full$database[[1]])
    pages <- split(seq_len(nrow(full)), ceiling(seq_len(nrow(full)) / terms_per_page))
    for (page in seq_along(pages)) {
      z <- full[pages[[page]], , drop = FALSE]
      size_breaks <- unique(as.numeric(stats::quantile(z$setSize, c(0, 0.5, 1), na.rm = TRUE)))
      z$fdr_hjust <- ifelse(z$NES >= 0, -0.15, 1.15)
      q <- ggplot2::ggplot(z, ggplot2::aes(y = stats::reorder(label, NES))) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = stats::reorder(label, NES), color = plot_direction_label), linewidth = .75, alpha = .7) +
        ggplot2::geom_point(ggplot2::aes(x = NES, size = setSize, color = plot_direction_label, shape = evidence_class), stroke = 1.1) +
        ggplot2::geom_text(ggplot2::aes(x = NES, label = plot_fdr_label, hjust = fdr_hjust), size = 2.7, color = "#303030") +
        ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = .45) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.22, .22))) +
        ggplot2::scale_size_continuous(breaks = size_breaks, range = c(2, 8)) +
        direction_scale + evidence_scale + common_theme +
        ggplot2::labs(title = paste0(database, " GSEA", title_suffix),
                      subtitle = paste0(comparison_prefix, "top ", top_n, " non-redundant terms per NES direction; q labels show FDR; displayed FDR <= ", format(gsea_fdr), "; page ", page, "/", length(pages)),
                      x = "Normalized enrichment score (NES)", y = NULL, size = "Gene-set size", shape = NULL) +
        ggplot2::guides(color = ggplot2::guide_legend(order = 1), size = ggplot2::guide_legend(order = 2),
                        shape = ggplot2::guide_legend(order = 3))
      suffix <- if (length(pages) > 1L) paste0("_page", page) else ""
      save_enrichment_plot(file.path(out, paste0("gsea_nes_", tolower(safe_name(database)), suffix)), q,
                           10.5, enrichment_detail_height(nrow(z)), plot_formats, plot_dpi)
    }
  }
  invisible(NULL)
}

classify_failure <- function(x) {
  if (grepl("Insufficient independent samples", x)) "skipped_low_replicates" else if (grepl("rank deficient|multiple conditions|not constant", x)) "invalid_design" else if (grepl("required|unavailable", x, ignore.case = TRUE)) "missing_dependency" else "failed"
}
