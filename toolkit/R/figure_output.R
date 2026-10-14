# One scientific figure format per run. Internal raster layers are not exports.
figure_format <- function(config) {
  get <- function(path) {
    x <- config
    for (key in strsplit(path, "\\.")[[1]]) {
      if (!is.list(x) || is.null(x[[key]])) return(NULL)
      x <- x[[key]]
    }
    x
  }
  fields <- c("output.figure_format", "plots.figure_format", "enrichment.plot_format", "visualization.figure_format")
  choices <- character()
  for (field in fields) {
    x <- get(field)
    if (is.null(x)) next
    if (!is.character(x) || length(x) != 1L || !tolower(x) %in% c("pdf", "png")) stop(field, " must select exactly one format: pdf or png")
    choices <- c(choices, tolower(x))
  }
  preview <- get("output.preview_png")
  if (!is.null(preview)) {
    if (!is.logical(preview) || length(preview) != 1L || is.na(preview)) stop("output.preview_png must be true or false")
    if (preview) choices <- c(choices, "png")
  }
  if (length(unique(choices)) > 1L) stop("Conflicting figure formats; select one format for the entire run")
  if (length(choices)) choices[[1]] else "pdf"
}

figure_device <- function(stem, width, height, config) {
  format <- figure_format(config)
  if (format == "pdf") grDevices::pdf(paste0(stem, ".pdf"), width=width, height=height, onefile=TRUE)
  else grDevices::png(paste0(stem, "_page%03d.png"), width=width, height=height, units="in", res=300)
}

# ggplot/Seurat point sizes are linear; matplotlib sizes are areas.
umap_point_size <- function(n_cells, override=NULL) {
  if(!is.null(override)) {
    if(!is.numeric(override) || length(override)!=1L || !is.finite(override) || override<=0) stop("UMAP point size must be a positive finite number")
    return(override)
  }
  if(length(n_cells)!=1L || !is.finite(n_cells) || n_cells<1) stop("UMAP requires at least one cell")
  min(1,max(.08,6/sqrt(n_cells)))
}

umap_palette <- function(labels, configured=NULL) {
  labels <- sort(unique(as.character(labels[!is.na(labels)])),method="radix")
  configured <- unlist(configured)
  key <- function(x) apply(grDevices::col2rgb(x,alpha=TRUE),2,paste,collapse=",")
  if(length(configured)) {
    if(is.null(names(configured)) || any(!nzchar(names(configured))) || anyDuplicated(names(configured))) stop("Categorical colour maps need unique label names")
    if(anyDuplicated(key(configured))) stop("Duplicate categorical colours, including equivalent colour aliases")
  }
  if(length(labels)>20L) warning("More than 20 categories: unique colours do not guarantee perceptual distinguishability; retain labels and inspect the plot")
  colors <- setNames(rep(NA_character_,length(labels)),labels)
  common <- intersect(labels,names(configured));colors[common] <- configured[common]
  used <- if(length(common)) unname(key(colors[common])) else character()
  for(label in setdiff(labels,common)) {
    hash <- 0
    for(code in utf8ToInt(enc2utf8(label))) hash <- (hash*31+code)%%104729
    attempt <- 0L
    repeat {
      color <- grDevices::hcl(h=(hash*137.508+attempt*37.1)%%360,c=65,l=58+(attempt%%3)*5,fixup=TRUE)
      signature <- unname(key(color))
      if(!signature %in% used) break
      attempt <- attempt+1L
      if(attempt>10000L) stop("Unable to construct unique categorical colours")
    }
    colors[label] <- color;used <- c(used,signature)
  }
  colors
}
