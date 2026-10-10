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
