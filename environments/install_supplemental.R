args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
manifest <- file.path(dirname(normalizePath(script)), "supplemental-r.json")
install <- "--install" %in% args
if (!requireNamespace("jsonlite", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
  stop("Deploy the locked Pixi R environment before supplemental packages")
}
specs <- jsonlite::fromJSON(manifest, simplifyVector = FALSE)$packages
option <- function(name, default = NULL) {
  values <- args[startsWith(args, paste0("--", name, "="))]
  if (length(values)) sub(paste0("^--", name, "="), "", values[[1]]) else default
}
library <- option("library", .Library)
cache <- option("cache")
if (install) dir.create(library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library, .libPaths()))
available <- function(spec) {
  description <- file.path(library, spec$package, "DESCRIPTION")
  if (!file.exists(description)) return(FALSE)
  d <- read.dcf(description)
  if (!identical(unname(d[1, "Version"]), spec$version)) return(FALSE)
  if (!is.null(spec$commit)) {
    stamp <- file.path(library, spec$package, ".scrna-source.json")
    revision <- if ("RemoteSha" %in% colnames(d)) d[1, "RemoteSha"] else if (file.exists(stamp)) jsonlite::fromJSON(stamp)$commit else ""
    if (!identical(unname(revision), spec$commit)) return(FALSE)
  }
  TRUE
}
failed <- character()
for (spec in specs) {
  if (!available(spec) && install) {
    archive <- tempfile(fileext = ".tar.gz")
    downloaded <- FALSE
    cached <- if (!is.null(cache)) file.path(cache, spec$archive) else ""
    if (nzchar(cached) && file.exists(cached)) {
      if (!identical(digest::digest(cached, algo = spec$checksum_algorithm, file = TRUE), spec$checksum)) stop("Cached source checksum mismatch: ", spec$package)
      file.copy(cached, archive)
      downloaded <- TRUE
    }
    for (url in spec$urls) {
      if (downloaded) break
      status <- tryCatch(download.file(url, archive, mode = "wb", quiet = TRUE), error = function(e) 1L)
      if (status == 0L && identical(digest::digest(archive, algo = spec$checksum_algorithm, file = TRUE), spec$checksum)) {
        downloaded <- TRUE
        break
      }
    }
    if (!downloaded) stop("No source archive passed checksum verification: ", spec$package)
    # 仅安装已校验的归档；依赖由 pixi.lock 和本清单按顺序提供。
    withCallingHandlers(install.packages(archive, repos = NULL, type = "source", lib = library, dependencies = FALSE),
      warning = function(w) stop(conditionMessage(w)))
    description <- file.path(library, spec$package, "DESCRIPTION")
    if (!file.exists(description) || !identical(unname(read.dcf(description)[1, "Version"]), spec$version)) stop("Installed package version mismatch: ", spec$package)
    jsonlite::write_json(spec, file.path(library, spec$package, ".scrna-source.json"), auto_unbox = TRUE, pretty = TRUE)
    unlink(archive)
  }
  ok <- available(spec)
  cat(spec$package, spec$version, if (ok) "verified" else "missing_or_different", "\n")
  if (!ok) failed <- c(failed, spec$package)
}
if (length(failed)) stop("Supplemental packages need explicit deployment: ", paste(failed, collapse = ", "))
