#!/usr/bin/env Rscript
# Download DIANA-microT 2023 bulk interaction files (miRBase).
# See: https://dianalab.e-ce.uth.gr/microt_webserver/

args <- commandArgs(trailingOnly = TRUE)
force_download <- "--force" %in% args
species_arg <- args[args != "--force"]

if (file.exists("R/config.R")) {
  project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)
} else if (file.exists("../R/config.R")) {
  project_root <- normalizePath("..", winslash = "/", mustWork = FALSE)
} else {
  stop("Run from project root: Rscript scripts/download_diana_bulk.R", call. = FALSE)
}
setwd(project_root)

source("R/config.R")
ensure_data_dirs()

download_gz <- function(url, dest_gz, force = FALSE) {
  if (file.exists(dest_gz) && !force && file.info(dest_gz)$size > 0L) {
    message("  skip (exists): ", basename(dest_gz))
    return(invisible(dest_gz))
  }

  message("  downloading: ", basename(dest_gz))
  tmp <- tempfile(fileext = ".gz")
  on.exit(unlink(tmp), add = TRUE)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, 7200))

  status <- tryCatch(
    download.file(url, tmp, mode = "wb", quiet = TRUE),
    error = function(e) e
  )
  if (inherits(status, "error") || !identical(status, 0L)) {
    stop("Download failed for ", url, call. = FALSE)
  }

  if (!dir.exists(dirname(dest_gz))) {
    dir.create(dirname(dest_gz), recursive = TRUE, showWarnings = FALSE)
  }
  file.copy(tmp, dest_gz, overwrite = TRUE)
  invisible(dest_gz)
}

download_species <- function(species, force = FALSE) {
  if (!species %in% names(DIANA$bulk_files)) {
    stop("Unknown species code: ", species, call. = FALSE)
  }

  dest_dir <- PATHS$diana
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }

  resource <- DIANA$default_mirna_resource
  filename <- DIANA$bulk_files[[species]][[resource]]
  url <- paste0(DIANA$downloads_base, "/", filename)
  dest_gz <- file.path(dest_dir, filename)

  message("\nDIANA ", species, " (", filename, ")")
  download_gz(url, dest_gz, force = force)

  con <- gzfile(dest_gz, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)

  list(species = species, path = dest_gz, header = header)
}

species_to_fetch <- if (length(species_arg) > 0L) {
  species_arg
} else {
  c("hsa", "mmu")
}

message("DIANA-microT bulk download (miRBase flat files)")
if (force_download) {
  message("Force mode: re-downloading files")
}

results <- lapply(species_to_fetch, download_species, force = force_download)

message("\nDone. Files saved under: ", PATHS$diana)
for (res in results) {
  message("  ", res$species, ": ", res$path)
  message("    header: ", res$header)
}
