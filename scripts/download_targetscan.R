#!/usr/bin/env Rscript
# Download TargetScan 8.0 bulk prediction files.
# See: https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi

args <- commandArgs(trailingOnly = TRUE)
force_download <- "--force" %in% args

if (file.exists("R/config.R")) {
  project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)
} else if (file.exists("../R/config.R")) {
  project_root <- normalizePath("..", winslash = "/", mustWork = FALSE)
} else {
  stop("Run from project root: Rscript scripts/download_targetscan.R", call. = FALSE)
}
setwd(project_root)

source("R/config.R")
ensure_data_dirs()

download_zip <- function(url, zip_path, extract_path, force = FALSE) {
  if (file.exists(extract_path) && !force && file.info(extract_path)$size > 0L) {
    message("  skip (exists): ", basename(extract_path))
    return(invisible(extract_path))
  }

  message("  downloading: ", basename(url))
  tmp_zip <- tempfile(fileext = ".zip")
  on.exit(unlink(tmp_zip), add = TRUE)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, 3600))

  status <- tryCatch(
    download.file(url, tmp_zip, mode = "wb", quiet = TRUE),
    error = function(e) e
  )
  if (inherits(status, "error") || !identical(status, 0L)) {
    stop("Download failed for ", url, call. = FALSE)
  }

  utils::unzip(tmp_zip, exdir = dirname(extract_path))
  if (!file.exists(extract_path)) {
    stop("Expected file not found after unzip: ", extract_path, call. = FALSE)
  }
  invisible(extract_path)
}

download_species <- function(db_key, force = FALSE) {
  base_url <- TARGETSCAN$base_url[[db_key]]
  if (is.null(base_url)) {
    stop("Unknown TargetScan database key: ", db_key, call. = FALSE)
  }

  dest_dir <- file.path(PATHS$targetscan, db_key)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }

  message("\nTargetScan ", db_key, " (", base_url, ")")

  pairs <- list(
    list(zip = TARGETSCAN$zip_files$summary_all, tsv = TARGETSCAN$files$summary_all),
    list(zip = TARGETSCAN$zip_files$summary_default, tsv = TARGETSCAN$files$summary_default),
    list(zip = TARGETSCAN$zip_files$gene_info, tsv = TARGETSCAN$files$gene_info),
    list(zip = TARGETSCAN$zip_files$family_info, tsv = TARGETSCAN$files$family_info)
  )

  for (pair in pairs) {
    url <- paste0(base_url, "/", pair$zip)
    extract_path <- file.path(dest_dir, pair$tsv)
    download_zip(url, file.path(dest_dir, pair$zip), extract_path, force = force)
  }

  summary_path <- file.path(dest_dir, TARGETSCAN$files$summary_all)
  row_count <- if (file.exists(summary_path)) {
    length(readLines(summary_path, warn = FALSE)) - 1L
  } else {
    NA_integer_
  }

  list(db = db_key, dest_dir = dest_dir, summary_rows = row_count)
}

message("TargetScan ", TARGETSCAN$version, " bulk download")
if (force_download) {
  message("Force mode: re-downloading all files")
}

results <- list(
  download_species("vert_80", force = force_download),
  download_species("mmu_80", force = force_download)
)

message("\nDone. Files saved under: ", PATHS$targetscan)
for (res in results) {
  message(
    "  ", res$db, ": ", res$dest_dir,
    if (!is.na(res$summary_rows)) paste0(" (~", res$summary_rows, " summary rows)") else ""
  )
}
