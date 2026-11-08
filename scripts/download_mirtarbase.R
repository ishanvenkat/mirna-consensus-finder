#!/usr/bin/env Rscript
# Download miRTarBase 2025 bulk CSV files (awi.cuhk.edu.cn release 11.0).
#
# Usage:
#   Rscript scripts/download_mirtarbase.R [--force] [--full]
#
# Default: reporter-assay subset (miRTarBase_SE_R.csv, ~2 MB).
# --full: additionally download full MTI (miRTarBase_MTI.csv, ~1.7 GB).

args <- commandArgs(trailingOnly = TRUE)
force_download <- "--force" %in% args
download_full <- "--full" %in% args

if (file.exists("R/config.R")) {
  project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)
} else if (file.exists("../R/config.R")) {
  project_root <- normalizePath("..", winslash = "/", mustWork = FALSE)
} else {
  stop("Run from project root: Rscript scripts/download_mirtarbase.R", call. = FALSE)
}
setwd(project_root)

source("R/config.R")
ensure_data_dirs()

download_bulk <- function(file_key, force = FALSE) {
  url <- mirtarbase_bulk_url(file_key)
  dest <- mirtarbase_bulk_path(file_key)
  if (file.exists(dest) && !force && file.info(dest)$size > 0L) {
    message("  skip (exists): ", basename(dest))
    return(invisible(dest))
  }

  message("  downloading: ", url)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, 7200))

  parent <- dirname(dest)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }

  status <- tryCatch({
    if (requireNamespace("httr2", quietly = TRUE)) {
      resp <- httr2::request(url) |>
        httr2::req_user_agent("mirna-consensus-finder/1.0 (Rscript)") |>
        httr2::req_perform()
      httr2::resp_check_status(resp)
      writeBin(httr2::resp_body_raw(resp), dest)
      0L
    } else {
      old <- options(
        download.file.method = "libcurl",
        HTTPUserAgent = "mirna-consensus-finder/1.0 (Rscript)"
      )
      on.exit(options(old), add = TRUE)
      download.file(url, dest, mode = "wb", quiet = TRUE)
    }
  }, error = function(e) e)
  if (inherits(status, "error") || !identical(status, 0L)) {
    stop("Download failed for ", url, ": ", conditionMessage(status), call. = FALSE)
  }
  if (!file.exists(dest) || file.info(dest)$size == 0L) {
    stop("Downloaded file missing or empty: ", dest, call. = FALSE)
  }
  message("  wrote: ", dest, " (", file.info(dest)$size, " bytes)")
  invisible(dest)
}

message("miRTarBase ", MIRTARBASE$release_label, " (DB files ", MIRTARBASE$db_version, ")")
download_bulk("se_reporter", force = force_download)
if (isTRUE(download_full)) {
  download_bulk("mti", force = force_download)
} else {
  message(
    "  tip: use --full to download complete MTI (large); ",
    "reporter subset is enough for reporter/Western strong-evidence filtering."
  )
}
