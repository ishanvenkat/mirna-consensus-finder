#!/usr/bin/env Rscript
# Download miRDB v6.0 bulk predicted targets (all species in one file).
# See: https://mirdb.org/download.html

args <- commandArgs(trailingOnly = TRUE)
force_download <- "--force" %in% args

if (file.exists("R/config.R")) {
  project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)
} else if (file.exists("../R/config.R")) {
  project_root <- normalizePath("..", winslash = "/", mustWork = FALSE)
} else {
  stop("Run from project root: Rscript scripts/download_mirdb.R", call. = FALSE)
}
setwd(project_root)

source("R/config.R")
ensure_data_dirs()

dest_gz <- file.path(PATHS$mirdb, MIRDB$bulk_file)
url <- MIRDB$download_url

if (file.exists(dest_gz) && !force_download && file.info(dest_gz)$size > 0L) {
  message("skip (exists): ", dest_gz)
} else {
  message("Downloading miRDB ", MIRDB$version, " bulk predictions (~large file)...")
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
  message("Saved: ", dest_gz)
}

con <- gzfile(dest_gz, open = "rt")
on.exit(close(con), add = TRUE)
sample <- readLines(con, n = 3L, warn = FALSE)

message("\nDone. miRDB bulk file: ", dest_gz)
message("Format (no header): miRNA \\t RefSeq_accession \\t MirTarget_score")
message("Query by HUGO symbol uses RefSeq mapping cache: data/cache/refseq_symbol_map.tsv")
if (length(sample) > 0L) {
  message("Sample row: ", sample[[1L]])
}
