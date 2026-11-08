#!/usr/bin/env Rscript
# Query validation databases and write cache outputs.
#
# Usage:
#   Rscript scripts/export_validation.R GENE [species] [flags]
#
# Flags:
#   --skip-download=true|false   (default: false — run download if bulk missing)
#   --sources=mirtarbase,tarbase (comma-separated)
#   --mirtarbase-full=true       (pass --full to download script)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_validation.R GENE [species] [flags]\n",
    "\n",
    "  Writes data/cache/<GENE>_<species>_validation.csv|.rds\n",
    "  Writes split: ..._validation_mirtarbase_strong.csv, ..._validation_tarbase.csv\n",
    "\n",
    "Flags:\n",
    "  --skip-download=true|false\n",
    "  --sources=mirtarbase,tarbase\n",
    "  --mirtarbase-full=true\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

skip_download <- flag_bool(flags, "skip-download", FALSE)
mirtarbase_full <- flag_bool(flags, "mirtarbase-full", FALSE)
sources_raw <- flag_character(flags, "sources", "mirtarbase,tarbase")
sources <- strsplit(sources_raw, ",", fixed = TRUE)[[1L]]
sources <- trimws(sources)
sources <- sources[nzchar(sources)]

gene <- validate_gene_symbol(parsed$gene)
species <- parsed$species

source("R/config.R")
source("R/module_validation.R")
source("R/cache_predictions.R")

ensure_data_dirs()

needs_bulk <- "mirtarbase" %in% sources
if (needs_bulk && !isTRUE(skip_download)) {
  if (is.null(resolve_mirtarbase_bulk_path(TRUE))) {
    message("miRTarBase bulk missing; running scripts/download_mirtarbase.R")
    dl_args <- character()
    if (isTRUE(mirtarbase_full)) {
      dl_args <- c(dl_args, "--full")
    }
    status <- system2(
      "Rscript",
      c("scripts/download_mirtarbase.R", dl_args),
      stdout = "",
      stderr = ""
    )
    if (!identical(status, 0L)) {
      stop("download_mirtarbase.R failed", call. = FALSE)
    }
  }
}

validation <- query_validation(
  gene,
  species = species,
  sources = sources,
  mirtarbase_use_local = TRUE,
  mirtarbase_strong_methods = "reporter_western"
)

paths <- write_validation_cache(validation, gene, species)
message("Wrote ", paths$csv, " (", nrow(validation), " rows)")
message("Wrote ", paths$rds)

if (nrow(validation) > 0L) {
  mt <- validation[grepl("miRTarBase", validation$validation_sources), , drop = FALSE]
  tb <- validation[grepl("TarBase", validation$validation_sources), , drop = FALSE]
  if (nrow(mt) > 0L) {
    mt_csv <- validation_cache_path(gene, species, "validation_mirtarbase_strong", "csv")
    utils::write.csv(mt, mt_csv, row.names = FALSE)
    message("Wrote ", mt_csv, " (", nrow(mt), " miRTarBase rows)")
  }
  if (nrow(tb) > 0L) {
    tb_csv <- validation_cache_path(gene, species, "validation_tarbase", "csv")
    utils::write.csv(tb, tb_csv, row.names = FALSE)
    message("Wrote ", tb_csv, " (", nrow(tb), " TarBase rows)")
  }
}

prov <- attr(validation, "provenance")
if (!is.null(prov)) {
  message("Provenance:")
  message("  mirtarbase_release: ", prov$mirtarbase_release %||% NA)
  message("  multimir_db_version: ", prov$multimir_db_version %||% NA)
  message("  tarbase_snapshot: ", prov$tarbase_snapshot %||% NA)
}
