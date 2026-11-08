#!/usr/bin/env Rscript
# Join integrated consensus cache with validation overlay.
#
# Usage:
#   Rscript scripts/export_integrated_annotated.R GENE [species] [flags]
#
# Flags:
#   --use-integrated-cache=true   (default: true)
#   --use-validation-cache=true   (default: true; else query_validation live)
#   --append-validated=true|false (default: false)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_integrated_annotated.R GENE [species] [flags]\n",
    "\n",
    "  Writes data/cache/<GENE>_<species>_integrated_annotated.csv|.rds\n",
    "\n",
    "Requires integrated cache from export_score_integrated.R unless disabled.\n",
    "Requires validation cache from export_validation.R when --use-validation-cache=true.\n",
    "\n",
    "Flags:\n",
    "  --use-integrated-cache=true|false\n",
    "  --use-validation-cache=true|false\n",
    "  --append-validated=true|false\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

use_integrated_cache <- flag_bool(flags, "use-integrated-cache", TRUE)
use_validation_cache <- flag_bool(flags, "use-validation-cache", TRUE)
append_validated <- flag_bool(flags, "append-validated", FALSE)

gene <- validate_gene_symbol(parsed$gene)
species <- parsed$species

source("R/config.R")
source("R/module_validation.R")
source("R/cache_predictions.R")

ensure_data_dirs()

if (isTRUE(use_integrated_cache)) {
  if (!integrated_cache_exists(gene, species)) {
    stop(
      "Integrated cache missing. Run: Rscript scripts/export_score_integrated.R ",
      gene, " ", species,
      call. = FALSE
    )
  }
  integrated <- read_integrated_cache(gene, species)
} else {
  stop("--use-integrated-cache=false is not supported yet", call. = FALSE)
}

if (isTRUE(use_validation_cache)) {
  if (!validation_cache_exists(gene, species)) {
    stop(
      "Validation cache missing. Run: Rscript scripts/export_validation.R ",
      gene, " ", species,
      call. = FALSE
    )
  }
  validation <- read_validation_cache(gene, species)
} else {
  validation <- query_validation(gene, species = species)
}

annotated <- annotate_integrated_with_validation(
  integrated,
  validation,
  append_validated_only = append_validated
)

write_integrated_annotated_cache(annotated, gene, species)

supported <- sum(annotated$experimental_support, na.rm = TRUE)
message(
  "Annotated rows: ", nrow(annotated),
  " (experimental_support=TRUE on ", supported, " rows)"
)
