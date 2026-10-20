#!/usr/bin/env Rscript
# Export DIANA-microT-CDS predictions to data/cache/<GENE>_<species>_diana.csv|rds
#
# Usage:
#   Rscript scripts/export_score_diana.R GENE [species] [flags]
#
# Flags:
#   --threshold=0.7
#   --prefer-api=true|false       (default: true)
#   --mirna-resource=mirbase|mirgenedb

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_score_diana.R GENE [species] [flags]\n",
    "Flags: --threshold=0.7 --prefer-api=true|false --mirna-resource=mirbase\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

threshold <- flag_numeric(flags, "threshold", THRESHOLDS$diana_microt)
prefer_api <- flag_bool(flags, "prefer-api", TRUE)
mirna_resource <- flag_character(flags, "mirna-resource", DIANA$default_mirna_resource)

out <- export_prediction_source(
  "diana",
  parsed$gene,
  parsed$species,
  SOURCE_REGISTRY,
  threshold = threshold,
  prefer_api = prefer_api,
  mirna_resource = mirna_resource
)

message("Wrote ", out$paths$csv, " (", nrow(out$predictions), " rows)")
message("Wrote ", out$paths$rds)
