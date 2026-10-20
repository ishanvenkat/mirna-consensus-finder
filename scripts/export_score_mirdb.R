#!/usr/bin/env Rscript
# Export miRDB predictions to data/cache/<GENE>_<species>_mirdb.csv|rds
#
# Usage:
#   Rscript scripts/export_score_mirdb.R GENE [species] [flags]
#
# Flags:
#   --min-score=60

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_score_mirdb.R GENE [species] [--min-score=60]\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags
min_score <- flag_numeric(flags, "min-score", THRESHOLDS$mirdb)

out <- export_prediction_source(
  "mirdb",
  parsed$gene,
  parsed$species,
  SOURCE_REGISTRY,
  min_score = min_score
)

message("Wrote ", out$paths$csv, " (", nrow(out$predictions), " rows)")
message("Wrote ", out$paths$rds)
