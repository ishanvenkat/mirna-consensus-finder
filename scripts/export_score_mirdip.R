#!/usr/bin/env Rscript
# Export mirDIP predictions to data/cache/<GENE>_<species>_mirdip.csv|rds
#
# Usage:
#   Rscript scripts/export_score_mirdip.R GENE [species] [flags]
#
# Flags:
#   --score-class=Very High|High|Medium|Low
#   --min-score=<numeric>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_score_mirdip.R GENE [species] [flags]\n",
    "Flags: --score-class=Low --min-score=N\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags
score_class <- flag_character(flags, "score-class", MIRDIP$default_score_class)
min_score <- flag_numeric(flags, "min-score", THRESHOLDS$mirdip_min_score)

out <- export_prediction_source(
  "mirdip",
  parsed$gene,
  parsed$species,
  SOURCE_REGISTRY,
  score_class = score_class,
  min_score = min_score
)

message("Wrote ", out$paths$csv, " (", nrow(out$predictions), " rows)")
message("Wrote ", out$paths$rds)
