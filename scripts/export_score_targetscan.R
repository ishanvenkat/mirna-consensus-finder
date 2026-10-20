#!/usr/bin/env Rscript
# Export TargetScan predictions to data/cache/<GENE>_<species>_targetscan.csv|rds
#
# Usage:
#   Rscript scripts/export_score_targetscan.R GENE [species] [flags]
#
# Flags:
#   --predictions=all|default     (default: all)
#   --score-name=cwcs|biochem_occupancy
#   --min-cwcs=<numeric>          CWCS <= value (more negative = stronger)
#   --mirna-prefix=<prefix|species>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_score_targetscan.R GENE [species] [flags]\n",
    "Flags: --predictions=all|default --score-name=cwcs|biochem_occupancy ",
    "--min-cwcs=N --mirna-prefix=hsa|species\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

predictions <- flag_character(flags, "predictions", "all")
score_name <- flag_character(flags, "score-name", "cwcs")
min_cwcs <- flag_numeric(flags, "min-cwcs", NULL)
mirna_prefix <- flag_character(flags, "mirna-prefix", NULL)

out <- export_prediction_source(
  "targetscan",
  parsed$gene,
  parsed$species,
  SOURCE_REGISTRY,
  predictions = predictions,
  score_name = score_name,
  min_cwcs = min_cwcs,
  mirna_prefix = mirna_prefix
)

message("Wrote ", out$paths$csv, " (", nrow(out$predictions), " rows)")
message("Wrote ", out$paths$rds)
