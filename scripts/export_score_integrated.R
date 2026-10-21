#!/usr/bin/env Rscript
# Run find_mirna_regulators with integration and write cache outputs.
#
# Usage:
#   Rscript scripts/export_score_integrated.R GENE [species] [flags]
#
# Flags:
#   --method=tiered|rank_fusion|percentile   (default: tiered)
#   --use-cache=true|false                   (default: false)
#   --sources=targetscan,diana,mirdb,mirdip  (comma-separated)
#   --write-combined=true|false              (default: true)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/export_score_integrated.R GENE [species] [flags]\n",
    "\n",
    "  Writes data/cache/<GENE>_<species>_integrated.csv|.rds\n",
    "  Writes data/cache/<GENE>_<species>_integrated_consensus.csv|.rds (n_sources >= 2)\n",
    "  Writes data/cache/<GENE>_<species>_integrated_single_source.csv|.rds (n_sources == 1)\n",
    "  Optionally data/cache/<GENE>_<species>_predictions_combined.csv\n",
    "\n",
    "Flags:\n",
    "  --method=tiered|rank_fusion|percentile\n",
    "  --use-cache=true|false\n",
    "  --sources=targetscan,diana,mirdb,mirdip\n",
    "  --write-combined=true|false\n",
    sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

method <- flag_character(flags, "method", "tiered")
use_cache <- flag_bool(flags, "use-cache", FALSE)
write_combined <- flag_bool(flags, "write-combined", TRUE)
sources_raw <- flag_character(flags, "sources", NULL)
sources <- if (is.null(sources_raw)) {
  c("targetscan", "diana", "mirdb", "mirdip")
} else {
  strsplit(sources_raw, ",", fixed = TRUE)[[1L]]
}
sources <- trimws(sources)
sources <- sources[nzchar(sources)]

gene <- validate_gene_symbol(parsed$gene)
species <- parsed$species

out <- find_mirna_regulators(
  gene,
  species = species,
  sources = sources,
  integrate = TRUE,
  integrate_method = method,
  use_cached_predictions = use_cache,
  include_validation = FALSE
)

if (is.null(out$integrated) || nrow(out$integrated) == 0L) {
  stop("Integration returned no results", call. = FALSE)
}

paths <- write_integrated_cache(out$integrated, gene, species)
message("Wrote ", paths$csv, " (", nrow(out$integrated), " rows)")
message("Wrote ", paths$rds)
n_multi <- sum(out$integrated$n_sources >= 2L)
n_single <- sum(out$integrated$n_sources == 1L)
message(
  "Wrote ", paths$consensus_csv, " (", n_multi, " rows, n_sources >= 2)"
)
message("Wrote ", paths$consensus_rds)
message(
  "Wrote ", paths$single_source_csv, " (", n_single, " rows, n_sources == 1)"
)
message("Wrote ", paths$single_source_rds)

if (isTRUE(write_combined)) {
  combined_path <- write_combined_predictions_cache(out$predictions, gene, species)
  message("Wrote ", combined_path, " (", nrow(out$predictions), " rows)")
}
