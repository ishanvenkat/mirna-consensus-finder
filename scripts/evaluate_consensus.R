#!/usr/bin/env Rscript
# Rank-based recall@k and PR curves: tiered consensus vs single tools vs mirDIP.
#
# Usage:
#   Rscript scripts/evaluate_consensus.R GENE [species] [flags]
#
# Prerequisites (cached only by default):
#   Rscript scripts/export_validation.R GENE [species]
#   Rscript scripts/export_score_integrated.R GENE [species] --use-cache=true --method=tiered
#   Per-source prediction caches from export_score_*.R for each tool in --sources
#
# Flags:
#   --support-level=strong|any           (default: strong)
#   --consensus-subset=full|multi_source|both (default: full)
#   --k=10,25,50,100                     (comma-separated)
#   --max-pr-depth=200
#   --sources=targetscan,diana,mirdb,mirdip
#   --plot=true|false
#   --plot-path=docs/presentation/assets/pten-pr-curve.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || args[[1L]] %in% c("-h", "--help")) {
  cat(
    "Usage: Rscript scripts/evaluate_consensus.R GENE [species] [flags]\n",
    "\n",
    "  Writes data/cache/<GENE>_<species>_evaluation_recall_at_k.csv\n",
    "  Writes data/cache/<GENE>_<species>_evaluation_pr_curve.csv\n",
    "  Writes data/cache/<GENE>_<species>_evaluation_summary.rds\n",
    "\n",
    "Flags:\n",
    "  --support-level=strong|any\n",
    "  --consensus-subset=full|multi_source|both\n",
    "  --k=10,25,50,100\n",
    "  --max-pr-depth=200\n",
    "  --sources=targetscan,diana,mirdb,mirdip\n",
    "  --plot=true|false\n",
    "  --plot-path=<path>\n",
  sep = ""
  )
  quit(status = if (length(args) == 0L) 1L else 0L)
}

source("scripts/lib/export_cli.R")
setup_export_project()
parsed <- positional_gene_species(args)
flags <- parsed$flags

gene <- validate_gene_symbol(parsed$gene)
species <- parsed$species

support_level <- flag_character(flags, "support-level", "strong")
if (!support_level %in% c("strong", "any")) {
  stop("--support-level must be 'strong' or 'any'", call. = FALSE)
}

consensus_subset <- flag_character(flags, "consensus-subset", "full")
if (!consensus_subset %in% c("full", "multi_source", "both")) {
  stop(
    "--consensus-subset must be full, multi_source, or both",
    call. = FALSE
  )
}

k_raw <- flag_character(flags, "k", "10,25,50,100")
k_values <- as.integer(strsplit(k_raw, ",", fixed = TRUE)[[1L]])
k_values <- k_values[!is.na(k_values)]
if (length(k_values) == 0L) {
  stop("Invalid --k values", call. = FALSE)
}

max_pr_depth <- flag_numeric(flags, "max-pr-depth", 200L)
if (is.null(max_pr_depth) || max_pr_depth < 1L) {
  stop("--max-pr-depth must be a positive integer", call. = FALSE)
}
max_pr_depth <- as.integer(max_pr_depth)

sources_raw <- flag_character(flags, "sources", "targetscan,diana,mirdb,mirdip")
prediction_sources <- strsplit(sources_raw, ",", fixed = TRUE)[[1L]]
prediction_sources <- trimws(prediction_sources)
prediction_sources <- prediction_sources[nzchar(prediction_sources)]

do_plot <- flag_bool(flags, "plot", FALSE)
default_plot <- file.path(
  "docs",
  "presentation",
  "assets",
  paste0(tolower(gene), "-pr-curve.png")
)
plot_path <- flag_character(flags, "plot-path", default_plot)

source("R/evaluate_consensus.R")

result <- evaluate_consensus_rankings(
  gene_symbol = gene,
  species = species,
  support_level = support_level,
  consensus_subset = consensus_subset,
  k_values = k_values,
  max_pr_depth = max_pr_depth,
  prediction_sources = prediction_sources
)

paths <- write_evaluation_cache(result, gene, species)
message("Wrote ", paths$recall_csv)
message("Wrote ", paths$pr_curve_csv)
message("Wrote ", paths$summary_rds)

message(
  "Ground-truth positives (",
  result$metadata$support_level,
  "): ",
  result$metadata$n_positives
)
for (mid in names(result$rank_lengths)) {
  message("  ", mid, ": ranked length ", result$rank_lengths[[mid]])
}

if (isTRUE(do_plot)) {
  plot_file <- plot_evaluation_pr_curves(result$pr_curve, plot_path)
  message("Wrote PR plot: ", plot_file)
}