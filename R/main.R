# Main entry point for miRNA target prediction.

source("R/config.R")
source("R/utils_gene_id.R")
source("R/module_targetscan.R")
source("R/module_diana.R")
source("R/module_mirdb.R")
source("R/module_mirdip.R")
source("R/module_validation.R")
source("R/integrate_scores.R")
source("R/cache_predictions.R")

ensure_data_dirs()

SOURCE_REGISTRY <- list(
  targetscan = query_targets_targetscan,
  diana      = query_targets_diana,
  mirdb      = query_targets_mirdb,
  mirdip     = query_targets_mirdip
)

find_mirna_regulators <- function(
  gene_symbol,
  species = "hsa",
  sources = c("targetscan", "diana", "mirdb", "mirdip"),
  include_validation = TRUE,
  integrate = FALSE,
  integrate_method = c("rank_fusion", "percentile", "tiered"),
  use_cached_predictions = FALSE,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  integrate_method <- match.arg(integrate_method)

  unknown <- setdiff(sources, names(SOURCE_REGISTRY))
  if (length(unknown) > 0L) {
    stop("Unknown sources: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  predictions <- load_predictions(
    gene_symbol = gene_symbol,
    species = species,
    sources = sources,
    use_cached_predictions = use_cached_predictions,
    source_registry = SOURCE_REGISTRY,
    ...
  )

  integrated <- NULL
  if (isTRUE(integrate)) {
    integrated <- integrate_scores(predictions, method = integrate_method)
    fusion_sources <- if (integrate_method == "tiered") {
      INDEPENDENT_SOURCES
    } else {
      unique(predictions$source)
    }
    fusion_sources <- intersect(fusion_sources, unique(predictions$source))
    integrated <- annotate_integrated(
      integrated,
      predictions,
      fusion_sources = fusion_sources
    )
  }

  validation <- NULL
  if (isTRUE(include_validation)) {
    validation <- query_validation(gene_symbol, species = species, ...)
  }

  list(
    predictions = predictions,
    integrated = integrated,
    validation = validation
  )
}
