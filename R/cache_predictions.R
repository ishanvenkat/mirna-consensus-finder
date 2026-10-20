# Prediction cache I/O — standardized paths under data/cache/.

source("R/config.R")
source("R/utils_gene_id.R")

PREDICTION_CACHE_SOURCES <- c("targetscan", "diana", "mirdb", "mirdip")

PREDICTION_CACHE_REQUIRED_COLS <- c(
  "gene_symbol", "mirna", "score_raw", "score_direction", "source"
)

validate_prediction_source <- function(source) {
  if (!source %in% PREDICTION_CACHE_SOURCES) {
    stop(
      "Unknown prediction source: ", source,
      ". Expected one of: ", paste(PREDICTION_CACHE_SOURCES, collapse = ", "),
      call. = FALSE
    )
  }
  source
}

validate_cache_species <- function(species) {
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, ". Use 'hsa' or 'mmu'.", call. = FALSE)
  }
  species
}

prediction_cache_path <- function(gene_symbol, species, source, ext = "csv") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  source <- validate_prediction_source(source)
  ext <- match.arg(ext, c("csv", "rds"))
  file.path(
    PATHS$cache,
    paste0(gene_symbol, "_", species, "_", source, ".", ext)
  )
}

combined_predictions_cache_path <- function(gene_symbol, species) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  file.path(PATHS$cache, paste0(gene_symbol, "_", species, "_predictions_combined.csv"))
}

integrated_cache_path <- function(gene_symbol, species, ext = "csv") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  ext <- match.arg(ext, c("csv", "rds"))
  file.path(PATHS$cache, paste0(gene_symbol, "_", species, "_integrated.", ext))
}

validate_prediction_frame <- function(df, path_hint = "predictions") {
  if (!is.data.frame(df)) {
    stop(path_hint, " must be a data.frame", call. = FALSE)
  }
  missing <- setdiff(PREDICTION_CACHE_REQUIRED_COLS, names(df))
  if (length(missing) > 0L) {
    stop(path_hint, " missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  df
}

write_prediction_cache <- function(df, gene_symbol, species, source) {
  df <- validate_prediction_frame(df, "prediction cache")
  source <- validate_prediction_source(source)
  if (!all(df$source == source)) {
    stop("prediction cache source column must match: ", source, call. = FALSE)
  }
  ensure_data_dirs()
  csv_path <- prediction_cache_path(gene_symbol, species, source, "csv")
  rds_path <- prediction_cache_path(gene_symbol, species, source, "rds")
  write.csv(df, csv_path, row.names = FALSE)
  saveRDS(df, rds_path)
  invisible(list(csv = csv_path, rds = rds_path))
}

read_prediction_cache <- function(gene_symbol, species, source) {
  rds_path <- prediction_cache_path(gene_symbol, species, source, "rds")
  csv_path <- prediction_cache_path(gene_symbol, species, source, "csv")
  if (file.exists(rds_path)) {
    df <- readRDS(rds_path)
  } else if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    stop("Prediction cache not found: ", csv_path, call. = FALSE)
  }
  validate_prediction_frame(df, csv_path)
}

prediction_cache_exists <- function(gene_symbol, species, source) {
  csv_path <- prediction_cache_path(gene_symbol, species, source, "csv")
  rds_path <- prediction_cache_path(gene_symbol, species, source, "rds")
  file.exists(csv_path) || file.exists(rds_path)
}

bind_prediction_frames <- function(parts) {
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    stop("No prediction results to combine", call. = FALSE)
  }
  if (requireNamespace("dplyr", quietly = TRUE)) {
    return(dplyr::bind_rows(parts))
  }
  do.call(rbind, parts)
}

load_predictions <- function(
  gene_symbol,
  species,
  sources,
  use_cached_predictions,
  source_registry,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  unknown <- setdiff(sources, names(source_registry))
  if (length(unknown) > 0L) {
    stop("Unknown sources: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  parts <- lapply(sources, function(src) {
    if (isTRUE(use_cached_predictions) && prediction_cache_exists(gene_symbol, species, src)) {
      return(read_prediction_cache(gene_symbol, species, src))
    }
    source_registry[[src]](gene_symbol, species = species, ...)
  })
  bind_prediction_frames(parts)
}

export_prediction_source <- function(source, gene_symbol, species, source_registry, ...) {
  source <- validate_prediction_source(source)
  if (!source %in% names(source_registry)) {
    stop("Unknown source in registry: ", source, call. = FALSE)
  }
  df <- source_registry[[source]](gene_symbol, species = species, ...)
  paths <- write_prediction_cache(df, gene_symbol, species, source)
  list(predictions = df, paths = paths)
}

write_combined_predictions_cache <- function(df, gene_symbol, species) {
  df <- validate_prediction_frame(df, "combined predictions")
  ensure_data_dirs()
  path <- combined_predictions_cache_path(gene_symbol, species)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

write_integrated_cache <- function(df, gene_symbol, species) {
  if (!is.data.frame(df)) {
    stop("integrated result must be a data.frame", call. = FALSE)
  }
  ensure_data_dirs()
  csv_path <- integrated_cache_path(gene_symbol, species, "csv")
  rds_path <- integrated_cache_path(gene_symbol, species, "rds")
  write.csv(df, csv_path, row.names = FALSE)
  saveRDS(df, rds_path)
  invisible(list(csv = csv_path, rds = rds_path))
}
