# Module Interface Reference

## Example Return Tibbles

### TargetScan

```r
# query_targets_targetscan(..., mirna_prefix = NULL, mirna = NULL, mirna_family = NULL)
# Summary_Counts scores are family-level; miR_Family_Info.txt expands to all mature members.
# mirna_prefix: filters expanded mature names (e.g. "hsa" or "species")
# score_name: "cwcs" (default) or "biochem_occupancy"
# score_direction: "lower_is_stronger" for CWCS; "higher_is_stronger" for biochem occupancy
# aggregate_pct: conservation annotation — not merged into score_raw
tibble(
  gene_symbol = "PTEN",
  mirna = "hsa-miR-20a-5p",
  score_raw = -0.65,
  score_name = "cwcs",
  score_direction = "lower_is_stronger",
  site_count = 1L,
  source = "targetscan",
  retrieved_at = Sys.time(),
  aggregate_pct = 0.10,
  mirna_family = "miR-17/20/93/106/519/526/1283",
  targetscan_representative_mirna = "hsa-miR-17-5p"
)
```

### DIANA-microT-CDS

```r
# score_name: "microT"
# score_direction: "higher_is_stronger" (range ~0.3–1.0)
tibble(
  gene_symbol = "PTEN",
  mirna = "hsa-miR-21-5p",
  score_raw = 0.82,
  score_name = "microT",
  score_direction = "higher_is_stronger",
  site_count = NA_integer_,
  source = "diana",
  retrieved_at = Sys.time()
)
```

### mirDIP

```r
# score_name: "integrated"
tibble(
  gene_symbol = "PTEN",
  mirna = "hsa-miR-21-5p",
  score_raw = 0.91,
  score_name = "integrated",
  score_direction = "higher_is_stronger",
  site_count = NA_integer_,
  source = "mirdip",
  retrieved_at = Sys.time()
)
```

## Prediction cache (data/cache/)

| Artifact | Path |
|----------|------|
| Per source | `<GENE>_<species>_<source>.csv` and `.rds` |
| Combined | `<GENE>_<species>_predictions_combined.csv` |
| Integrated (full) | `<GENE>_<species>_integrated.csv` and `.rds` |
| Integrated consensus | `<GENE>_<species>_integrated_consensus.csv` and `.rds` (`n_sources >= 2`) |
| Integrated single-source | `<GENE>_<species>_integrated_single_source.csv` and `.rds` (`n_sources == 1`) |

Helpers live in `R/cache_predictions.R`: `write_prediction_cache`, `read_prediction_cache`, `load_predictions`, `annotate_integrated`, `write_integrated_cache`.

CLI: `scripts/export_score_{targetscan,diana,mirdb,mirdip}.R`, `scripts/export_score_integrated.R` (shared `scripts/lib/export_cli.R`).

## Main Entry Point

```r
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
  # 1. Validate gene_symbol
  # 2. load_predictions() per source (live, or cache + live fallback)
  # 3. Optional integrate_scores()
  # 4. Optional query_validation()
  # 5. Return list(predictions, integrated, validation)
}
```

## Merge Pattern

```r
predictions <- purrr::map_dfr(sources, function(src) {
  query_fn <- get(paste0("query_targets_", src), envir = asNamespace("mirnaTargets"))
  query_fn(gene_symbol, species = species)
})
```

Each module file should define `query_targets_<source>` and re-export as `query_targets` when sourced individually, or register in a module list in `main.R`.
