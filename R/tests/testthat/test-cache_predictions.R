project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source_cache <- function() {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  source("R/config.R")
  source("R/utils_gene_id.R")
  source("R/cache_predictions.R")
}

sample_prediction_frame <- function(source = "targetscan") {
  data.frame(
    gene_symbol = "PTEN",
    mirna = "hsa-miR-21-5p",
    score_raw = if (source == "targetscan") -1.0 else 0.9,
    score_direction = if (source == "targetscan") {
      "lower_is_stronger"
    } else {
      "higher_is_stronger"
    },
    source = source,
    stringsAsFactors = FALSE
  )
}

with_temp_cache <- function(code) {
  source_cache()
  tmp_cache <- tempfile("mirna_cache_")
  dir.create(tmp_cache)
  old_cache <- PATHS$cache
  PATHS$cache <<- tmp_cache
  on.exit({
    PATHS$cache <<- old_cache
    unlink(tmp_cache, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

test_that("prediction_cache_path uses GENE_species_source naming", {
  source_cache()
  expect_equal(
    basename(prediction_cache_path("PTEN", "hsa", "targetscan", "csv")),
    "PTEN_hsa_targetscan.csv"
  )
  expect_equal(
    basename(integrated_cache_path("PTEN", "hsa", "rds")),
    "PTEN_hsa_integrated.rds"
  )
  expect_equal(
    basename(combined_predictions_cache_path("PTEN", "hsa")),
    "PTEN_hsa_predictions_combined.csv"
  )
})

test_that("write_prediction_cache round-trips CSV and RDS", {
  with_temp_cache({
    df <- sample_prediction_frame("diana")
    paths <- write_prediction_cache(df, "PTEN", "hsa", "diana")
    expect_true(file.exists(paths$csv))
    expect_true(file.exists(paths$rds))
    back <- read_prediction_cache("PTEN", "hsa", "diana")
    expect_equal(back$mirna, df$mirna)
    expect_equal(back$source, "diana")
  })
})

test_that("load_predictions uses cache when present and live fallback otherwise", {
  with_temp_cache({
    write_prediction_cache(
      sample_prediction_frame("targetscan"),
      "PTEN",
      "hsa",
      "targetscan"
    )
    live_called <- FALSE
    registry <- list(
      targetscan = function(g, ...) {
        stop("should not call live targetscan")
      },
      diana = function(g, ...) {
        live_called <<- TRUE
        sample_prediction_frame("diana")
      }
    )
    combined <- load_predictions(
      "PTEN",
      "hsa",
      sources = c("targetscan", "diana"),
      use_cached_predictions = TRUE,
      source_registry = registry
    )
    expect_true(live_called)
    expect_equal(nrow(combined), 2L)
    expect_setequal(combined$source, c("targetscan", "diana"))
  })
})

test_that("load_predictions calls live when use_cached_predictions is FALSE", {
  with_temp_cache({
    write_prediction_cache(
      sample_prediction_frame("targetscan"),
      "PTEN",
      "hsa",
      "targetscan"
    )
    live_ts <- FALSE
    registry <- list(
      targetscan = function(g, ...) {
        live_ts <<- TRUE
        sample_prediction_frame("targetscan")
      }
    )
    load_predictions(
      "PTEN",
      "hsa",
      sources = "targetscan",
      use_cached_predictions = FALSE,
      source_registry = registry
    )
    expect_true(live_ts)
  })
})
