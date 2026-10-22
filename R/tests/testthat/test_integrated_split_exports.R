project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source_project <- function() {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  source("R/config.R")
  source("R/utils_gene_id.R")
  source("R/integrate_scores.R")
  source("R/cache_predictions.R")
}

with_temp_cache <- function(code) {
  source_project()
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

test_that("annotate_integrated sets sources_hit and consensus_tier", {
  source_project()

  predictions <- data.frame(
    gene_symbol = rep("PTEN", 4L),
    mirna = c(
      "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
      "hsa-miR-100-3p"
    ),
    score_raw = c(-1.0, 0.9, 85, 0.5),
    score_direction = c(
      "lower_is_stronger", "higher_is_stronger", "higher_is_stronger",
      "higher_is_stronger"
    ),
    source = c("targetscan", "diana", "mirdb", "diana"),
    stringsAsFactors = FALSE
  )

  integrated <- integrate_scores(predictions, method = "tiered")
  out <- annotate_integrated(
    integrated,
    predictions,
    fusion_sources = c("targetscan", "diana", "mirdb")
  )

  hit21 <- out[out$mirna == "hsa-miR-21-5p", , drop = FALSE]
  expect_equal(hit21$sources_hit[[1L]], "targetscan,diana,mirdb")
  expect_equal(hit21$consensus_tier[[1L]], "multi_source")
  expect_equal(hit21$n_sources, 3L)

  hit100 <- out[out$mirna == "hsa-miR-100-3p", , drop = FALSE]
  expect_equal(hit100$sources_hit[[1L]], "diana")
  expect_equal(hit100$consensus_tier[[1L]], "single_source")
  expect_equal(hit100$n_sources, 1L)
})

test_that("write_integrated_cache writes split primary and secondary files", {
  with_temp_cache({
    predictions <- data.frame(
      gene_symbol = rep("PTEN", 4L),
      mirna = c(
        "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
        "hsa-miR-100-3p"
      ),
      score_raw = c(-1.0, 0.9, 85, 0.5),
      score_direction = c(
        "lower_is_stronger", "higher_is_stronger", "higher_is_stronger",
        "higher_is_stronger"
      ),
      source = c("targetscan", "diana", "mirdb", "diana"),
      stringsAsFactors = FALSE
    )

    integrated <- integrate_scores(predictions, method = "tiered")
    integrated <- annotate_integrated(
      integrated,
      predictions,
      fusion_sources = c("targetscan", "diana", "mirdb")
    )
    paths <- write_integrated_cache(integrated, "PTEN", "hsa")

    expect_true(file.exists(paths$csv))
    expect_true(file.exists(paths$consensus_csv))
    expect_true(file.exists(paths$single_source_csv))

    full <- read.csv(paths$csv, stringsAsFactors = FALSE)
    multi <- read.csv(paths$consensus_csv, stringsAsFactors = FALSE)
    single <- read.csv(paths$single_source_csv, stringsAsFactors = FALSE)

    expect_equal(nrow(full), 2L)
    expect_equal(nrow(multi), 1L)
    expect_equal(nrow(single), 1L)
    expect_equal(multi$mirna[[1L]], "hsa-miR-21-5p")
    expect_equal(multi$consensus_rank[[1L]], 1L)
    expect_equal(single$mirna[[1L]], "hsa-miR-100-3p")
    expect_equal(single$single_source_rank[[1L]], 1L)
    expect_false("hsa-miR-100-3p" %in% multi$mirna)
  })
})
