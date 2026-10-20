project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

fixture_root <- normalizePath(
  file.path(project_root, "data", "cache", "fixtures", "diana"),
  winslash = "/",
  mustWork = FALSE
)
fixture_ensembl_cache <- file.path(fixture_root, "ensembl_symbol_map.tsv")
fixture_api_json <- file.path(fixture_root, "api_pten_response.json")

use_diana_fixture <- function() {
  old_wd <- getwd()
  setwd(project_root)
  if (requireNamespace("withr", quietly = TRUE)) {
    withr::defer(setwd(old_wd), envir = parent.frame())
  } else {
    on.exit(setwd(old_wd), add = TRUE, envir = parent.frame())
  }
  source("R/config.R")
  source("R/module_diana.R")
  options(
    mirna.diana.data_root = fixture_root,
    mirna.ensembl.cache_path = fixture_ensembl_cache,
    mirna.diana.disable_api = TRUE
  )
}

with_project_wd <- function(expr) {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

test_that("query_targets_diana bulk fallback returns expected columns for PTEN", {
  use_diana_fixture()
  result <- query_targets_diana("PTEN", species = "hsa", threshold = 0.7)

  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c(
      "gene_symbol", "mirna", "score_raw", "score_name", "score_direction",
      "site_count", "source", "retrieved_at"
    )
  )
  expect_equal(nrow(result), 3L)
  expect_true(all(result$score_name == "microT"))
  expect_true(all(result$source == "diana"))
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
  expect_true(all(result$score_raw >= 0.7))
})

test_that("query_targets_diana applies threshold on bulk data", {
  use_diana_fixture()
  result <- query_targets_diana("PTEN", species = "hsa", threshold = 0.8)
  expect_equal(nrow(result), 1L)
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
})

test_that("parse API fixture via mock response path", {
  use_diana_fixture()
  options(mirna.diana.disable_api = FALSE, mirna.diana.api_response_path = fixture_api_json)
  on.exit(options(mirna.diana.disable_api = TRUE, mirna.diana.api_response_path = NULL), add = TRUE)

  result <- query_targets_diana("PTEN", species = "hsa", threshold = 0.7)
  expect_equal(nrow(result), 2L)
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
})

test_that("missing bulk data throws actionable error", {
  with_project_wd({
    source("R/config.R")
    source("R/module_diana.R")
    options(
      mirna.diana.data_root = file.path(tempdir(), "missing_diana"),
      mirna.ensembl.cache_path = fixture_ensembl_cache,
      mirna.diana.disable_api = TRUE
    )
    on.exit(options(
      mirna.diana.data_root = NULL,
      mirna.ensembl.cache_path = NULL,
      mirna.diana.disable_api = NULL
    ), add = TRUE)

    expect_error(
      query_targets_diana("PTEN", species = "hsa"),
      "Run: Rscript scripts/download_diana_bulk.R"
    )
  })
})

test_that("symbol_to_ensembl reads fixture cache", {
  with_project_wd({
    source("R/config.R")
    source("R/utils_gene_id.R")
    options(mirna.ensembl.cache_path = fixture_ensembl_cache)
    on.exit(options(mirna.ensembl.cache_path = NULL), add = TRUE)

    expect_equal(symbol_to_ensembl("PTEN", "hsa"), "ENSG00000171862")
  })
})
