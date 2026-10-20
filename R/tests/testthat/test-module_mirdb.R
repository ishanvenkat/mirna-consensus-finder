project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

fixture_root <- normalizePath(
  file.path(project_root, "data", "cache", "fixtures", "mirdb"),
  winslash = "/",
  mustWork = FALSE
)
fixture_refseq_cache <- file.path(fixture_root, "refseq_symbol_map.tsv")

use_mirdb_fixture <- function() {
  old_wd <- getwd()
  setwd(project_root)
  if (requireNamespace("withr", quietly = TRUE)) {
    withr::defer(setwd(old_wd), envir = parent.frame())
  } else {
    on.exit(setwd(old_wd), add = TRUE, envir = parent.frame())
  }
  source("R/config.R")
  source("R/module_mirdb.R")
  options(
    mirna.mirdb.data_root = fixture_root,
    mirna.refseq.cache_path = fixture_refseq_cache
  )
}

with_project_wd <- function(expr) {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

test_that("query_targets_mirdb bulk returns expected columns for PTEN", {
  use_mirdb_fixture()
  result <- query_targets_mirdb("PTEN", species = "hsa", min_score = 60)

  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c(
      "gene_symbol", "mirna", "score_raw", "score_name", "score_direction",
      "site_count", "source", "retrieved_at"
    )
  )
  expect_equal(nrow(result), 3L)
  expect_true(all(result$score_name == "MirTarget"))
  expect_true(all(result$source == "mirdb"))
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
  expect_true(all(result$score_raw >= 60))
})

test_that("query_targets_mirdb applies min_score threshold", {
  use_mirdb_fixture()
  result <- query_targets_mirdb("PTEN", species = "hsa", min_score = 80)
  expect_equal(nrow(result), 1L)
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
})

test_that("missing bulk data throws actionable error", {
  with_project_wd({
    source("R/config.R")
    source("R/module_mirdb.R")
    options(
      mirna.mirdb.data_root = file.path(tempdir(), "missing_mirdb"),
      mirna.refseq.cache_path = fixture_refseq_cache
    )
    on.exit(options(
      mirna.mirdb.data_root = NULL,
      mirna.refseq.cache_path = NULL
    ), add = TRUE)

    expect_error(
      query_targets_mirdb("PTEN", species = "hsa"),
      "Run: Rscript scripts/download_mirdb.R"
    )
  })
})
