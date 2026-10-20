project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

fixture_api <- normalizePath(
  file.path(project_root, "data", "cache", "fixtures", "mirdip", "api_pten_response.txt"),
  winslash = "/",
  mustWork = FALSE
)

use_mirdip_fixture <- function() {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  source("R/config.R")
  source("R/module_mirdip.R")
  options(
    mirna.mirdip.disable_api = TRUE,
    mirna.mirdip.api_response_path = fixture_api
  )
}

test_that("query_targets_mirdip parses API fixture for PTEN", {
  use_mirdip_fixture()
  result <- query_targets_mirdip("PTEN", species = "hsa", score_class = "Low")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_true(all(result$source == "mirdip"))
  expect_equal(result$mirna[[1]], "hsa-miR-21-5p")
  expect_equal(result$score_class[[1]], "Very High")
})

test_that("mirDIP rejects non-human species", {
  use_mirdip_fixture()
  expect_error(query_targets_mirdip("PTEN", species = "mmu"), "human")
})
