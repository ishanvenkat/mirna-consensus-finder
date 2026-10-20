project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

fixture_path <- normalizePath(
  file.path(project_root, "data", "cache", "fixtures", "validation", "pten_validated.json"),
  winslash = "/",
  mustWork = FALSE
)

test_that("query_validation reads fixture overlay for PTEN", {
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")
  options(mirna.validation.fixture_path = fixture_path)
  on.exit(options(mirna.validation.fixture_path = NULL), add = TRUE)

  result <- query_validation("PTEN", species = "hsa")
  expect_equal(nrow(result), 2L)
  expect_true(all(result$experimental_support))
  expect_equal(result$support_level[[1]], "strong")
})
