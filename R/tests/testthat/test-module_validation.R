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

snippet_path <- normalizePath(
  file.path(
    project_root,
    "data",
    "cache",
    "fixtures",
    "validation",
    "mirtarbase_pten_snippet.csv"
  ),
  winslash = "/",
  mustWork = FALSE
)

test_that("is_strong_mirtarbase_method detects reporter and Western", {
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  expect_true(is_strong_mirtarbase_method("Luciferase reporter assay"))
  expect_true(is_strong_mirtarbase_method("Western blot"))
  expect_false(is_strong_mirtarbase_method("qRT-PCR"))
  expect_false(is_strong_mirtarbase_method("PAR-CLIP"))
})

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

test_that("mirtarbase local filter keeps reporter/Western PTEN rows", {
  skip_if_not(file.exists(snippet_path))
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  options(mirna.mirtarbase.bulk_path = snippet_path)
  on.exit(options(mirna.mirtarbase.bulk_path = NULL), add = TRUE)

  rows <- query_validation_mirtarbase_local(
    "PTEN",
    species = "hsa",
    mirtarbase_strong_methods = "reporter_western"
  )
  mirnas <- vapply(rows, function(r) r$mirna, character(1L))
  expect_true("hsa-miR-21-5p" %in% mirnas)
  expect_true("hsa-miR-214-3p" %in% mirnas)
  expect_true("hsa-miR-93-5p" %in% mirnas)
  expect_false("hsa-miR-19b-3p" %in% mirnas)
  expect_false("hsa-miR-221-3p" %in% mirnas)
})

test_that("query_validation_tarbase_multimir returns PTEN rows when multiMiR installed", {
  skip_if_not(requireNamespace("multiMiR", quietly = TRUE))
  skip_if(isTRUE(getOption("mirna.validation.disable_multimir")))

  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  rows <- query_validation_tarbase_multimir("PTEN", "hsa")
  expect_type(rows, "list")
  if (length(rows) > 0L) {
    expect_equal(toupper(unique(vapply(rows, function(r) r$sources, character(1L)))), "TARBASE")
  }
})
