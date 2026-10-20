project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source_utils <- function() {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  source("R/config.R")
  source("R/utils_gene_id.R")
}

test_that("validate_gene_symbol accepts PTEN", {
  source_utils()
  expect_equal(validate_gene_symbol("pten"), "PTEN")
})

test_that("validate_gene_symbol rejects invalid input", {
  source_utils()
  expect_error(validate_gene_symbol(""), "Invalid gene symbol")
  expect_error(validate_gene_symbol(c("A", "B")), "single character")
})
