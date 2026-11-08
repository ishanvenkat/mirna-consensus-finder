project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

make_integrated <- function() {
  data.frame(
    gene_symbol = c("PTEN", "PTEN", "PTEN"),
    mirna = c("hsa-miR-21-5p", "hsa-miR-100-5p", "hsa-miR-1-3p"),
    consensus_score = c(3.5, 2.0, 1.0),
    consensus_rank = c(1L, 2L, 3L),
    n_sources = c(3L, 2L, 1L),
    integration_method = "tiered",
    sources_hit = c("targetscan;diana", "targetscan", "mirdb"),
    consensus_tier = c("multi_source", "multi_source", "single_source"),
  stringsAsFactors = FALSE
  )
}

make_validation <- function() {
  data.frame(
    gene_symbol = c("PTEN", "PTEN"),
    mirna = c("hsa-miR-21-5p", "hsa-miR-214-3p"),
    experimental_support = c(TRUE, TRUE),
    support_level = c("strong", "strong"),
    validation_methods = c("Luciferase reporter assay", "Western blot"),
    validation_sources = c("miRTarBase", "miRTarBase"),
    binding_site_reported = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

test_that("annotate_integrated_with_validation leaves ranks and scores unchanged", {
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  integrated <- make_integrated()
  validation <- make_validation()
  out <- annotate_integrated_with_validation(integrated, validation)

  expect_equal(out$consensus_rank, integrated$consensus_rank)
  expect_equal(out$consensus_score, integrated$consensus_score)
  expect_true(out$experimental_support[[1]])
  expect_false(out$experimental_support[[2]])
  expect_equal(out$row_type, c("predicted", "predicted", "predicted"))
})

test_that("append_validated_only adds validation-only miRNAs", {
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  integrated <- make_integrated()
  validation <- make_validation()
  out <- annotate_integrated_with_validation(
    integrated,
    validation,
    append_validated_only = TRUE
  )

  expect_equal(nrow(out), 4L)
  extra <- out[out$row_type == "validated_only", , drop = FALSE]
  expect_equal(nrow(extra), 1L)
  expect_equal(extra$mirna[[1]], "hsa-miR-214-3p")
  expect_true(extra$experimental_support[[1]])
  expect_true(is.na(extra$consensus_rank[[1]]))
})

test_that("empty validation marks experimental_support FALSE", {
  setwd(project_root)
  source("R/config.R")
  source("R/module_validation.R")

  integrated <- make_integrated()
  validation <- empty_validation_result()
  out <- annotate_integrated_with_validation(integrated, validation)
  expect_false(any(out$experimental_support))
})
