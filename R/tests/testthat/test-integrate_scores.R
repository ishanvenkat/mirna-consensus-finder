project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

test_that("integrate_scores rank_fusion aggregates independent sources", {
  setwd(project_root)
  source("R/integrate_scores.R")

  predictions <- data.frame(
    gene_symbol = rep("PTEN", 5L),
    mirna = c(
      "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
      "hsa-miR-100-3p", "hsa-miR-21-5p"
    ),
    score_raw = c(-1.8, 0.85, 85, 0.5, 0.91),
    score_direction = c(
      "lower_is_stronger", "higher_is_stronger", "higher_is_stronger",
      "higher_is_stronger", "higher_is_stronger"
    ),
    source = c("targetscan", "diana", "mirdb", "diana", "mirdip"),
    score_class = c(NA, NA, NA, NA, "Very High"),
    stringsAsFactors = FALSE
  )

  fused <- integrate_scores(predictions, method = "rank_fusion")
  expect_equal(fused$gene_symbol[[1]], "PTEN")
  expect_true("hsa-miR-21-5p" %in% fused$mirna)
  expect_equal(fused$mirna[[1]], "hsa-miR-21-5p")
  expect_equal(fused$mirdip_score_class[[1]], "Very High")
})

test_that("integrate_scores tiered excludes mirDIP from fusion sources", {
  setwd(project_root)
  source("R/integrate_scores.R")

  predictions <- data.frame(
    gene_symbol = c("PTEN", "PTEN"),
    mirna = c("hsa-miR-21-5p", "hsa-miR-21-5p"),
    score_raw = c(10, 0.99),
    score_direction = rep("higher_is_stronger", 2L),
    source = c("mirdb", "mirdip"),
    stringsAsFactors = FALSE
  )

  fused <- integrate_scores(predictions, method = "tiered")
  expect_equal(nrow(fused), 1L)
  expect_equal(fused$n_sources, 1L)
})

test_that("integrate_scores ranks TargetScan CWCS with lower_is_stronger", {
  setwd(project_root)
  source("R/integrate_scores.R")

  predictions <- data.frame(
    gene_symbol = rep("PTEN", 2L),
    mirna = c("hsa-miR-21-5p", "hsa-miR-100-3p"),
    score_raw = c(-0.002, -0.799),
    score_direction = rep("lower_is_stronger", 2L),
    source = rep("targetscan", 2L),
    stringsAsFactors = FALSE
  )

  fused <- integrate_scores(predictions, method = "rank_fusion")
  expect_equal(fused$mirna[[1]], "hsa-miR-100-3p")
})
