project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source_eval <- function() {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  source("R/evaluate_consensus.R")
}

test_that("validation_positives filters strong experimental support", {
  source_eval()
  validation <- data.frame(
    mirna = c("hsa-miR-21-5p", "hsa-miR-214-3p", "hsa-miR-19b-3p"),
    experimental_support = c(TRUE, TRUE, TRUE),
    support_level = c("strong", "strong", "weak"),
    stringsAsFactors = FALSE
  )
  strong <- validation_positives(validation, support_level = "strong")
  expect_equal(sort(strong$mirnas), c("hsa-miR-21-5p", "hsa-miR-214-3p"))
  any_level <- validation_positives(validation, support_level = "any")
  expect_equal(length(any_level$mirnas), 3L)
})

test_that("recall_at_k matches hits in top k", {
  source_eval()
  positives <- c("hsa-miR-21-5p", "hsa-miR-214-3p")
  ranked <- c(
    "hsa-miR-1-3p", "hsa-miR-2-3p", "hsa-miR-3-3p", "hsa-miR-4-3p",
    "hsa-miR-21-5p", "hsa-miR-6-3p", "hsa-miR-7-3p", "hsa-miR-214-3p"
  )
  tbl <- recall_at_k_table(ranked, positives, c(10L), "test_method")
  expect_equal(tbl$hits[[1L]], 2L)
  expect_equal(tbl$recall[[1L]], 1)
  expect_equal(tbl$k_effective[[1L]], 8L)
  expect_equal(tbl$precision[[1L]], 2 / 8)

  ranked_one <- ranked
  ranked_one[[5L]] <- "hsa-miR-5-3p"
  ranked_one[[8L]] <- "hsa-miR-214-3p"
  tbl_one <- recall_at_k_table(ranked_one, positives, c(10L), "test_method")
  expect_equal(tbl_one$hits[[1L]], 1L)
  expect_equal(tbl_one$recall[[1L]], 0.5)
})

test_that("pr_curve_from_ranking at depth 10 matches cumulative hits", {
  source_eval()
  positives <- c("hsa-miR-21-5p", "hsa-miR-214-3p")
  ranked <- c(
    "hsa-miR-1-3p", "hsa-miR-2-3p", "hsa-miR-3-3p", "hsa-miR-4-3p",
    "hsa-miR-21-5p", "hsa-miR-6-3p", "hsa-miR-7-3p", "hsa-miR-214-3p",
    "hsa-miR-9-3p", "hsa-miR-10-3p"
  )
  pr <- pr_curve_from_ranking(ranked, positives, max_depth = 10L)
  row10 <- pr[pr$depth == 10L, , drop = FALSE]
  expect_equal(row10$hits_cumulative[[1L]], 2L)
  expect_equal(row10$precision[[1L]], 0.2)
  expect_equal(row10$recall[[1L]], 1)
})

test_that("rank_predictions_by_source respects lower_is_stronger", {
  source_eval()
  predictions <- data.frame(
    gene_symbol = rep("PTEN", 2L),
    mirna = c("hsa-miR-21-5p", "hsa-miR-100-3p"),
    score_raw = c(-0.002, -0.799),
    score_direction = rep("lower_is_stronger", 2L),
    source = rep("targetscan", 2L),
    stringsAsFactors = FALSE
  )
  ranked <- rank_predictions_by_source(predictions, "targetscan")
  expect_equal(ranked[[1L]], "hsa-miR-100-3p")
  expect_equal(ranked[[2L]], "hsa-miR-21-5p")
})

test_that("ranked_list_from_integrated orders by consensus_rank", {
  source_eval()
  integrated <- data.frame(
    mirna = c("hsa-miR-b", "hsa-miR-a"),
    consensus_rank = c(2L, 1L),
    consensus_score = c(10, 5),
    n_sources = c(2L, 1L),
    stringsAsFactors = FALSE
  )
  full <- ranked_list_from_integrated(integrated, subset = "full")
  expect_equal(full, c("hsa-miR-a", "hsa-miR-b"))

  multi <- ranked_list_from_integrated(integrated, subset = "multi_source")
  expect_equal(multi, "hsa-miR-b")
})