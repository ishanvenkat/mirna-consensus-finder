project_root <- normalizePath(
  file.path(getwd(), "..", "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "config.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

fixture_root <- normalizePath(
  file.path(project_root, "data", "cache", "fixtures", "targetscan"),
  winslash = "/",
  mustWork = FALSE
)

source_project <- function(relpath) {
  source(file.path(project_root, relpath))
}

clear_targetscan_cache <- function() {
  if (exists(".targetscan_cache", envir = .GlobalEnv)) {
    rm(list = ls(envir = .targetscan_cache), envir = .targetscan_cache)
  }
}

use_targetscan_fixture <- function() {
  withr::local_dir(project_root, .local_envir = parent.frame())
  source("R/config.R")
  source("R/module_targetscan.R")
  clear_targetscan_cache()
  withr::local_options(
    list(mirna.targetscan.data_root = fixture_root),
    .local_envir = parent.frame()
  )
  withr::defer(clear_targetscan_cache(), envir = parent.frame())
}

with_project_wd <- function(expr) {
  old_wd <- getwd()
  setwd(project_root)
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

test_that("query_targets_targetscan returns expected columns for PTEN", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("PTEN", species = "hsa")

  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c(
      "gene_symbol", "mirna", "score_raw", "score_name", "score_direction",
      "site_count", "source", "retrieved_at", "aggregate_pct",
      "mirna_family", "targetscan_representative_mirna"
    )
  )
  expect_equal(nrow(result), 8L)
  expect_true(all(result$score_direction == "lower_is_stronger"))
  expect_true(all(result$source == "targetscan"))
  expect_equal(result$mirna[[1]], "hsa-let-7a-5p")
  expect_equal(result$score_raw[[1]], -1.85)
  expect_lte(result$score_raw[[1]], result$score_raw[[2]])
})

test_that("query_targets_targetscan filters to representative transcript only", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("PTEN", species = "hsa")
  expect_equal(nrow(result), 8L)
})

test_that("unknown gene returns empty tibble", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("ZZZGENE", species = "hsa")
  expect_equal(nrow(result), 0L)
})

test_that("missing bulk data throws actionable error", {
  with_project_wd({
    source("R/config.R")
    source("R/module_targetscan.R")
    clear_targetscan_cache()
    old_root <- getOption("mirna.targetscan.data_root")
    options(mirna.targetscan.data_root = file.path(tempdir(), "nonexistent_targetscan"))
    on.exit({
      options(mirna.targetscan.data_root = old_root)
      clear_targetscan_cache()
    }, add = TRUE)

    expect_error(
      query_targets_targetscan("PTEN", species = "hsa"),
      "Run: Rscript scripts/download_targetscan.R"
    )
  })
})

test_that("validate_gene_symbol rejects bad input before file access", {
  with_project_wd({
    source("R/utils_gene_id.R")
    source("R/module_targetscan.R")
    expect_error(query_targets_targetscan(""), "Invalid gene symbol")
  })
})

test_that("predictions = default uses conservative summary file", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("PTEN", species = "hsa", predictions = "default")
  expect_equal(nrow(result), 2L)
  expect_true(all(c("hsa-let-7a-5p", "hsa-let-7b-5p") %in% result$mirna))
})

test_that("score_name biochem_occupancy selects occupancy column", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("PTEN", species = "hsa", score_name = "biochem_occupancy")
  expect_true(all(result$score_name == "biochem_occupancy"))
  let7 <- result[result$mirna == "hsa-let-7a-5p", , drop = FALSE]
  expect_equal(let7$score_raw[[1]], 0.12)
  expect_equal(result$mirna[[1]], "hsa-miR-494-3p")
  expect_equal(result$score_raw[[1]], 1.73)
})

test_that("mirna_prefix filters expanded mature miRNA names", {
  use_targetscan_fixture()
  all_rows <- query_targets_targetscan("PTEN", species = "hsa")
  expect_equal(nrow(all_rows), 8L)

  hsa_only <- query_targets_targetscan("PTEN", species = "hsa", mirna_prefix = "hsa")
  expect_equal(nrow(hsa_only), 8L)
  expect_true(all(grepl("^hsa-", hsa_only$mirna)))

  species_prefix <- query_targets_targetscan("PTEN", species = "hsa", mirna_prefix = "species")
  expect_equal(nrow(species_prefix), 8L)

  mmu_only <- query_targets_targetscan("PTEN", species = "hsa", mirna_prefix = "mmu")
  expect_equal(nrow(mmu_only), 0L)
})

test_that("mirna and mirna_family filter summary rows before expansion", {
  use_targetscan_fixture()
  one_mirna <- query_targets_targetscan("PTEN", mirna = "hsa-miR-21-5p")
  expect_equal(nrow(one_mirna), 1L)
  expect_equal(one_mirna$score_raw[[1]], -1.42)

  family <- query_targets_targetscan("PTEN", mirna_family = "let-7/98")
  expect_equal(nrow(family), 2L)
  expect_setequal(family$mirna, c("hsa-let-7a-5p", "hsa-let-7b-5p"))
})

test_that("family expansion emits all seed-family members for miR-17 cluster", {
  use_targetscan_fixture()
  result <- query_targets_targetscan("PTEN", species = "hsa")
  cluster <- c("hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p")
  expect_true(all(cluster %in% result$mirna))
  cluster_rows <- result[result$mirna %in% cluster, , drop = FALSE]
  expect_equal(length(unique(cluster_rows$score_raw)), 1L)
  expect_equal(unique(cluster_rows$targetscan_representative_mirna), "hsa-miR-17-5p")
})

test_that("min_cwcs retains stronger (more negative) CWCS rows", {
  use_targetscan_fixture()
  filtered <- query_targets_targetscan("PTEN", min_cwcs = -0.7)
  expect_equal(nrow(filtered), 4L)
  expect_true(all(filtered$score_raw <= -0.7))
})

test_that("hsa-miR-494-3p outranks hsa-miR-145-5p for PTEN by CWCS", {
  use_targetscan_fixture()
  # CWCS anchors from TargetScan 8.0 PTEN bulk (-0.799 vs -0.002); ranking regression.
  result <- query_targets_targetscan("PTEN", species = "hsa", mirna_prefix = "species")
  expect_true(all(c("hsa-miR-494-3p", "hsa-miR-145-5p") %in% result$mirna))
  i494 <- match("hsa-miR-494-3p", result$mirna)
  i145 <- match("hsa-miR-145-5p", result$mirna)
  expect_lt(i494, i145)
  expect_lt(result$score_raw[[i494]], result$score_raw[[i145]])
})

test_that("CWCS uses lower_is_stronger and biochem uses higher_is_stronger", {
  use_targetscan_fixture()
  cwcs <- query_targets_targetscan("PTEN", score_name = "cwcs")
  biochem <- query_targets_targetscan("PTEN", score_name = "biochem_occupancy")
  expect_equal(unique(cwcs$score_direction), "lower_is_stronger")
  expect_equal(unique(biochem$score_direction), "higher_is_stronger")
})
