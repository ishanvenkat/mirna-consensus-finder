#!/usr/bin/env Rscript
# Smoke test: query PTEN across available sources.

args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) >= 1L) args[[1L]] else "PTEN"

source("R/main.R")

cat("Querying miRNA regulators for:", gene, "\n\n")

tryCatch(
  {
    result <- find_mirna_regulators(gene, sources = c("targetscan"))
    print(result)
  },
  error = function(e) {
    cat("Expected during scaffolding — modules not yet implemented:\n")
    cat(conditionMessage(e), "\n")
  }
)
