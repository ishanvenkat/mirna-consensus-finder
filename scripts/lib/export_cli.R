# Shared helpers for export_score_* CLI scripts (sourced, not run directly).

if (!file.exists("R/config.R") && file.exists("../R/config.R")) {
  setwd(normalizePath("..", winslash = "/", mustWork = FALSE))
}

resolve_project_root <- function() {
  if (file.exists("R/config.R")) {
    return(normalizePath(".", winslash = "/", mustWork = FALSE))
  }
  if (file.exists("../R/config.R")) {
    return(normalizePath("..", winslash = "/", mustWork = FALSE))
  }
  stop("Run from project root: Rscript scripts/export_score_<source>.R", call. = FALSE)
}

parse_flag_args <- function(args) {
  flags <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) {
      next
    }
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      parts <- strsplit(kv, "=", fixed = TRUE)[[1L]]
      key <- parts[[1L]]
      val <- paste(parts[-1L], collapse = "=")
    } else {
      key <- kv
      val <- "TRUE"
    }
    flags[[key]] <- val
  }
  flags
}

positional_gene_species <- function(args) {
  pos <- args[!grepl("^--", args)]
  if (length(pos) < 1L) {
    stop("GENE argument is required", call. = FALSE)
  }
  gene <- pos[[1L]]
  species <- if (length(pos) >= 2L) pos[[2L]] else "hsa"
  flag_args <- args[grepl("^--", args)]
  list(gene = gene, species = species, flags = parse_flag_args(flag_args))
}

flag_bool <- function(flags, key, default = FALSE) {
  if (is.null(flags[[key]])) {
    return(default)
  }
  val <- tolower(as.character(flags[[key]]))
  val %in% c("1", "true", "yes", "y")
}

flag_numeric <- function(flags, key, default = NULL) {
  if (is.null(flags[[key]])) {
    return(default)
  }
  as.numeric(flags[[key]])
}

flag_character <- function(flags, key, default = NULL) {
  if (is.null(flags[[key]])) {
    return(default)
  }
  as.character(flags[[key]])
}

setup_export_project <- function() {
  root <- resolve_project_root()
  setwd(root)
  source("R/main.R")
  root
}
