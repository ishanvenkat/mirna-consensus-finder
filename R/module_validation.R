# Experimental validation overlay — miRTarBase and DIANA-TarBase.
# Read-only reference layer; never mixed into prediction scores.

source("R/config.R")
source("R/utils_gene_id.R")

empty_validation_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    experimental_support = logical(),
    support_level = character(),
    validation_methods = character(),
    validation_sources = character(),
    binding_site_reported = logical(),
    retrieved_at = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )
}

normalize_validation_rows <- function(gene_symbol, rows) {
  if (length(rows) == 0L) {
    return(empty_validation_result())
  }

  mirna <- vapply(rows, function(r) r$mirna, character(1L))
  support_level <- vapply(rows, function(r) r$support_level %||% NA_character_, character(1L))
  methods <- vapply(rows, function(r) r$methods %||% "", character(1L))
  sources <- vapply(rows, function(r) r$sources, character(1L))
  binding <- vapply(rows, function(r) isTRUE(r$binding_site_reported), logical(1L))

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = mirna,
    experimental_support = TRUE,
    support_level = support_level,
    validation_methods = methods,
    validation_sources = sources,
    binding_site_reported = binding,
    retrieved_at = Sys.time(),
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out[, c("gene_symbol", "mirna")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_validation_fixture <- function(path, gene_symbol) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read validation fixtures", call. = FALSE)
  }
  payload <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  rows <- list()
  for (item in payload) {
    if (toupper(item$gene_symbol) != gene_symbol) {
      next
    }
    rows[[length(rows) + 1L]] <- list(
      mirna = item$mirna,
      support_level = item$support_level,
      methods = paste(item$methods, collapse = "; "),
      sources = paste(item$sources, collapse = ";"),
      binding_site_reported = item$binding_site_reported
    )
  }
  rows
}

query_validation_multimir <- function(gene_symbol, species) {
  if (!requireNamespace("multiMiR", quietly = TRUE)) {
    stop(
      "Package 'multiMiR' is required for validation lookup. ",
      "Install with BiocManager::install('multiMiR').",
      call. = FALSE
    )
  }

  org <- species
  mm <- multiMiR::get_multimir(
    org = org,
    target = gene_symbol,
    table = "validated",
    mirtarbase = TRUE,
    tarbase = TRUE,
    summary = FALSE
  )

  rows <- list()
  if (!is.null(mm@data$mirtarbase) && nrow(mm@data$mirtarbase) > 0L) {
    mt <- mm@data$mirtarbase
    mirna_col <- if ("mature_mirna_id" %in% names(mt)) "mature_mirna_id" else "miRNA"
    target_col <- if ("target_symbol" %in% names(mt)) "target_symbol" else "target"
    support_col <- if ("support_type" %in% names(mt)) "support_type" else "support_type"
    method_col <- if ("experiments" %in% names(mt)) "experiments" else "experiment"

    for (i in seq_len(nrow(mt))) {
      if (toupper(trimws(mt[[target_col]][[i]])) != gene_symbol) {
        next
      }
      support_raw <- if (support_col %in% names(mt)) mt[[support_col]][[i]] else NA_character_
      support_level <- if (is.na(support_raw)) {
        NA_character_
      } else if (grepl("strong", support_raw, ignore.case = TRUE)) {
        "strong"
      } else {
        "weak"
      }
      rows[[length(rows) + 1L]] <- list(
        mirna = mt[[mirna_col]][[i]],
        support_level = support_level,
        methods = if (method_col %in% names(mt)) as.character(mt[[method_col]][[i]]) else "",
        sources = "miRTarBase",
        binding_site_reported = FALSE
      )
    }
  }

  if (!is.null(mm@data$tarbase) && nrow(mm@data$tarbase) > 0L) {
    tb <- mm@data$tarbase
    mirna_col <- if ("mature_mirna_id" %in% names(tb)) "mature_mirna_id" else "mirna"
    target_col <- if ("target_symbol" %in% names(tb)) "target_symbol" else "gene"
    method_col <- if ("experiment" %in% names(tb)) "experiment" else "method"
    site_col <- if ("binding_site" %in% names(tb)) "binding_site" else NULL

    for (i in seq_len(nrow(tb))) {
      if (toupper(trimws(tb[[target_col]][[i]])) != gene_symbol) {
        next
      }
      binding <- FALSE
      if (!is.null(site_col) && site_col %in% names(tb)) {
        binding <- nzchar(trimws(as.character(tb[[site_col]][[i]])))
      }
      rows[[length(rows) + 1L]] <- list(
        mirna = tb[[mirna_col]][[i]],
        support_level = NA_character_,
        methods = if (method_col %in% names(tb)) as.character(tb[[method_col]][[i]]) else "",
        sources = "TarBase",
        binding_site_reported = binding
      )
    }
  }

  rows
}

merge_validation_rows <- function(rows) {
  if (length(rows) == 0L) {
    return(list())
  }
  df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(
      mirna = r$mirna,
      support_level = r$support_level %||% NA_character_,
      methods = r$methods %||% "",
      sources = r$sources,
      binding_site_reported = isTRUE(r$binding_site_reported),
      stringsAsFactors = FALSE
    )
  }))

  split_rows <- split(df, df$mirna)
  lapply(split_rows, function(part) {
    levels <- part$support_level[!is.na(part$support_level)]
    support_level <- if (any(levels == "strong")) {
      "strong"
    } else if (length(levels) > 0L) {
      levels[[1L]]
    } else {
      NA_character_
    }
    list(
      mirna = part$mirna[[1L]],
      support_level = support_level,
      methods = paste(unique(part$methods[part$methods != ""]), collapse = "; "),
      sources = paste(unique(part$sources), collapse = ";"),
      binding_site_reported = any(part$binding_site_reported)
    )
  })
}

query_validation <- function(gene_symbol, species = "hsa", ...) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }

  fixture_path <- getOption("mirna.validation.fixture_path")
  rows <- if (!is.null(fixture_path) && nzchar(fixture_path)) {
    parse_validation_fixture(fixture_path, gene_symbol)
  } else if (isTRUE(getOption("mirna.validation.disable_multimir"))) {
    stop(
      "Validation lookup disabled (mirna.validation.disable_multimir = TRUE)",
      call. = FALSE
    )
  } else {
    query_validation_multimir(gene_symbol, species)
  }

  rows <- merge_validation_rows(rows)
  normalize_validation_rows(gene_symbol, rows)
}
