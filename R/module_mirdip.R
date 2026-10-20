# mirDIP integrative prediction module — Http_U API with response cache.

source("R/config.R")
source("R/utils_gene_id.R")

mirdip_api_url <- function() {
  paste0(API$mirdip_base, MIRDIP$api_path)
}

mirdip_cache_path <- function(gene_symbol, score_class) {
  safe_class <- gsub("[^A-Za-z0-9]+", "_", score_class)
  file.path(PATHS$cache, "mirdip", paste0(gene_symbol, "_", safe_class, ".txt"))
}

empty_mirdip_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    score_raw = numeric(),
    score_name = character(),
    score_direction = character(),
    site_count = integer(),
    source = character(),
    retrieved_at = as.POSIXct(character()),
    n_sources = integer(),
    score_class = character(),
    stringsAsFactors = FALSE
  )
}

validate_mirdip_score_class <- function(score_class) {
  if (!is.character(score_class) || length(score_class) != 1L || !nzchar(score_class)) {
    stop("score_class must be a single non-empty string", call. = FALSE)
  }
  if (!score_class %in% names(MIRDIP$score_class_codes)) {
    stop(
      "score_class must be one of: ",
      paste(names(MIRDIP$score_class_codes), collapse = ", "),
      call. = FALSE
    )
  }
  score_class
}

parse_mirdip_response <- function(body, gene_symbol) {
  body <- trimws(body)
  if (!nzchar(body)) {
    return(empty_mirdip_result())
  }

  parts <- strsplit(body, "\001", fixed = TRUE)[[1L]]
  fields <- list()
  for (part in parts) {
    kv <- strsplit(part, "\002", fixed = TRUE)[[1L]]
    if (length(kv) >= 2L && nzchar(kv[[1L]])) {
      fields[[kv[[1L]]]] <- paste(kv[-1L], collapse = "\002")
    }
  }

  if (is.null(fields$results)) {
    return(empty_mirdip_result())
  }

  results_text <- trimws(fields$results)
  if (!nzchar(results_text)) {
    return(empty_mirdip_result())
  }

  tg <- read.table(
    text = results_text,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
  if (nrow(tg) == 0L) {
    return(empty_mirdip_result())
  }

  gene_col <- if ("Gene Symbol" %in% names(tg)) {
    "Gene Symbol"
  } else if ("Gene.Symbol" %in% names(tg)) {
    "Gene.Symbol"
  } else {
    stop("mirDIP response missing gene symbol column", call. = FALSE)
  }
  mirna_col <- if ("MicroRNA" %in% names(tg)) "MicroRNA" else "microRNA"
  score_col <- if ("Integrated Score" %in% names(tg)) {
    "Integrated Score"
  } else {
    "Integrated.Score"
  }
  nsrc_col <- if ("Number of Sources" %in% names(tg)) {
    "Number of Sources"
  } else {
    "Number.of.Sources"
  }
  class_col <- if ("Score Class" %in% names(tg)) "Score Class" else "Score.Class"

  tg <- tg[toupper(trimws(tg[[gene_col]])) == gene_symbol, , drop = FALSE]
  if (nrow(tg) == 0L) {
    return(empty_mirdip_result())
  }

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = tg[[mirna_col]],
    score_raw = as.numeric(tg[[score_col]]),
    score_name = MIRDIP$score_name,
    score_direction = "higher_is_stronger",
    site_count = NA_integer_,
    source = "mirdip",
    retrieved_at = Sys.time(),
    n_sources = as.integer(tg[[nsrc_col]]),
    score_class = tg[[class_col]],
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$score_raw) & nzchar(out$mirna), , drop = FALSE]
  ord <- order(out$score_raw, out$mirna, decreasing = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

write_mirdip_cache <- function(path, body) {
  cache_dir <- dirname(path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  writeLines(body, path, useBytes = TRUE)
  invisible(path)
}

fetch_mirdip_api <- function(gene_symbol, score_class) {
  mock_path <- getOption("mirna.mirdip.api_response_path")
  if (!is.null(mock_path) && nzchar(mock_path)) {
    if (!file.exists(mock_path)) {
      stop("mirDIP mock response not found: ", mock_path, call. = FALSE)
    }
    return(paste(readLines(mock_path, warn = FALSE), collapse = "\n"))
  }

  if (isTRUE(getOption("mirna.mirdip.disable_api"))) {
    stop("mirDIP API disabled (mirna.mirdip.disable_api = TRUE)", call. = FALSE)
  }

  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop(
      "Package 'httr2' is required for mirDIP API access. ",
      "Install with: install.packages('httr2')",
      call. = FALSE
    )
  }

  cache_path <- mirdip_cache_path(gene_symbol, score_class)
  if (file.exists(cache_path)) {
    return(paste(readLines(cache_path, warn = FALSE), collapse = "\n"))
  }

  code <- MIRDIP$score_class_codes[[score_class]]
  req <- httr2::request(mirdip_api_url())
  req <- httr2::req_method(req, "POST")
  req <- httr2::req_body_form(
    req,
    genesymbol = gene_symbol,
    microrna = "",
    scoreClass = code
  )
  req <- httr2::req_timeout(req, 180)

  resp <- tryCatch(
    httr2::req_perform(req),
    error = function(e) {
      stop("mirDIP API request failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  body <- httr2::resp_body_string(resp)
  write_mirdip_cache(cache_path, body)
  body
}

query_targets_mirdip <- function(
  gene_symbol,
  species = "hsa",
  score_class = MIRDIP$default_score_class,
  min_score = THRESHOLDS$mirdip_min_score,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  score_class <- validate_mirdip_score_class(score_class)

  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }
  if (isTRUE(MIRDIP$human_only) && species != "hsa") {
    stop(
      "mirDIP supports human (hsa) only; got species = ", species,
      call. = FALSE
    )
  }

  body <- fetch_mirdip_api(gene_symbol, score_class)
  result <- parse_mirdip_response(body, gene_symbol)

  if (!is.null(min_score)) {
    if (!is.numeric(min_score) || length(min_score) != 1L) {
      stop("min_score must be a single numeric value or NULL", call. = FALSE)
    }
    result <- result[!is.na(result$score_raw) & result$score_raw >= min_score, , drop = FALSE]
    rownames(result) <- NULL
  }

  result
}

query_targets <- query_targets_mirdip
