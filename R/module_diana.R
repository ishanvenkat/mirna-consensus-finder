# DIANA-microT-CDS prediction module — REST API (2023 webserver + legacy) with bulk fallback.

source("R/config.R")
source("R/utils_gene_id.R")

resolve_diana_bulk_key <- function(species, mirna_resource = DIANA$default_mirna_resource) {
  if (!species %in% names(DIANA$bulk_files)) {
    stop("Unsupported species for DIANA bulk: ", species, call. = FALSE)
  }
  files <- DIANA$bulk_files[[species]]
  if (is.null(files[[mirna_resource]])) {
    stop("Unknown DIANA miRNA resource key: ", mirna_resource, call. = FALSE)
  }
  files[[mirna_resource]]
}

diana_bulk_path <- function(species, mirna_resource = DIANA$default_mirna_resource) {
  file.path(get_diana_root(), resolve_diana_bulk_key(species, mirna_resource))
}

diana_api_cache_path <- function(species, gene_symbol, threshold) {
  file.path(
    PATHS$cache,
    "diana",
    paste0(species, "_", gene_symbol, "_", threshold, ".json")
  )
}

empty_diana_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    score_raw = numeric(),
    score_name = character(),
    score_direction = character(),
    site_count = integer(),
    source = character(),
    retrieved_at = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )
}

normalize_diana_rows <- function(gene_symbol, rows) {
  if (length(rows) == 0L) {
    return(empty_diana_result())
  }

  mirna <- vapply(rows, function(row) row$mirna, character(1L))
  scores <- vapply(rows, function(row) as.numeric(row$score), numeric(1L))

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = mirna,
    score_raw = scores,
    score_name = DIANA$score_name,
    score_direction = "higher_is_stronger",
    site_count = NA_integer_,
    source = "diana",
    retrieved_at = Sys.time(),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$score_raw), , drop = FALSE]
  ord <- order(out$score_raw, out$mirna, decreasing = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_diana_legacy_response <- function(body, ensembl_id, threshold) {
  if (grepl("<!DOCTYPE html", body, ignore.case = TRUE)) {
    stop("DIANA legacy API returned HTML (service unavailable)", call. = FALSE)
  }
  body <- trimws(body)
  if (nchar(body) == 0L) {
    return(list())
  }

  lines <- strsplit(body, "\n", fixed = TRUE)[[1L]]
  lines <- lines[nchar(trimws(lines)) > 0L]
  rows <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 3L) {
      next
    }
    score <- suppressWarnings(as.numeric(parts[[length(parts)]]))
    if (is.na(score) || score < threshold) {
      next
    }
    rows[[length(rows) + 1L]] <- list(
      mirna = parts[[1L]],
      score = score
    )
  }
  rows
}

parse_diana_webserver_json <- function(payload, threshold) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to parse DIANA API responses", call. = FALSE)
  }
  parsed <- jsonlite::fromJSON(payload, simplifyVector = FALSE)
  data <- parsed$data
  if (is.null(data) || length(data) == 0L) {
    return(list())
  }

  rows <- list()
  for (item in data) {
    score <- suppressWarnings(as.numeric(item$interaction_score))
    if (is.na(score) || score < threshold) {
      next
    }
    mirna <- item$mirna_name
    if (is.null(mirna) || !nzchar(mirna)) {
      next
    }
    rows[[length(rows) + 1L]] <- list(mirna = mirna, score = score)
  }
  rows
}

diana_api_form_fields <- function(gene_symbol, species, threshold) {
  list(
    mirnas = "",
    genes = gene_symbol,
    mre_region = "",
    mirna_exp_tissue = "",
    gene_exp_tissue = "",
    "filter[microt_version]" = DIANA$microt_version_label,
    "filter[mirna_resource]" = DIANA$mirna_resource_label,
    "filter[mirna_confidence]" = "",
    "filter[species]" = SPECIES[[species]]$diana_species,
    "filter[interaction_score]" = as.character(threshold),
    "filter[additional_support]" = ""
  )
}

fetch_diana_legacy_api <- function(ensembl_id, threshold) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for DIANA REST access", call. = FALSE)
  }

  req <- httr2::request(API$diana_microt_cds)
  req <- httr2::req_url_query(
    req,
    genes = ensembl_id,
    threshold = threshold
  )
  req <- httr2::req_timeout(req, 30)
  resp <- httr2::req_perform(req)
  httr2::resp_body_string(resp)
}

fetch_diana_webserver_api <- function(gene_symbol, species, threshold) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for DIANA REST access", call. = FALSE)
  }

  fields <- diana_api_form_fields(gene_symbol, species, threshold)
  req <- httr2::request(paste0(API$diana_microt_webserver, "/interactions_query"))
  req <- httr2::req_timeout(req, 120)
  req <- httr2::req_body_form(req, !!!fields)
  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_string(resp)

  rows <- parse_diana_webserver_json(body, threshold)
  list(body = body, rows = rows, parsed = jsonlite::fromJSON(body, simplifyVector = FALSE))
}

fetch_diana_webserver_pages <- function(gene_symbol, species, threshold) {
  first <- fetch_diana_webserver_api(gene_symbol, species, threshold)
  rows <- first$rows
  parsed <- first$parsed
  last_page <- parsed$last_page
  if (!is.null(last_page) && is.numeric(last_page) && last_page > 1L) {
    fields <- diana_api_form_fields(gene_symbol, species, threshold)
    for (page in 2:last_page) {
      req <- httr2::request(paste0(API$diana_microt_webserver, "/interactions"))
      req <- httr2::req_timeout(req, 120)
      req <- httr2::req_url_query(req, !!!fields, page = page)
      resp <- httr2::req_perform(req)
      page_rows <- parse_diana_webserver_json(httr2::resp_body_string(resp), threshold)
      rows <- c(rows, page_rows)
    }
  }
  rows
}

read_diana_api_cache <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  payload <- readLines(path, warn = FALSE)
  payload <- paste(payload, collapse = "\n")
  list(body = payload, mtime = file.info(path)$mtime)
}

write_diana_api_cache <- function(path, body) {
  cache_dir <- dirname(path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  writeLines(body, path, useBytes = TRUE)
  invisible(path)
}

query_diana_api <- function(gene_symbol, species, threshold) {
  mock_path <- getOption("mirna.diana.api_response_path")
  if (!is.null(mock_path)) {
    cached <- read_diana_api_cache(mock_path)
    rows <- parse_diana_webserver_json(cached$body, threshold)
    return(normalize_diana_rows(gene_symbol, rows))
  }

  cache_path <- diana_api_cache_path(species, gene_symbol, threshold)
  cached <- read_diana_api_cache(cache_path)
  if (!is.null(cached)) {
    rows <- tryCatch(
      parse_diana_webserver_json(cached$body, threshold),
      error = function(e) NULL
    )
    if (!is.null(rows)) {
      return(normalize_diana_rows(gene_symbol, rows))
    }
  }

  rows <- list()
  cache_body <- NULL

  tryCatch({
    ensembl_id <- symbol_to_ensembl(gene_symbol, species)
    legacy_body <- fetch_diana_legacy_api(ensembl_id, threshold)
    rows <- parse_diana_legacy_response(legacy_body, ensembl_id, threshold)
  }, error = function(e) {
    rows <<- list()
  })

  if (length(rows) == 0L) {
    tryCatch({
      rows <- fetch_diana_webserver_pages(gene_symbol, species, threshold)
      if (length(rows) > 0L) {
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
          stop("Package 'jsonlite' is required to cache DIANA API responses", call. = FALSE)
        }
        cache_body <- jsonlite::toJSON(
          list(data = lapply(rows, function(r) {
            list(mirna_name = r$mirna, interaction_score = as.character(r$score))
          })),
          auto_unbox = TRUE
        )
      }
    }, error = function(e) {
      stop("DIANA REST API failed: ", conditionMessage(e), call. = FALSE)
    })
  }

  if (!is.null(cache_body) && nzchar(cache_body)) {
    write_diana_api_cache(cache_path, cache_body)
  }

  normalize_diana_rows(gene_symbol, rows)
}

filter_diana_bulk_gz <- function(path, ensembl_id, threshold) {
  if (!file.exists(path)) {
    stop("DIANA bulk file not found: ", path, call. = FALSE)
  }

  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1L, warn = FALSE)
  if (length(header) == 0L) {
    return(list())
  }

  cols <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  mirna_col <- match(DIANA$cols$mirna, cols)
  ensembl_col <- match(DIANA$cols$ensembl_gene_id, cols)
  score_col <- match(DIANA$cols$interaction_score, cols)
  if (any(is.na(c(mirna_col, ensembl_col, score_col)))) {
    stop("DIANA bulk file has unexpected columns", call. = FALSE)
  }

  needle <- paste0("\t", ensembl_id, "\t")
  rows <- list()
  repeat {
    chunk <- readLines(con, n = 50000L, warn = FALSE)
    if (length(chunk) == 0L) {
      break
    }
    hit_idx <- grep(needle, chunk, fixed = TRUE)
    if (length(hit_idx) == 0L) {
      next
    }
    for (line in chunk[hit_idx]) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      score <- suppressWarnings(as.numeric(parts[[score_col]]))
      if (is.na(score) || score < threshold) {
        next
      }
      rows[[length(rows) + 1L]] <- list(
        mirna = parts[[mirna_col]],
        score = score
      )
    }
  }
  rows
}

query_diana_bulk <- function(gene_symbol, species, threshold, mirna_resource = DIANA$default_mirna_resource) {
  ensembl_id <- symbol_to_ensembl(gene_symbol, species)
  path <- diana_bulk_path(species, mirna_resource)
  rows <- filter_diana_bulk_gz(path, ensembl_id, threshold)
  normalize_diana_rows(gene_symbol, rows)
}

query_targets_diana <- function(
  gene_symbol,
  species = "hsa",
  threshold = THRESHOLDS$diana_microt,
  mirna_resource = DIANA$default_mirna_resource,
  prefer_api = TRUE,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, ". Use 'hsa' or 'mmu'.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1) {
    stop("threshold must be a single numeric value between 0 and 1", call. = FALSE)
  }

  use_api <- isTRUE(prefer_api) && !isTRUE(getOption("mirna.diana.disable_api"))

  if (use_api) {
    api_result <- tryCatch(
      query_diana_api(gene_symbol, species, threshold),
      error = function(e) e
    )
    if (!inherits(api_result, "error")) {
      return(api_result)
    }
  }

  bulk_result <- tryCatch(
    query_diana_bulk(gene_symbol, species, threshold, mirna_resource = mirna_resource),
    error = function(e) e
  )
  if (!inherits(bulk_result, "error")) {
    return(bulk_result)
  }

  api_msg <- if (exists("api_result") && inherits(api_result, "error")) {
    conditionMessage(api_result)
  } else {
    "API disabled or unavailable"
  }
  bulk_msg <- conditionMessage(bulk_result)
  stop(
    "DIANA query failed.\n",
    "  API: ", api_msg, "\n",
    "  Bulk: ", bulk_msg, "\n",
    "Run: Rscript scripts/download_diana_bulk.R",
    call. = FALSE
  )
}

query_targets <- query_targets_diana
