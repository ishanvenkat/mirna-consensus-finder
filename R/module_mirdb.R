# miRDB prediction module — bulk MirTarget file filter by gene symbol.
# Bulk column 2 is RefSeq accession (NM_/XM_), not HUGO symbol — map via symbol_to_refseq().

source("R/config.R")
source("R/utils_gene_id.R")

mirdb_bulk_path <- function() {
  file.path(get_mirdb_root(), MIRDB$bulk_file)
}

empty_mirdb_result <- function() {
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

resolve_mirdb_mirna_prefix <- function(species) {
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, ". Use 'hsa' or 'mmu'.", call. = FALSE)
  }
  SPECIES[[species]]$mirna_prefix
}

filter_mirdb_bulk <- function(path, refseq_ids, mirna_prefix, min_score) {
  if (!file.exists(path)) {
    stop("miRDB bulk file not found: ", path, call. = FALSE)
  }
  if (length(refseq_ids) == 0L) {
    return(list())
  }

  refseq_set <- unique(refseq_ids)
  open_con <- function(p) {
    if (grepl("\\.gz$", p, ignore.case = TRUE)) {
      gzfile(p, open = "rt")
    } else {
      file(p, open = "rt")
    }
  }

  con <- open_con(path)
  on.exit(close(con), add = TRUE)

  prefix_pattern <- paste0("^", mirna_prefix, "-")
  rows <- list()

  repeat {
    chunk <- readLines(con, n = 50000L, warn = FALSE)
    if (length(chunk) == 0L) {
      break
    }
    hit_idx <- integer()
    for (refseq_id in refseq_set) {
      needle <- paste0("\t", refseq_id, "\t")
      hit_idx <- c(hit_idx, grep(needle, chunk, fixed = TRUE))
    }
    hit_idx <- unique(hit_idx)
    if (length(hit_idx) == 0L) {
      next
    }
    for (line in chunk[hit_idx]) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(parts) < 3L) {
        next
      }
      mirna <- parts[[1L]]
      target <- parts[[2L]]
      if (!target %in% refseq_set) {
        next
      }
      if (!grepl(prefix_pattern, mirna, ignore.case = TRUE)) {
        next
      }
      score <- suppressWarnings(as.numeric(parts[[3L]]))
      if (is.na(score) || score < min_score) {
        next
      }
      rows[[length(rows) + 1L]] <- list(mirna = mirna, score = score)
    }
  }
  rows
}

aggregate_mirdb_mirna_scores <- function(rows) {
  if (length(rows) == 0L) {
    return(rows)
  }
  mirna <- vapply(rows, function(row) row$mirna, character(1L))
  scores <- vapply(rows, function(row) as.numeric(row$score), numeric(1L))
  tapply(scores, mirna, max, na.rm = TRUE)
}

normalize_mirdb_rows <- function(gene_symbol, rows) {
  if (length(rows) == 0L) {
    return(empty_mirdb_result())
  }

  best_scores <- aggregate_mirdb_mirna_scores(rows)
  mirna <- names(best_scores)
  scores <- as.numeric(best_scores)

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = mirna,
    score_raw = scores,
    score_name = MIRDB$score_name,
    score_direction = "higher_is_stronger",
    site_count = NA_integer_,
    source = "mirdb",
    retrieved_at = Sys.time(),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$score_raw), , drop = FALSE]
  ord <- order(out$score_raw, out$mirna, decreasing = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

query_targets_mirdb <- function(
  gene_symbol,
  species = "hsa",
  min_score = THRESHOLDS$mirdb,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  mirna_prefix <- resolve_mirdb_mirna_prefix(species)

  if (!is.numeric(min_score) || length(min_score) != 1L) {
    stop("min_score must be a single numeric value", call. = FALSE)
  }
  if (min_score < MIRDB$display_min_score) {
    stop(
      "min_score must be at least ", MIRDB$display_min_score,
      " (miRDB display minimum)",
      call. = FALSE
    )
  }

  path <- mirdb_bulk_path()
  if (!file.exists(path)) {
    stop(
      "miRDB bulk file not found: ", path, "\n",
      "Run: Rscript scripts/download_mirdb.R",
      call. = FALSE
    )
  }

  refseq_ids <- symbol_to_refseq(gene_symbol, species)
  rows <- filter_mirdb_bulk(path, refseq_ids, mirna_prefix, min_score)
  normalize_mirdb_rows(gene_symbol, rows)
}

query_targets <- query_targets_mirdb
