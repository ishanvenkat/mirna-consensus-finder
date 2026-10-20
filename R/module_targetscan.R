# TargetScan prediction module — bulk Summary Counts filter by gene symbol.

source("R/config.R")
source("R/utils_gene_id.R")

.targetscan_cache <- new.env(parent = emptyenv())

resolve_targetscan_db <- function(species) {
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, ". Use 'hsa' or 'mmu'.", call. = FALSE)
  }
  SPECIES[[species]]$targetscan_db
}

targetscan_data_path <- function(species, file_key) {
  db_key <- resolve_targetscan_db(species)
  filename <- TARGETSCAN$files[[file_key]]
  if (is.null(filename)) {
    stop("Unknown TargetScan file key: ", file_key, call. = FALSE)
  }

  path <- file.path(get_targetscan_root(), db_key, filename)
  if (!file.exists(path)) {
    stop(
      "TargetScan data file not found: ", path, "\n",
      "Run: Rscript scripts/download_targetscan.R",
      call. = FALSE
    )
  }
  path
}

read_targetscan_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(
      path,
      sep = "\t",
      header = TRUE,
      showProgress = FALSE
    ), stringsAsFactors = FALSE))
  }
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

is_representative_transcript <- function(value) {
  normalized <- toupper(trimws(as.character(value)))
  normalized %in% c("1", "Y", "YES", "TRUE")
}

filter_representative_transcript <- function(summary_df, gene_info_df) {
  gi_cols <- TARGETSCAN$gene_info_cols
  tx_col <- gi_cols$transcript_id
  rep_col <- gi_cols$representative

  if (!all(c(tx_col, rep_col) %in% names(gene_info_df))) {
    stop("Gene info file missing expected columns", call. = FALSE)
  }

  rep_tx <- gene_info_df[is_representative_transcript(gene_info_df[[rep_col]]), , drop = FALSE]
  rep_ids <- unique(rep_tx[[tx_col]])

  summary_df[summary_df[[TARGETSCAN$cols$transcript_id]] %in% rep_ids, , drop = FALSE]
}

resolve_targetscan_score_column <- function(hits, preferred, fallbacks = character()) {
  candidates <- c(preferred, fallbacks)
  for (col in candidates) {
    if (col %in% names(hits)) {
      return(col)
    }
  }
  stop("Summary Counts file missing score column: ", preferred, call. = FALSE)
}

build_targetscan_family_index <- function(family_df, species_tax_id) {
  fcols <- TARGETSCAN$family_info_cols
  required <- c(fcols$family, fcols$seed, fcols$species_id, fcols$mirna)
  missing <- setdiff(required, names(family_df))
  if (length(missing) > 0L) {
    stop("miR_Family_Info missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  sub <- family_df[family_df[[fcols$species_id]] == species_tax_id, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(list(by_seed = list(), by_name = list()))
  }

  mirna_vals <- trimws(as.character(sub[[fcols$mirna]]))
  keep <- nzchar(mirna_vals)
  sub <- sub[keep, , drop = FALSE]
  mirna_vals <- mirna_vals[keep]

  by_seed <- split(mirna_vals, paste(sub[[fcols$seed]], sub[[fcols$species_id]], sep = "\t"))
  by_name <- split(mirna_vals, paste(sub[[fcols$family]], sub[[fcols$species_id]], sep = "\t"))
  list(
    by_seed = lapply(by_seed, unique),
    by_name = lapply(by_name, unique)
  )
}

lookup_targetscan_family_members <- function(family_key, species_tax_id, index) {
  key <- paste(family_key, species_tax_id, sep = "\t")
  if (!is.null(index$by_seed[[key]])) {
    return(index$by_seed[[key]])
  }
  if (!is.null(index$by_name[[key]])) {
    return(index$by_name[[key]])
  }
  NULL
}

expand_targetscan_family_rows <- function(hits, species_tax_id, family_index, cols) {
  if (nrow(hits) == 0L) {
    return(hits)
  }
  fam_col <- cols$mirna_family
  rep_col <- cols$representative_mirna
  mirna_out <- character()
  row_idx <- integer()
  for (i in seq_len(nrow(hits))) {
    family_key <- as.character(hits[[fam_col]][[i]])
    representative <- as.character(hits[[rep_col]][[i]])
    members <- lookup_targetscan_family_members(family_key, species_tax_id, family_index)
    if (is.null(members) || length(members) == 0L) {
      members <- representative
    }
    mirna_out <- c(mirna_out, members)
    row_idx <- c(row_idx, rep(i, length(members)))
  }
  cbind(hits[row_idx, , drop = FALSE], mirna = mirna_out, stringsAsFactors = FALSE)
}

dedupe_targetscan_by_mirna <- function(result, score_direction) {
  if (nrow(result) == 0L) {
    return(result)
  }
  mirnas <- unique(result$mirna)
  rows <- lapply(mirnas, function(m) {
    sub <- result[result$mirna == m, , drop = FALSE]
    if (nrow(sub) == 1L) {
      return(sub)
    }
    if (identical(score_direction, "lower_is_stronger")) {
      best <- which.min(sub$score_raw)
    } else {
      best <- which.max(sub$score_raw)
    }
    sub[best, , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

load_targetscan_tables <- function(species) {
  cache_key <- resolve_targetscan_db(species)
  if (!is.null(.targetscan_cache[[cache_key]])) {
    return(.targetscan_cache[[cache_key]])
  }

  gene_info <- read_targetscan_table(targetscan_data_path(species, "gene_info"))
  summary_all <- read_targetscan_table(targetscan_data_path(species, "summary_all"))
  summary_default <- read_targetscan_table(targetscan_data_path(species, "summary_default"))
  family_info <- read_targetscan_table(targetscan_data_path(species, "family_info"))
  species_tax_id <- SPECIES[[species]]$ncbi_tax_id
  family_index <- build_targetscan_family_index(family_info, species_tax_id)

  tables <- list(
    gene_info = gene_info,
    summary_all = summary_all,
    summary_default = summary_default,
    family_info = family_info,
    family_index = family_index,
    species_tax_id = species_tax_id
  )
  .targetscan_cache[[cache_key]] <- tables
  tables
}

validate_mirna_prefix <- function(prefix) {
  if (is.null(prefix)) {
    return(NULL)
  }
  if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix)) {
    stop("mirna_prefix must be a single non-empty string or NULL", call. = FALSE)
  }
  if (!grepl("^[A-Za-z0-9]+$", prefix)) {
    stop("mirna_prefix must contain only letters and digits", call. = FALSE)
  }
  prefix
}

validate_mirna_query <- function(mirna) {
  if (is.null(mirna)) {
    return(NULL)
  }
  if (!is.character(mirna) || length(mirna) != 1L || !nzchar(trimws(mirna))) {
    stop("mirna must be a single non-empty string or NULL", call. = FALSE)
  }
  trimws(mirna)
}

validate_mirna_family_query <- function(mirna_family) {
  if (is.null(mirna_family)) {
    return(NULL)
  }
  if (!is.character(mirna_family) || length(mirna_family) != 1L || !nzchar(trimws(mirna_family))) {
    stop("mirna_family must be a single non-empty string or NULL", call. = FALSE)
  }
  trimws(mirna_family)
}

resolve_mirna_prefix <- function(mirna_prefix, species) {
  prefix <- validate_mirna_prefix(mirna_prefix)
  if (is.null(prefix)) {
    return(NULL)
  }
  if (identical(tolower(prefix), "species")) {
    if (!species %in% names(SPECIES)) {
      stop("Unsupported species: ", species, call. = FALSE)
    }
    return(SPECIES[[species]]$mirna_prefix)
  }
  prefix
}

filter_targetscan_summary_mirna <- function(hits, cols, mirna, mirna_family) {
  if (!is.null(mirna)) {
    col <- hits[[cols$representative_mirna]]
    hits <- hits[col == mirna, , drop = FALSE]
  }
  if (!is.null(mirna_family)) {
    fam_col <- cols$mirna_family
    if (!fam_col %in% names(hits)) {
      stop("Summary Counts file missing miRNA family column: ", fam_col, call. = FALSE)
    }
    hits <- hits[hits[[fam_col]] == mirna_family, , drop = FALSE]
  }
  hits
}

filter_targetscan_mirna_members <- function(result, mirna_prefix, mirna) {
  if (!is.null(mirna_prefix)) {
    pattern <- paste0("^", mirna_prefix, "-")
    result <- result[grepl(pattern, result$mirna, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(mirna)) {
    result <- result[result$mirna == mirna, , drop = FALSE]
  }
  result
}

targetscan_score_direction <- function(score_name) {
  if (score_name == "cwcs") {
    return("lower_is_stronger")
  }
  "higher_is_stronger"
}

apply_targetscan_score_threshold <- function(result, min_cwcs, score_name) {
  if (is.null(min_cwcs)) {
    return(result)
  }
  if (score_name == "cwcs") {
    # CWCS: more negative = stronger; min_cwcs is the weakest (highest) value kept.
    return(result[!is.na(result$score_raw) & result$score_raw <= min_cwcs, , drop = FALSE])
  }
  result[!is.na(result$score_raw) & result$score_raw >= min_cwcs, , drop = FALSE]
}

sort_targetscan_result <- function(result) {
  if (nrow(result) == 0L) {
    return(result)
  }
  direction <- result$score_direction[[1L]]
  if (identical(direction, "lower_is_stronger")) {
    ord <- order(result$score_raw, result$mirna, decreasing = FALSE)
  } else {
    ord <- order(result$score_raw, result$mirna, decreasing = TRUE)
  }
  result <- result[ord, , drop = FALSE]
  rownames(result) <- NULL
  result
}

empty_targetscan_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    score_raw = numeric(),
    score_name = character(),
    score_direction = character(),
    site_count = integer(),
    source = character(),
    retrieved_at = as.POSIXct(character()),
    aggregate_pct = numeric(),
    mirna_family = character(),
    targetscan_representative_mirna = character(),
    stringsAsFactors = FALSE
  )
}

query_targets_targetscan <- function(
  gene_symbol,
  species = "hsa",
  predictions = c("all", "default"),
  score_name = c("cwcs", "biochem_occupancy"),
  min_cwcs = NULL,
  mirna_prefix = NULL,
  mirna = NULL,
  mirna_family = NULL,
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  predictions <- match.arg(predictions)
  score_name <- match.arg(score_name)
  mirna_prefix <- resolve_mirna_prefix(mirna_prefix, species)
  mirna <- validate_mirna_query(mirna)
  mirna_family <- validate_mirna_family_query(mirna_family)

  if (!is.null(min_cwcs) && (!is.numeric(min_cwcs) || length(min_cwcs) != 1L)) {
    stop("min_cwcs must be a single numeric value or NULL", call. = FALSE)
  }
  if (!is.null(min_cwcs) && score_name == "biochem_occupancy") {
    stop(
      "min_cwcs applies to CWCS only; use score_name = \"cwcs\" or filter biochem scores after query",
      call. = FALSE
    )
  }

  tables <- load_targetscan_tables(species)
  summary_key <- if (predictions == "all") "summary_all" else "summary_default"
  summary_df <- filter_representative_transcript(
    tables[[summary_key]],
    tables$gene_info
  )

  cols <- TARGETSCAN$cols
  species_tax_id <- tables$species_tax_id
  species_col <- cols$species_id
  if (!species_col %in% names(summary_df)) {
    stop("Summary Counts file missing species column: ", species_col, call. = FALSE)
  }

  hits <- summary_df[summary_df[[cols$gene_symbol]] == gene_symbol, , drop = FALSE]
  hits <- hits[hits[[species_col]] == species_tax_id, , drop = FALSE]
  hits <- filter_targetscan_summary_mirna(hits, cols, mirna, mirna_family)
  if (nrow(hits) == 0L) {
    return(empty_targetscan_result())
  }

  if (score_name == "cwcs") {
    score_col <- resolve_targetscan_score_column(hits, cols$cwcs)
  } else {
    score_col <- resolve_targetscan_score_column(
      hits,
      cols$biochem_occupancy,
      fallbacks = c(
        "Predicted occupancy - transfected miRNA",
        "Predicted occupancy (transfected miRNA)"
      )
    )
  }

  score_direction <- targetscan_score_direction(score_name)
  expanded <- expand_targetscan_family_rows(
    hits,
    species_tax_id,
    tables$family_index,
    cols
  )
  site_count <- expanded[[cols$conserved_sites]] + expanded[[cols$nonconserved_sites]]
  result <- data.frame(
    gene_symbol = gene_symbol,
    mirna = as.character(expanded$mirna),
    score_raw = as.numeric(expanded[[score_col]]),
    score_name = score_name,
    score_direction = score_direction,
    site_count = as.integer(site_count),
    source = "targetscan",
    retrieved_at = Sys.time(),
    aggregate_pct = as.numeric(expanded[[cols$aggregate_pct]]),
    mirna_family = as.character(expanded[[cols$mirna_family]]),
    targetscan_representative_mirna = as.character(expanded[[cols$representative_mirna]]),
    stringsAsFactors = FALSE
  )

  result <- dedupe_targetscan_by_mirna(result, score_direction)
  result <- filter_targetscan_mirna_members(result, mirna_prefix, mirna)
  result <- apply_targetscan_score_threshold(result, min_cwcs, score_name)
  sort_targetscan_result(result)
}

query_targets <- query_targets_targetscan
