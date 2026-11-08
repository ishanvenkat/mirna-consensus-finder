# Experimental validation overlay — miRTarBase and DIANA-TarBase.
# Read-only reference layer; never mixed into prediction scores.

source("R/config.R")
source("R/utils_gene_id.R")

VALIDATION_SOURCES <- c("mirtarbase", "tarbase")

empty_validation_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    experimental_support = logical(),
    support_level = character(),
    validation_methods = character(),
    validation_sources = character(),
    binding_site_reported = logical(),
    mirtarbase_release = character(),
    multimir_db_version = character(),
    tarbase_snapshot = character(),
    retrieved_at = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

is_strong_mirtarbase_method <- function(text) {
  if (length(text) == 0L || is.na(text)) {
    return(FALSE)
  }
  text <- as.character(text)
  if (!nzchar(trimws(text))) {
    return(FALSE)
  }
  grepl(MIRTARBASE$strong_method_pattern, text, ignore.case = TRUE)
}

read_mirtarbase_table <- function(path) {
  if (!file.exists(path)) {
    stop("miRTarBase bulk file not found: ", path, call. = FALSE)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(
      path,
      data.table = FALSE,
      showProgress = FALSE
    )))
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mirtarbase_species_code <- function(species) {
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }
  SPECIES[[species]]$mirna_prefix
}

mirtarbase_row_strong <- function(row) {
  exp_col <- MIRTARBASE$experiments_col
  rep_col <- MIRTARBASE$reporter_flag_col
  wb_col <- MIRTARBASE$western_flag_col

  from_text <- if (exp_col %in% names(row)) {
    is_strong_mirtarbase_method(row[[exp_col]])
  } else {
    FALSE
  }
  from_flags <- FALSE
  if (rep_col %in% names(row) && !is.na(row[[rep_col]])) {
    from_flags <- from_flags || identical(as.character(row[[rep_col]]), "1")
  }
  if (wb_col %in% names(row) && !is.na(row[[wb_col]])) {
    from_flags <- from_flags || identical(as.character(row[[wb_col]]), "1")
  }
  from_text || from_flags
}

resolve_mirtarbase_bulk_path <- function(mirtarbase_use_local) {
  if (!isTRUE(mirtarbase_use_local)) {
    return(NULL)
  }
  override <- getOption("mirna.mirtarbase.bulk_path")
  if (!is.null(override) && nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  full_path <- mirtarbase_bulk_path("mti")
  if (file.exists(full_path) && file.info(full_path)$size > 0L) {
    return(full_path)
  }
  se_path <- mirtarbase_bulk_path("se_reporter")
  if (file.exists(se_path) && file.info(se_path)$size > 0L) {
    return(se_path)
  }
  NULL
}

query_validation_mirtarbase_local <- function(
  gene_symbol,
  species = "hsa",
  mirtarbase_strong_methods = c("reporter_western", "all")
) {
  mirtarbase_strong_methods <- match.arg(mirtarbase_strong_methods)
  bulk_path <- resolve_mirtarbase_bulk_path(TRUE)
  if (is.null(bulk_path)) {
    stop(
      "miRTarBase local bulk not found. Run: ",
      "Rscript scripts/download_mirtarbase.R",
      call. = FALSE
    )
  }

  tbl <- read_mirtarbase_table(bulk_path)
  target_col <- MIRTARBASE$target_col
  species_col <- MIRTARBASE$species_col
  mirna_col <- MIRTARBASE$mirna_col
  exp_col <- MIRTARBASE$experiments_col

  if (!target_col %in% names(tbl)) {
    stop("Unexpected miRTarBase columns in ", basename(bulk_path), call. = FALSE)
  }

  org <- mirtarbase_species_code(species)
  hits <- tbl[toupper(trimws(tbl[[target_col]])) == gene_symbol, , drop = FALSE]
  if (species_col %in% names(hits)) {
    hits <- hits[
      tolower(trimws(as.character(hits[[species_col]]))) == org,
      ,
      drop = FALSE
    ]
  }

  rows <- list()
  if (nrow(hits) == 0L) {
    return(rows)
  }

  for (i in seq_len(nrow(hits))) {
    row <- hits[i, , drop = FALSE]
    if (mirtarbase_strong_methods == "reporter_western" && !mirtarbase_row_strong(row)) {
      next
    }
    methods <- if (exp_col %in% names(row)) {
      as.character(row[[exp_col]][[1L]])
    } else {
      ""
    }
    rows[[length(rows) + 1L]] <- list(
      mirna = as.character(row[[mirna_col]][[1L]]),
      support_level = "strong",
      methods = methods,
      sources = "miRTarBase",
      binding_site_reported = FALSE
    )
  }
  rows
}

multimir_db_metadata <- function() {
  if (!requireNamespace("multiMiR", quietly = TRUE)) {
    return(list(multimir_db_version = NA_character_, tarbase_snapshot = NA_character_))
  }
  multimir_db_version <- NA_character_
  tarbase_snapshot <- NA_character_
  info <- tryCatch(multiMiR::multimir_dbInfo(), error = function(e) NULL)
  if (!is.null(info) && is.data.frame(info) && nrow(info) > 0L) {
    name_col <- if ("source_name" %in% names(info)) {
      "source_name"
    } else if ("database" %in% names(info)) {
      "database"
    } else {
      NULL
    }
    ver_col <- if ("source_version" %in% names(info)) {
      "source_version"
    } else if ("version" %in% names(info)) {
      "version"
    } else {
      NULL
    }
    map_col <- if ("map_name" %in% names(info)) "map_name" else name_col
    if (!is.null(map_col) && !is.null(ver_col)) {
      tb <- info[tolower(info[[map_col]]) == "tarbase", , drop = FALSE]
      if (nrow(tb) > 0L) {
        label <- if (!is.null(name_col)) tb[[name_col]][[1L]] else "TarBase"
        tarbase_snapshot <- paste0(label, " ", tb[[ver_col]][[1L]])
      }
    }
  }
  ver <- tryCatch(multiMiR::multimir_dbInfoVersions(), error = function(e) NULL)
  if (!is.null(ver) && is.data.frame(ver) && nrow(ver) > 0L) {
    if ("VERSION" %in% names(ver)) {
      multimir_db_version <- as.character(ver$VERSION[[1L]])
    } else if ("version" %in% names(ver)) {
      multimir_db_version <- as.character(ver$version[[1L]])
    } else if ("db_version" %in% names(ver)) {
      multimir_db_version <- as.character(ver$db_version[[1L]])
    }
  }
  list(multimir_db_version = multimir_db_version, tarbase_snapshot = tarbase_snapshot)
}

query_validation_tarbase_multimir <- function(gene_symbol, species) {
  if (!requireNamespace("multiMiR", quietly = TRUE)) {
    stop(
      "Package 'multiMiR' is required for TarBase lookup. ",
      "Install with BiocManager::install('multiMiR').",
      call. = FALSE
    )
  }

  mm <- multiMiR::get_multimir(
    org = species,
    target = gene_symbol,
    table = "tarbase",
    summary = FALSE
  )

  rows <- list()
  if (is.null(mm@data$tarbase) || nrow(mm@data$tarbase) == 0L) {
    return(rows)
  }

  tb <- mm@data$tarbase
  mirna_col <- if ("mature_mirna_id" %in% names(tb)) {
    "mature_mirna_id"
  } else if ("mirna" %in% names(tb)) {
    "mirna"
  } else {
    "mature_mirna_acc"
  }
  target_col <- if ("target_symbol" %in% names(tb)) {
    "target_symbol"
  } else {
    "gene"
  }
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
      methods = if (method_col %in% names(tb)) {
        as.character(tb[[method_col]][[i]])
      } else {
        ""
      },
      sources = "TarBase",
      binding_site_reported = binding
    )
  }
  rows
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

normalize_validation_rows <- function(
  gene_symbol,
  rows,
  provenance = list()
) {
  if (length(rows) == 0L) {
    out <- empty_validation_result()
    attr(out, "provenance") <- provenance
    return(out)
  }

  mirna <- vapply(rows, function(r) r$mirna, character(1L))
  support_level <- vapply(rows, function(r) r$support_level %||% NA_character_, character(1L))
  methods <- vapply(rows, function(r) r$methods %||% "", character(1L))
  sources <- vapply(rows, function(r) r$sources, character(1L))
  binding <- vapply(rows, function(r) isTRUE(r$binding_site_reported), logical(1L))

  mirtarbase_release <- provenance$mirtarbase_release %||% NA_character_
  multimir_db_version <- provenance$multimir_db_version %||% NA_character_
  tarbase_snapshot <- provenance$tarbase_snapshot %||% NA_character_

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = mirna,
    experimental_support = TRUE,
    support_level = support_level,
    validation_methods = methods,
    validation_sources = sources,
    binding_site_reported = binding,
    mirtarbase_release = mirtarbase_release,
    multimir_db_version = multimir_db_version,
    tarbase_snapshot = tarbase_snapshot,
    retrieved_at = Sys.time(),
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out[, c("gene_symbol", "mirna")]), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "provenance") <- provenance
  out
}

query_validation <- function(
  gene_symbol,
  species = "hsa",
  sources = VALIDATION_SOURCES,
  mirtarbase_use_local = TRUE,
  mirtarbase_strong_methods = "reporter_western",
  ...
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }
  sources <- intersect(sources, VALIDATION_SOURCES)
  if (length(sources) == 0L) {
    stop("No validation sources selected", call. = FALSE)
  }

  fixture_path <- getOption("mirna.validation.fixture_path")
  if (!is.null(fixture_path) && nzchar(fixture_path)) {
    rows <- parse_validation_fixture(fixture_path, gene_symbol)
    rows <- merge_validation_rows(rows)
    prov <- list(
      mirtarbase_release = "fixture",
      multimir_db_version = NA_character_,
      tarbase_snapshot = NA_character_
    )
    return(normalize_validation_rows(gene_symbol, rows, provenance = prov))
  }

  if (isTRUE(getOption("mirna.validation.disable_multimir"))) {
    if ("tarbase" %in% sources) {
      stop(
        "Validation lookup disabled (mirna.validation.disable_multimir = TRUE)",
        call. = FALSE
      )
    }
  }

  rows <- list()
  provenance <- list(
    mirtarbase_release = paste0(
      MIRTARBASE$release_label,
      " / ",
      MIRTARBASE$db_version
    ),
    multimir_db_version = NA_character_,
    tarbase_snapshot = NA_character_
  )

  if ("mirtarbase" %in% sources) {
    mirtarbase_rows <- query_validation_mirtarbase_local(
      gene_symbol,
      species = species,
      mirtarbase_strong_methods = mirtarbase_strong_methods
    )
    rows <- c(rows, mirtarbase_rows)
  }

  if ("tarbase" %in% sources) {
    if (isTRUE(getOption("mirna.validation.disable_multimir"))) {
      stop(
        "TarBase requires multiMiR; mirna.validation.disable_multimir is TRUE",
        call. = FALSE
      )
    }
    tarbase_rows <- query_validation_tarbase_multimir(gene_symbol, species)
    rows <- c(rows, tarbase_rows)
    meta <- multimir_db_metadata()
    provenance$multimir_db_version <- meta$multimir_db_version
    provenance$tarbase_snapshot <- meta$tarbase_snapshot
  }

  rows <- merge_validation_rows(rows)
  normalize_validation_rows(gene_symbol, rows, provenance = provenance)
}

VALIDATION_OVERLAY_COLS <- c(
  "experimental_support",
  "support_level",
  "validation_methods",
  "validation_sources",
  "binding_site_reported"
)

init_validation_overlay_columns <- function(df) {
  n <- nrow(df)
  if (n == 0L) {
    df$experimental_support <- logical()
    df$support_level <- character()
    df$validation_methods <- character()
    df$validation_sources <- character()
    df$binding_site_reported <- logical()
    df$row_type <- character()
    return(df)
  }
  df$experimental_support <- FALSE
  df$support_level <- NA_character_
  df$validation_methods <- NA_character_
  df$validation_sources <- NA_character_
  df$binding_site_reported <- FALSE
  df$row_type <- "predicted"
  df
}

validation_lookup_by_mirna <- function(validation) {
  if (nrow(validation) == 0L) {
    return(validation)
  }
  cols <- c("mirna", VALIDATION_OVERLAY_COLS)
  cols <- intersect(cols, names(validation))
  val <- validation[, cols, drop = FALSE]
  val[!duplicated(val$mirna), , drop = FALSE]
}

annotate_integrated_with_validation <- function(
  integrated,
  validation,
  append_validated_only = FALSE
) {
  if (!is.data.frame(integrated)) {
    stop("integrated must be a data.frame", call. = FALSE)
  }
  if (!is.data.frame(validation)) {
    stop("validation must be a data.frame", call. = FALSE)
  }
  if (!"gene_symbol" %in% names(integrated) || !"mirna" %in% names(integrated)) {
    stop("integrated must include gene_symbol and mirna", call. = FALSE)
  }

  provenance <- attr(validation, "provenance")

  if (nrow(integrated) == 0L) {
    out <- init_validation_overlay_columns(integrated)
    attr(out, "validation_provenance") <- provenance
    return(out)
  }

  genes_int <- unique(integrated$gene_symbol)
  if (length(genes_int) != 1L) {
    stop("integrated must have a single gene_symbol", call. = FALSE)
  }
  gene_int <- toupper(trimws(genes_int[[1L]]))

  if (nrow(validation) > 0L) {
    if (!"gene_symbol" %in% names(validation)) {
      stop("validation must include gene_symbol", call. = FALSE)
    }
    genes_val <- unique(validation$gene_symbol)
    if (length(genes_val) != 1L) {
      stop("validation must have a single gene_symbol", call. = FALSE)
    }
    if (toupper(trimws(genes_val[[1L]])) != gene_int) {
      stop("validation gene_symbol must match integrated", call. = FALSE)
    }
  }

  out <- integrated
  rank_before <- out$consensus_rank
  score_before <- if ("consensus_score" %in% names(out)) out$consensus_score else NULL

  out <- init_validation_overlay_columns(out)

  val_lookup <- validation_lookup_by_mirna(validation)
  if (nrow(validation) == 0L) {
    attr(out, "validation_provenance") <- provenance
    return(out)
  }

  val_idx <- match(out$mirna, val_lookup$mirna)
  matched <- !is.na(val_idx)

  out$experimental_support <- matched
  out$support_level <- ifelse(matched, val_lookup$support_level[val_idx], NA_character_)
  out$validation_methods <- ifelse(
    matched,
    val_lookup$validation_methods[val_idx],
    NA_character_
  )
  out$validation_sources <- ifelse(
    matched,
    val_lookup$validation_sources[val_idx],
    NA_character_
  )
  out$binding_site_reported <- ifelse(
    matched,
    val_lookup$binding_site_reported[val_idx],
    FALSE
  )

  if (!identical(out$consensus_rank, rank_before)) {
    stop("consensus_rank must not change during validation overlay", call. = FALSE)
  }
  if ("consensus_score" %in% names(out) &&
      !identical(out$consensus_score, score_before)) {
    stop("consensus_score must not change during validation overlay", call. = FALSE)
  }

  if (isTRUE(append_validated_only)) {
    extra <- setdiff(val_lookup$mirna, out$mirna)
    if (length(extra) > 0L) {
      append_rows <- lapply(extra, function(m) {
        hit <- val_lookup[val_lookup$mirna == m, , drop = FALSE][1L, , drop = FALSE]
        row <- as.list(out[1L, , drop = FALSE])
        row$gene_symbol <- integrated$gene_symbol[[1L]]
        row$mirna <- m
        row$experimental_support <- TRUE
        row$support_level <- hit$support_level[[1L]]
        row$validation_methods <- hit$validation_methods[[1L]]
        row$validation_sources <- hit$validation_sources[[1L]]
        row$binding_site_reported <- hit$binding_site_reported[[1L]]
        row$row_type <- "validated_only"
        consensus_cols <- c(
          "consensus_score", "consensus_rank", "n_sources", "integration_method",
          "mirdip_score_raw", "mirdip_score_class", "sources_hit", "consensus_tier"
        )
        for (col in consensus_cols) {
          if (!col %in% names(row)) {
            next
          }
          val <- row[[col]]
          if (is.integer(val)) {
            row[[col]] <- NA_integer_
          } else if (is.numeric(val)) {
            row[[col]] <- NA_real_
          } else if (is.logical(val)) {
            row[[col]] <- NA
          } else {
            row[[col]] <- NA_character_
          }
        }
        as.data.frame(row, stringsAsFactors = FALSE)
      })
      append_df <- do.call(rbind, append_rows)
      rownames(append_df) <- NULL
      out <- rbind(out, append_df)
    }
  }

  attr(out, "validation_provenance") <- provenance
  out
}
