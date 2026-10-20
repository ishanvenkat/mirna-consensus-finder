# Score normalization and consensus fusion (rank-based; no raw score averaging).

INDEPENDENT_SOURCES <- c("targetscan", "diana", "mirdb")
META_SOURCES <- c("mirdip")

empty_integrated_result <- function() {
  data.frame(
    gene_symbol = character(),
    mirna = character(),
    consensus_score = numeric(),
    consensus_rank = integer(),
    n_sources = integer(),
    integration_method = character(),
    mirdip_score_raw = numeric(),
    mirdip_score_class = character(),
    stringsAsFactors = FALSE
  )
}

source_rank_within <- function(df) {
  if (nrow(df) == 0L) {
    return(df)
  }
  direction <- df$score_direction[[1L]]
  if (identical(direction, "lower_is_stronger")) {
    ord <- order(df$score_raw, df$mirna, decreasing = FALSE)
  } else {
    ord <- order(df$score_raw, df$mirna, decreasing = TRUE)
  }
  df <- df[ord, , drop = FALSE]
  df$source_rank <- seq_len(nrow(df))
  df
}

integrate_scores <- function(
  predictions,
  method = c("rank_fusion", "percentile", "tiered")
) {
  method <- match.arg(method)
  required <- c(
    "gene_symbol", "mirna", "score_raw", "score_direction", "source"
  )
  if (!is.data.frame(predictions)) {
    stop("predictions must be a data.frame", call. = FALSE)
  }
  missing <- setdiff(required, names(predictions))
  if (length(missing) > 0L) {
    stop("predictions missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(predictions) == 0L) {
    return(empty_integrated_result())
  }

  genes <- unique(predictions$gene_symbol)
  if (length(genes) != 1L) {
    stop(
      "integrate_scores expects a single gene_symbol; found ",
      length(genes),
      call. = FALSE
    )
  }
  gene_symbol <- genes[[1L]]

  fusion_sources <- if (method == "tiered") {
    INDEPENDENT_SOURCES
  } else {
    unique(predictions$source)
  }
  fusion_sources <- intersect(fusion_sources, unique(predictions$source))

  ranked <- lapply(fusion_sources, function(src) {
    sub <- predictions[predictions$source == src, , drop = FALSE]
    if (nrow(sub) == 0L) {
      return(NULL)
    }
    sub <- sub[!is.na(sub$score_raw), , drop = FALSE]
    if (nrow(sub) == 0L) {
      return(NULL)
    }
    sub <- sub[!duplicated(sub$mirna), , drop = FALSE]
    source_rank_within(sub)
  })
  ranked <- ranked[!vapply(ranked, is.null, logical(1L))]
  if (length(ranked) == 0L) {
    return(empty_integrated_result())
  }

  mirnas <- unique(unlist(lapply(ranked, function(df) df$mirna)))
  n_sources <- vapply(mirnas, function(m) {
    sum(vapply(ranked, function(df) m %in% df$mirna, logical(1L)))
  }, integer(1L))

  if (method == "percentile") {
    consensus <- vapply(mirnas, function(m) {
      vals <- vapply(ranked, function(df) {
        hit <- df[df$mirna == m, , drop = FALSE]
        if (nrow(hit) == 0L) {
          return(NA_real_)
        }
        n <- nrow(df)
        if (n <= 1L) {
          return(100)
        }
        100 * (1 - (hit$source_rank[[1L]] - 1) / (n - 1))
      }, numeric(1L))
      mean(vals, na.rm = TRUE)
    }, numeric(1L))
  } else {
    consensus <- vapply(mirnas, function(m) {
      sum(vapply(ranked, function(df) {
        hit <- df[df$mirna == m, , drop = FALSE]
        if (nrow(hit) == 0L) {
          return(0)
        }
        as.numeric(nrow(df) - hit$source_rank[[1L]] + 1L)
      }, numeric(1L)))
    }, numeric(1L))
  }

  out <- data.frame(
    gene_symbol = gene_symbol,
    mirna = mirnas,
    consensus_score = consensus,
    n_sources = n_sources,
    integration_method = method,
    stringsAsFactors = FALSE
  )

  mirdip <- predictions[predictions$source == "mirdip", , drop = FALSE]
  if (nrow(mirdip) > 0L) {
    mirdip <- mirdip[!duplicated(mirdip$mirna), , drop = FALSE]
    score_class <- if ("score_class" %in% names(mirdip)) {
      mirdip$score_class
    } else {
      NA_character_
    }
    out$mirdip_score_raw <- vapply(out$mirna, function(m) {
      hit <- mirdip[mirdip$mirna == m, , drop = FALSE]
      if (nrow(hit) == 0L) {
        return(NA_real_)
      }
      hit$score_raw[[1L]]
    }, numeric(1L))
    out$mirdip_score_class <- vapply(out$mirna, function(m) {
      hit <- mirdip[mirdip$mirna == m, , drop = FALSE]
      if (nrow(hit) == 0L || !"score_class" %in% names(hit)) {
        return(NA_character_)
      }
      hit$score_class[[1L]]
    }, character(1L))
  } else {
    out$mirdip_score_raw <- NA_real_
    out$mirdip_score_class <- NA_character_
  }

  ord <- order(out$consensus_score, out$n_sources, out$mirna, decreasing = TRUE)
  out <- out[ord, , drop = FALSE]
  out$consensus_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
