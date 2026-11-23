# Rank-based evaluation of consensus vs single-tool predictions against validation.

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

source("R/config.R")
source("R/utils_gene_id.R")
source("R/cache_predictions.R")
source("R/integrate_scores.R")

EVALUATION_METHODS <- c(
  "consensus_tiered",
  "targetscan",
  "diana",
  "mirdb",
  "mirdip"
)

EVALUATION_METHOD_LABELS <- c(
  consensus_tiered = "Tiered consensus",
  consensus_tiered_multi = "Tiered consensus (n_sources >= 2)",
  targetscan = "TargetScan",
  diana = "DIANA-microT-CDS",
  mirdb = "miRDB",
  mirdip = "mirDIP"
)

evaluation_recall_cache_path <- function(gene_symbol, species, ext = "csv") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  ext <- match.arg(ext, c("csv", "rds"))
  file.path(
    PATHS$cache,
    paste0(gene_symbol, "_", species, "_evaluation_recall_at_k.", ext)
  )
}

evaluation_pr_curve_cache_path <- function(gene_symbol, species, ext = "csv") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  ext <- match.arg(ext, c("csv", "rds"))
  file.path(
    PATHS$cache,
    paste0(gene_symbol, "_", species, "_evaluation_pr_curve.", ext)
  )
}

evaluation_summary_cache_path <- function(gene_symbol, species) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  file.path(
    PATHS$cache,
    paste0(gene_symbol, "_", species, "_evaluation_summary.rds")
  )
}

validation_positives <- function(
  validation_df,
  support_level = c("strong", "any")
) {
  support_level <- match.arg(support_level)
  if (!is.data.frame(validation_df)) {
    stop("validation_df must be a data.frame", call. = FALSE)
  }
  if (nrow(validation_df) == 0L) {
    return(list(
      mirnas = character(),
      n_positives = 0L,
      support_level = support_level,
      n_rows = 0L
    ))
  }
  req <- c("mirna", "experimental_support")
  missing <- setdiff(req, names(validation_df))
  if (length(missing) > 0L) {
    stop(
      "validation_df missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  rows <- validation_df[
    validation_df$experimental_support %in% TRUE,
    ,
    drop = FALSE
  ]
  if (support_level == "strong") {
    if (!"support_level" %in% names(rows)) {
      stop(
        "support_level column required when support_level='strong'",
        call. = FALSE
      )
    }
    rows <- rows[rows$support_level == "strong", , drop = FALSE]
  }
  mirnas <- unique(rows$mirna)
  mirnas <- mirnas[!is.na(mirnas) & nzchar(mirnas)]
  list(
    mirnas = mirnas,
    n_positives = length(mirnas),
    support_level = support_level,
    n_rows = nrow(rows)
  )
}

rank_predictions_by_source <- function(predictions_df, source) {
  if (!is.data.frame(predictions_df)) {
    stop("predictions_df must be a data.frame", call. = FALSE)
  }
  source <- validate_prediction_source(source)
  sub <- predictions_df[predictions_df$source == source, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(character())
  }
  sub <- sub[!is.na(sub$score_raw), , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(character())
  }
  sub <- sub[!duplicated(sub$mirna), , drop = FALSE]
  ranked <- source_rank_within(sub)
  ranked$mirna
}

ranked_list_from_integrated <- function(
  integrated_df,
  subset = c("full", "multi_source")
) {
  subset <- match.arg(subset)
  if (!is.data.frame(integrated_df)) {
    stop("integrated_df must be a data.frame", call. = FALSE)
  }
  if (nrow(integrated_df) == 0L) {
    return(character())
  }
  if (!"consensus_rank" %in% names(integrated_df)) {
    stop("integrated_df must include consensus_rank", call. = FALSE)
  }
  df <- integrated_df
  if (subset == "multi_source") {
    if (!"n_sources" %in% names(df)) {
      stop(
        "integrated_df must include n_sources for multi_source subset",
        call. = FALSE
      )
    }
    df <- df[df$n_sources >= 2L, , drop = FALSE]
    if (nrow(df) == 0L) {
      return(character())
    }
    df <- rank_integrated_consensus_subset(df)
  }
  ord <- order(df$consensus_rank, df$mirna, decreasing = FALSE)
  df <- df[ord, , drop = FALSE]
  df$mirna
}

count_hits_in_prefix <- function(ranked_mirnas, positives, k) {
  k_eff <- min(k, length(ranked_mirnas))
  if (k_eff == 0L || length(positives) == 0L) {
    return(0L)
  }
  prefix <- ranked_mirnas[seq_len(k_eff)]
  sum(prefix %in% positives)
}

recall_at_k_table <- function(
  ranked_mirnas,
  positives,
  k_values,
  method_id
) {
  if (length(k_values) == 0L) {
    return(data.frame(
      method = character(),
      k = integer(),
      k_effective = integer(),
      hits = integer(),
      n_positives = integer(),
      recall = numeric(),
      precision = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  n_pos <- length(positives)
  k_values <- unique(as.integer(k_values))
  k_values <- k_values[k_values > 0L]
  if (length(k_values) == 0L) {
    stop("k_values must contain positive integers", call. = FALSE)
  }
  rows <- lapply(k_values, function(k) {
    k_eff <- min(k, length(ranked_mirnas))
    hits <- count_hits_in_prefix(ranked_mirnas, positives, k)
    recall <- if (n_pos > 0L) hits / n_pos else NA_real_
    precision <- if (k_eff > 0L) hits / k_eff else NA_real_
    data.frame(
      method = method_id,
      k = k,
      k_effective = k_eff,
      hits = hits,
      n_positives = n_pos,
      recall = recall,
      precision = precision,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pr_curve_from_ranking <- function(
  ranked_mirnas,
  positives,
  max_depth = NULL
) {
  n_pos <- length(positives)
  n_ranked <- length(ranked_mirnas)
  if (n_ranked == 0L) {
    return(data.frame(
      depth = integer(),
      hits_cumulative = integer(),
      n_positives = integer(),
      recall = numeric(),
      precision = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  depth_limit <- if (is.null(max_depth)) {
    n_ranked
  } else {
    min(as.integer(max_depth), n_ranked)
  }
  if (depth_limit < 1L) {
    return(data.frame(
      depth = integer(),
      hits_cumulative = integer(),
      n_positives = integer(),
      recall = numeric(),
      precision = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  depths <- seq_len(depth_limit)
  hits_cum <- integer(depth_limit)
  found <- character()
  for (i in depths) {
    m <- ranked_mirnas[[i]]
    if (m %in% positives && !m %in% found) {
      found <- c(found, m)
    }
    hits_cum[[i]] <- length(found)
  }
  data.frame(
    depth = depths,
    hits_cumulative = hits_cum,
    n_positives = n_pos,
    recall = if (n_pos > 0L) hits_cum / n_pos else NA_real_,
    precision = if (depth_limit > 0L) hits_cum / depths else NA_real_,
    stringsAsFactors = FALSE
  )
}

assert_tiered_integrated <- function(integrated_df) {
  if (!"integration_method" %in% names(integrated_df)) {
    stop(
      "integrated cache missing integration_method; re-run export_score_integrated.R",
      call. = FALSE
    )
  }
  methods <- unique(integrated_df$integration_method)
  if (length(methods) != 1L || methods[[1L]] != "tiered") {
    stop(
      "integrated cache integration_method must be 'tiered' for corrected consensus; found: ",
      paste(methods, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

evaluation_method_rankings <- function(
  integrated_df,
  gene_symbol,
  species,
  consensus_subset = c("full", "multi_source", "both"),
  prediction_sources = EVALUATION_METHODS[-1L]
) {
  consensus_subset <- match.arg(consensus_subset)
  assert_tiered_integrated(integrated_df)

  rankings <- list()

  if (consensus_subset %in% c("full", "both")) {
    rankings$consensus_tiered <- ranked_list_from_integrated(
      integrated_df,
      subset = "full"
    )
  }
  if (consensus_subset %in% c("multi_source", "both")) {
    rankings$consensus_tiered_multi <- ranked_list_from_integrated(
      integrated_df,
      subset = "multi_source"
    )
  }

  for (src in prediction_sources) {
    if (!prediction_cache_exists(gene_symbol, species, src)) {
      stop(
        "Prediction cache missing for ",
        src,
        ". Run export_score_",
        src,
        ".R",
        call. = FALSE
      )
    }
    preds <- read_prediction_cache(gene_symbol, species, src)
    rankings[[src]] <- rank_predictions_by_source(preds, src)
  }

  rankings
}

evaluate_consensus_rankings <- function(
  gene_symbol,
  species = "hsa",
  support_level = c("strong", "any"),
  consensus_subset = c("full", "multi_source", "both"),
  k_values = c(10L, 25L, 50L, 100L),
  max_pr_depth = 200L,
  prediction_sources = c("targetscan", "diana", "mirdb", "mirdip")
) {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  species <- validate_cache_species(species)
  support_level <- match.arg(support_level)
  consensus_subset <- match.arg(consensus_subset)

  if (!validation_cache_exists(gene_symbol, species)) {
    stop(
      "Validation cache missing. Run: Rscript scripts/export_validation.R ",
      gene_symbol,
      " ",
      species,
      call. = FALSE
    )
  }
  if (!integrated_cache_exists(gene_symbol, species)) {
    stop(
      "Integrated cache missing. Run: Rscript scripts/export_score_integrated.R ",
      gene_symbol,
      " ",
      species,
      call. = FALSE
    )
  }

  validation <- read_validation_cache(gene_symbol, species)
  integrated <- read_integrated_cache(gene_symbol, species)
  positives <- validation_positives(validation, support_level = support_level)

  rankings <- evaluation_method_rankings(
    integrated,
    gene_symbol,
    species,
    consensus_subset = consensus_subset,
    prediction_sources = prediction_sources
  )

  recall_parts <- lapply(names(rankings), function(mid) {
    recall_at_k_table(
      rankings[[mid]],
      positives$mirnas,
      k_values,
      mid
    )
  })
  recall_at_k <- do.call(rbind, recall_parts)
  rownames(recall_at_k) <- NULL

  pr_parts <- lapply(names(rankings), function(mid) {
    pr <- pr_curve_from_ranking(
      rankings[[mid]],
      positives$mirnas,
      max_depth = max_pr_depth
    )
    if (nrow(pr) > 0L) {
      pr$method <- mid
    }
    pr
  })
  pr_curve <- do.call(rbind, pr_parts)
  rownames(pr_curve) <- NULL

  rank_lengths <- vapply(rankings, length, integer(1L))

  prov_cols <- c(
    "mirtarbase_release",
    "multimir_db_version",
    "tarbase_snapshot"
  )
  prov <- list()
  for (col in prov_cols) {
    if (col %in% names(validation) && nrow(validation) > 0L) {
      vals <- unique(validation[[col]])
      vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
      if (length(vals) > 0L) {
        prov[[col]] <- vals[[1L]]
      }
    }
  }

  metadata <- c(
    list(
      gene_symbol = gene_symbol,
      species = species,
      support_level = support_level,
      consensus_subset = consensus_subset,
      k_values = k_values,
      max_pr_depth = max_pr_depth,
      n_positives = positives$n_positives,
      integration_method = unique(integrated$integration_method),
      evaluated_at = Sys.time()
    ),
    prov
  )

  list(
    positives = positives,
    rankings = rankings,
    recall_at_k = recall_at_k,
    pr_curve = pr_curve,
    rank_lengths = rank_lengths,
    metadata = metadata
  )
}

write_evaluation_cache <- function(result, gene_symbol, species) {
  if (!is.list(result) || !"recall_at_k" %in% names(result)) {
    stop("result must be output of evaluate_consensus_rankings()", call. = FALSE)
  }
  ensure_data_dirs()
  recall_csv <- evaluation_recall_cache_path(gene_symbol, species, "csv")
  pr_csv <- evaluation_pr_curve_cache_path(gene_symbol, species, "csv")
  summary_rds <- evaluation_summary_cache_path(gene_symbol, species)
  write.csv(result$recall_at_k, recall_csv, row.names = FALSE)
  write.csv(result$pr_curve, pr_csv, row.names = FALSE)
  saveRDS(result, summary_rds)
  invisible(list(
    recall_csv = recall_csv,
    pr_curve_csv = pr_csv,
    summary_rds = summary_rds
  ))
}

plot_evaluation_pr_curves <- function(
  pr_curve,
  outfile,
  method_labels = EVALUATION_METHOD_LABELS
) {
  if (!is.data.frame(pr_curve) || !"method" %in% names(pr_curve)) {
    stop("pr_curve must include a method column", call. = FALSE)
  }
  if (nrow(pr_curve) == 0L) {
    stop("pr_curve has no rows to plot", call. = FALSE)
  }
  plot_path <- outfile
  parent <- dirname(normalizePath(plot_path, winslash = "/", mustWork = FALSE))
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }

  methods <- unique(pr_curve$method)
  pal <- c(
    consensus_tiered = "#2563eb",
    consensus_tiered_multi = "#1d4ed8",
    targetscan = "#0ea5e9",
    diana = "#06b6d4",
    mirdb = "#9333ea",
    mirdip = "#ea580c"
  )

  # Open PNG before any graphics calls; par() without a device opens Rplots.pdf in batch mode.
  png(plot_path, width = 800, height = 600, res = 120)
  on.exit({
    dev.off()
  }, add = TRUE)

  xlim <- range(pr_curve$recall, na.rm = TRUE)
  ylim <- range(pr_curve$precision, na.rm = TRUE)
  if (all(!is.finite(xlim))) {
    xlim <- c(0, 1)
  }
  if (all(!is.finite(ylim))) {
    ylim <- c(0, 1)
  }
  xlim[1] <- min(0, xlim[1])
  ylim[1] <- min(0, ylim[1])
  xlim[2] <- max(1, xlim[2])
  ylim[2] <- max(ylim[2], 0.05)

  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "Recall",
    ylab = "Precision",
    main = "Precision–recall by rank depth"
  )
  for (mid in methods) {
    sub <- pr_curve[pr_curve$method == mid, , drop = FALSE]
    ord <- order(sub$recall, sub$depth)
    sub <- sub[ord, , drop = FALSE]
    col <- pal[[mid]] %||% "#64748b"
    lab <- method_labels[[mid]] %||% mid
    lines(sub$recall, sub$precision, col = col, lwd = 2)
  }
  legend(
    "topright",
    legend = vapply(
      methods,
      function(m) method_labels[[m]] %||% m,
      character(1L)
    ),
    col = vapply(methods, function(m) pal[[m]] %||% "#64748b", character(1L)),
    lwd = 2,
    bty = "n",
    cex = 0.85
  )
  invisible(plot_path)
}