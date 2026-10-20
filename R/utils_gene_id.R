# Gene symbol validation and Ensembl ID mapping.
# DIANA-microT-CDS legacy REST API requires Ensembl gene IDs; bulk files use ENSG as well.

HUGO_SYMBOL_PATTERN <- "^[A-Z0-9][A-Z0-9-]{0,38}$"
ENSEMBL_GENE_PATTERN <- "^ENSG[0-9]+$"
REFSEQ_ACCESSION_PATTERN <- "^N[MRXY][_][0-9A-Z]+$"

validate_gene_symbol <- function(gene_symbol) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L) {
    stop("gene_symbol must be a single character string", call. = FALSE)
  }
  symbol <- toupper(trimws(gene_symbol))
  if (!grepl(HUGO_SYMBOL_PATTERN, symbol)) {
    stop("Invalid gene symbol format: ", gene_symbol, call. = FALSE)
  }
  symbol
}

validate_ensembl_gene_id <- function(ensembl_id) {
  if (!is.character(ensembl_id) || length(ensembl_id) != 1L) {
    stop("ensembl_id must be a single character string", call. = FALSE)
  }
  id <- trimws(ensembl_id)
  if (!grepl(ENSEMBL_GENE_PATTERN, id)) {
    stop("Invalid Ensembl gene ID format: ", ensembl_id, call. = FALSE)
  }
  id
}

read_ensembl_cache <- function() {
  path <- get_ensembl_cache_path()
  if (!file.exists(path)) {
    return(data.frame(
      species = character(),
      gene_symbol = character(),
      ensembl_gene_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("species", "gene_symbol", "ensembl_gene_id")
  if (!all(needed %in% names(df))) {
    stop("Ensembl cache missing columns: ", paste(setdiff(needed, names(df)), collapse = ", "), call. = FALSE)
  }
  df$gene_symbol <- toupper(trimws(df$gene_symbol))
  df
}

append_ensembl_cache <- function(species, gene_symbol, ensembl_gene_id) {
  path <- get_ensembl_cache_path()
  cache_dir <- dirname(path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  row <- data.frame(
    species = species,
    gene_symbol = gene_symbol,
    ensembl_gene_id = ensembl_gene_id,
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    existing <- read_ensembl_cache()
    existing <- existing[!(existing$species == species & existing$gene_symbol == gene_symbol), , drop = FALSE]
    row <- rbind(existing, row)
  }
  write.table(row, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(path)
}

lookup_ensembl_biomart <- function(gene_symbol, species) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package 'biomaRt' is required to map ", gene_symbol, " to Ensembl. ",
      "Install with BiocManager::install('biomaRt') or use DIANA bulk fallback.",
      call. = FALSE
    )
  }

  dataset <- switch(
    species,
    hsa = "hsapiens_gene_ensembl",
    mmu = "mmusculus_gene_ensembl",
    stop("Unsupported species for Ensembl mapping: ", species, call. = FALSE)
  )
  symbol_attr <- if (species == "hsa") "hgnc_symbol" else "mgi_symbol"

  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  hits <- biomaRt::getBM(
    attributes = c(symbol_attr, "ensembl_gene_id", "gene_biotype"),
    filters = symbol_attr,
    values = gene_symbol,
    mart = mart
  )
  if (nrow(hits) == 0L) {
    stop("No Ensembl gene ID found for symbol ", gene_symbol, " (", species, ")", call. = FALSE)
  }

  exact <- hits[toupper(trimws(hits[[symbol_attr]])) == gene_symbol, , drop = FALSE]
  if (nrow(exact) == 0L) {
    stop("No Ensembl gene ID found for symbol ", gene_symbol, " (", species, ")", call. = FALSE)
  }

  coding <- exact[exact$gene_biotype == "protein_coding", , drop = FALSE]
  pick <- if (nrow(coding) > 0L) coding else exact
  validate_ensembl_gene_id(pick$ensembl_gene_id[[1L]])
}

validate_refseq_accession <- function(refseq_id) {
  if (!is.character(refseq_id) || length(refseq_id) != 1L) {
    stop("refseq_id must be a single character string", call. = FALSE)
  }
  id <- trimws(refseq_id)
  if (!grepl(REFSEQ_ACCESSION_PATTERN, id)) {
    stop("Invalid RefSeq accession format: ", refseq_id, call. = FALSE)
  }
  id
}

read_refseq_cache <- function() {
  path <- get_refseq_cache_path()
  if (!file.exists(path)) {
    return(data.frame(
      species = character(),
      gene_symbol = character(),
      refseq_accession = character(),
      stringsAsFactors = FALSE
    ))
  }
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("species", "gene_symbol", "refseq_accession")
  if (!all(needed %in% names(df))) {
    stop(
      "RefSeq cache missing columns: ",
      paste(setdiff(needed, names(df)), collapse = ", "),
      call. = FALSE
    )
  }
  df$gene_symbol <- toupper(trimws(df$gene_symbol))
  df$refseq_accession <- trimws(df$refseq_accession)
  df
}

append_refseq_cache <- function(species, gene_symbol, refseq_ids) {
  path <- get_refseq_cache_path()
  cache_dir <- dirname(path)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  refseq_ids <- unique(vapply(refseq_ids, validate_refseq_accession, character(1L)))
  rows <- data.frame(
    species = species,
    gene_symbol = gene_symbol,
    refseq_accession = refseq_ids,
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    existing <- read_refseq_cache()
    existing <- existing[!(existing$species == species & existing$gene_symbol == gene_symbol), , drop = FALSE]
    rows <- rbind(existing, rows)
  }
  write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(path)
}

lookup_refseq_orgdb <- function(gene_symbol, species) {
  org_pkg <- switch(
    species,
    hsa = "org.Hs.eg.db",
    mmu = "org.Mm.eg.db",
    stop("Unsupported species for RefSeq mapping: ", species, call. = FALSE)
  )
  if (!requireNamespace(org_pkg, quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(NULL)
  }
  org_db <- getExportedValue(org_pkg, org_pkg)
  keytype <- if (species == "hsa") "SYMBOL" else "SYMBOL"
  hits <- AnnotationDbi::select(
    org_db,
    keys = gene_symbol,
    columns = "REFSEQ",
    keytype = keytype
  )
  if (is.null(hits) || nrow(hits) == 0L) {
    return(character())
  }
  refseq <- unique(trimws(hits$REFSEQ))
  refseq <- refseq[!is.na(refseq) & nzchar(refseq)]
  refseq[grepl(REFSEQ_ACCESSION_PATTERN, refseq)]
}

lookup_refseq_biomart <- function(gene_symbol, species) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    return(character())
  }

  dataset <- switch(
    species,
    hsa = "hsapiens_gene_ensembl",
    mmu = "mmusculus_gene_ensembl",
    stop("Unsupported species for RefSeq mapping: ", species, call. = FALSE)
  )
  symbol_attr <- if (species == "hsa") "hgnc_symbol" else "mgi_symbol"

  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  hits <- biomaRt::getBM(
    attributes = c(symbol_attr, "refseq_mrna"),
    filters = symbol_attr,
    values = gene_symbol,
    mart = mart
  )
  if (nrow(hits) == 0L) {
    return(character())
  }

  exact <- hits[toupper(trimws(hits[[symbol_attr]])) == gene_symbol, , drop = FALSE]
  if (nrow(exact) == 0L) {
    return(character())
  }
  refseq <- unique(trimws(exact$refseq_mrna))
  refseq <- refseq[!is.na(refseq) & nzchar(refseq)]
  refseq[grepl(REFSEQ_ACCESSION_PATTERN, refseq)]
}

symbol_to_refseq <- function(gene_symbol, species = "hsa") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }

  cache <- read_refseq_cache()
  hit <- cache[cache$species == species & cache$gene_symbol == gene_symbol, , drop = FALSE]
  if (nrow(hit) > 0L) {
    return(unique(vapply(hit$refseq_accession, validate_refseq_accession, character(1L))))
  }

  refseq_ids <- lookup_refseq_orgdb(gene_symbol, species)
  if (length(refseq_ids) == 0L) {
    refseq_ids <- lookup_refseq_biomart(gene_symbol, species)
  }
  if (length(refseq_ids) == 0L) {
    stop(
      "No RefSeq accessions found for symbol ", gene_symbol, " (", species, "). ",
      "Install org.Hs.eg.db / biomaRt or add rows to ", get_refseq_cache_path(),
      call. = FALSE
    )
  }

  append_refseq_cache(species, gene_symbol, refseq_ids)
  refseq_ids
}

symbol_to_ensembl <- function(gene_symbol, species = "hsa") {
  gene_symbol <- validate_gene_symbol(gene_symbol)
  if (!species %in% names(SPECIES)) {
    stop("Unsupported species: ", species, call. = FALSE)
  }

  cache <- read_ensembl_cache()
  hit <- cache[cache$species == species & cache$gene_symbol == gene_symbol, , drop = FALSE]
  if (nrow(hit) > 0L) {
    return(validate_ensembl_gene_id(hit$ensembl_gene_id[[1L]]))
  }

  ensembl_id <- lookup_ensembl_biomart(gene_symbol, species)
  append_ensembl_cache(species, gene_symbol, ensembl_id)
  ensembl_id
}
