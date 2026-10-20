# Shared configuration for miRNA target prediction modules.
# All paths and thresholds live here — do not hardcode in module files.

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = FALSE)

PATHS <- list(
  raw     = file.path(PROJECT_ROOT, "data", "raw"),
  cache   = file.path(PROJECT_ROOT, "data", "cache"),
  targetscan = file.path(PROJECT_ROOT, "data", "raw", "targetscan"),
  diana   = file.path(PROJECT_ROOT, "data", "raw", "diana"),
  mirdb   = file.path(PROJECT_ROOT, "data", "raw", "mirdb")
)

SPECIES <- list(
  hsa = list(
    label = "Homo sapiens",
    targetscan_db = "vert_80",
    diana_species = "human",
    mirna_prefix = "hsa",
    ncbi_tax_id = 9606L
  ),
  mmu = list(
    label = "Mus musculus",
    targetscan_db = "mmu_80",
    diana_species = "mouse",
    mirna_prefix = "mmu",
    ncbi_tax_id = 10090L
  )
)

THRESHOLDS <- list(
  diana_microt = 0.7,
  mirdb = 60,         # default mining cutoff (miRDB FAQ); display minimum is 50
  mirdip_min_score = NULL,
  targetscan_min_cwcs = NULL  # optional filter: keep rows with CWCS <= this value (more negative = stronger)
)

TARGETSCAN <- list(
  version = "8.0",
  files = list(
    summary_all = "Summary_Counts.all_predictions.txt",
    summary_default = "Summary_Counts.default_predictions.txt",
    gene_info = "Gene_info.txt",
    family_info = "miR_Family_Info.txt"
  ),
  zip_files = list(
    summary_all = "Summary_Counts.all_predictions.txt.zip",
    summary_default = "Summary_Counts.default_predictions.txt.zip",
    gene_info = "Gene_info.txt.zip",
    family_info = "miR_Family_Info.txt.zip"
  ),
  base_url = list(
    vert_80 = "https://www.targetscan.org/vert_80/vert_80_data_download",
    mmu_80 = "https://www.targetscan.org/mmu_80/mmu_80_data_download"
  ),
  cols = list(
    gene_symbol = "Gene Symbol",
    mirna_family = "miRNA family",
    representative_mirna = "Representative miRNA",
    transcript_id = "Transcript ID",
    cwcs = "Cumulative weighted context++ score",
    aggregate_pct = "Aggregate PCT",
    total_context = "Total context++ score",
    biochem_occupancy = "Predicted occupancy (transfected miRNA)",
    conserved_sites = "Total num conserved sites",
    nonconserved_sites = "Total num nonconserved sites",
    species_id = "Species ID"
  ),
  family_info_cols = list(
    family = "miR family",
    seed = "Seed+m8",
    species_id = "Species ID",
    mirna = "MiRBase ID"
  ),
  gene_info_cols = list(
    transcript_id = "Transcript ID",
    gene_symbol = "Gene symbol",
    representative = "Representative transcript?"
  )
)

get_targetscan_root <- function() {
  override <- getOption("mirna.targetscan.data_root")
  if (!is.null(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  PATHS$targetscan
}

DIANA <- list(
  version = "2023",
  downloads_base = paste0(
    "https://dianalab.e-ce.uth.gr/microt_webserver_downloads/flat_files"
  ),
  bulk_files = list(
    hsa = list(
      mirbase = "interactions_human.microT.mirbase.txt.gz",
      mirgenedb = "interactions_human.microT.mirgenedb.txt.gz"
    ),
    mmu = list(
      mirbase = "interactions_mouse.microT.mirbase.txt.gz",
      mirgenedb = "interactions_mouse.microT.mirgenedb.txt.gz"
    )
  ),
  default_mirna_resource = "mirbase",
  mirna_resource_label = "miRBase v22.1",
  microt_version_label = "microT-CDS 2022",
  cols = list(
    mirna = "mirna",
    ensembl_gene_id = "ensembl_gene_id",
    interaction_score = "interaction_score"
  ),
  score_name = "microT",
  ensembl_cache_file = "ensembl_symbol_map.tsv"
)

get_diana_root <- function() {
  override <- getOption("mirna.diana.data_root")
  if (!is.null(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  PATHS$diana
}

MIRDB <- list(
  version = "6.0",
  download_url = "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
  bulk_file = "miRDB_v6.0_prediction_result.txt.gz",
  target_id_type = "refseq",
  cols = c("mirna", "refseq_accession", "score"),
  score_name = "MirTarget",
  display_min_score = 50,
  refseq_cache_file = "refseq_symbol_map.tsv"
)

get_refseq_cache_path <- function() {
  override <- getOption("mirna.refseq.cache_path")
  if (!is.null(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  file.path(PATHS$cache, MIRDB$refseq_cache_file)
}

get_mirdb_root <- function() {
  override <- getOption("mirna.mirdb.data_root")
  if (!is.null(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  PATHS$mirdb
}

MIRDIP <- list(
  api_path = "/mirDIP/Http_U",
  score_class_codes = c(
    "Very High" = "0",
    "High" = "1",
    "Medium" = "2",
    "Low" = "3"
  ),
  default_score_class = "Low",
  score_name = "integrated",
  human_only = TRUE
)

get_ensembl_cache_path <- function() {
  override <- getOption("mirna.ensembl.cache_path")
  if (!is.null(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  file.path(PATHS$cache, DIANA$ensembl_cache_file)
}

API <- list(
  diana_microt_cds = "http://62.217.127.8/DianaTools/microT_CDSApi",
  diana_microt_webserver = paste0(
    "https://dianalab.e-ce.uth.gr/microt_webserver_backend/api"
  ),
  mirdip_base = "https://ophid.utoronto.ca"
)

ensure_data_dirs <- function() {
  for (path in unlist(PATHS)) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }
}
