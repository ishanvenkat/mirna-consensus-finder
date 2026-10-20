---
name: mirna-r-modules
description: >-
  Implement modular R predictors for miRNA target lookup with a shared
  query_targets interface. Use when writing R/module_*.R, scripts/download_*.R,
  integrate_scores.R, querying TargetScan/DIANA/mirDIP/miRDB APIs or bulk data,
  or running the PTEN smoke test.
---

# miRNA R Modules

## When to Use

Apply when implementing or extending prediction modules, download scripts, score integration, or the main entry point.

## Module Contract

Every predictor module exports:

```r
query_targets <- function(gene_symbol, species = "hsa", ...) { ... }
```

Return a tibble with columns:

| Column | Type | Description |
|--------|------|-------------|
| `gene_symbol` | character | HUGO symbol (validated) |
| `mirna` | character | miRNA name (miRBase format) |
| `score_raw` | numeric | Score as returned by source |
| `score_name` | character | e.g. `"context++"`, `"microT"`, `"integrated"` |
| `score_direction` | character | `"higher_is_stronger"` or `"lower_is_stronger"` |
| `site_count` | integer | Number of sites if applicable; `NA` otherwise |
| `source` | character | `"targetscan"`, `"diana"`, `"mirdb"`, `"mirdip"` |
| `retrieved_at` | POSIXct | Timestamp of query |

See [module-interface.md](module-interface.md) for examples and merge pattern.

## File Layout

| File | Role |
|------|------|
| `R/config.R` | Species, paths, thresholds — single source of truth |
| `R/utils_gene_id.R` | Symbol ↔ Ensembl mapping (DIANA requires Ensembl) |
| `R/module_<source>.R` | One file per data source |
| `R/module_validation.R` | miRTarBase/TarBase overlay (read-only) |
| `R/integrate_scores.R` | Normalization + fusion (after research complete) |
| `R/main.R` | `find_mirna_regulators()` entry point |

## Data Access

| Source | Method | Fallback |
|--------|--------|----------|
| TargetScan | Bulk TSV in `data/raw/` + local filter | Re-download via `scripts/download_targetscan.R` |
| DIANA | REST API via `httr2` | Bulk gz from DIANA webserver |
| mirDIP | mirDIP API; cache in `data/cache/` | — |
| miRDB | Bulk download + local filter | — |
| Validation | `multiMiR` or TarBase REST | — |

**Never scrape HTML result tables.**

## Integration Constraints

- Do **not** implement simple score averaging until `docs/research/score-comparison-matrix.md` is complete
- Prefer rank fusion or percentile normalization within source
- mirDIP is a **meta-layer** — do not double-count TargetScan/DIANA/miRDB inside consensus weights
- Validation evidence is an **overlay flag** — never mixed into prediction scores

## Testing

- Smoke test each module with `"PTEN"` before marking done
- Unit tests in `R/tests/testthat/` using cached fixtures in `data/cache/`
- Validate gene symbols before API calls or file filters

## Security

- Never pass user gene input to shell commands or unvalidated file paths
- Use allowlisted species codes and sanitized symbols only

## Additional Resources

- Interface examples: [module-interface.md](module-interface.md)
- Config defaults: [`R/config.R`](../../../R/config.R)
