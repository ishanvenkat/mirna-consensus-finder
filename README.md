# miRNA Target-Prediction Consensus Tool

Gene-in → ranked miRNA-out tool that queries multiple miRNA target-prediction resources and prepares for consensus scoring validated against experimental databases.

**Example:** Enter `PTEN` → ranked list of miRNAs predicted to regulate PTEN with per-source scores.

## Architecture

```
Gene → TargetScan + DIANA + miRDB + mirDIP → normalization → consensus score → validation (miRTarBase/TarBase)
```

## Project Layout

| Path | Purpose |
|------|---------|
| `docs/research/` | Algorithm profiles and score comparison matrix |
| `docs/presentation/` | Slideshow content and bibliography |
| `R/` | Modular R predictors and integration |
| `scripts/` | Bulk data download helpers |
| `data/raw/` | Bulk prediction downloads (gitignored) |
| `data/cache/` | Per-gene API query cache |
| `.cursor/skills/` | Agent skills for research, R modules, presentation |
| `.cursor/rules/` | Persistent project constraints |

## Quick Start

Download TargetScan 8.0 bulk data (one-time, ~650 MB total for human + mouse):

```bash
Rscript scripts/download_targetscan.R
```

Download miRDB v6.0 bulk predictions (one-time, large file):

```bash
Rscript scripts/download_mirdb.R
```

Optional R packages: `data.table` (faster reads), `dplyr` (multi-source column binding).

Required for DIANA/mirDIP API access (bulk fallbacks work without it for TargetScan/miRDB):

```r
install.packages("httr2")
```

```r
source("R/main.R")

find_mirna_regulators("PTEN", sources = "targetscan")
# TargetScan-only query with conservation annotation:
# query_targets_targetscan("PTEN", predictions = "all", score_name = "cwcs")  # CWCS: more negative = stronger
# query_targets_targetscan("PTEN", min_cwcs = -0.5)  # keep CWCS <= -0.5 (stronger repression)
# Human mature miRNAs only (after family expansion from miR_Family_Info.txt):
# query_targets_targetscan("PTEN", mirna_prefix = "species")
```

Export per-source predictions (live query → cache):

```bash
Rscript scripts/export_score_targetscan.R PTEN hsa --mirna-prefix=species
Rscript scripts/export_score_diana.R PTEN hsa
Rscript scripts/export_score_mirdb.R PTEN hsa
Rscript scripts/export_score_mirdip.R PTEN hsa
# → data/cache/<GENE>_<species>_<source>.csv and .rds
```

Integrate via `find_mirna_regulators()` (live queries by default, or read cache with fallback):

```bash
Rscript scripts/export_score_integrated.R PTEN hsa --method=tiered
Rscript scripts/export_score_integrated.R PTEN hsa --use-cache=true --method=tiered
# → data/cache/PTEN_hsa_integrated.csv|.rds
# → data/cache/PTEN_hsa_predictions_combined.csv (default)
```

In R, `find_mirna_regulators("PTEN", integrate = TRUE, use_cached_predictions = TRUE)` uses cached per-source files when present.

## Status

| Phase | Status |
|-------|--------|
| 0 — Skills & rules | Complete |
| 1 — Scaffolding | Complete |
| 2 — Algorithm research | Complete — see `docs/research/` |
| 3 — Slideshow | Complete — see `docs/presentation/` |
| 4 — R modules | TargetScan, DIANA, miRDB, mirDIP implemented; `integrate_scores()` wired in `main.R` |
| 5 — Validation demo | `query_validation()` via multiMiR (fixture-backed tests) |

## Key Principle

**Do not blindly average scores** from different tools. Each score has different semantics — see `docs/research/score-comparison-matrix.md`.

## Prediction vs Validation

- **Predictors:** TargetScan, DIANA-microT-CDS, miRDB, mirDIP
- **Validation:** miRTarBase, DIANA-TarBase (experimental evidence only)

## Reference

Project context: [`miRNA project.pdf`](miRNA%20project.pdf)
