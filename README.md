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
| `scripts/` | Bulk data download and per-gene export CLIs |
| `data/raw/` | Bulk prediction downloads (gitignored) |
| `data/cache/` | Per-gene exports, API caches, symbol maps (gitignored) |
| `.cursor/skills/` | Agent skills for research, R modules, presentation |
| `.cursor/rules/` | Persistent project constraints |

## Prerequisites

Run all commands below from the **project root** (the directory that contains `R/` and `scripts/`).

**R packages**

| Package | Used for |
|---------|----------|
| `httr2` | DIANA and mirDIP live API queries |
| `data.table` | Optional; faster reads of large bulk files |
| `dplyr` | Optional; binding columns when integrating sources |

```r
install.packages(c("httr2", "data.table", "dplyr"))
```

Scripts create `data/raw/` and `data/cache/` automatically if they are missing.

## Data not in Git

`data/raw/` and `data/cache/` are listed in `.gitignore`. After cloning the repo you need to **download bulk files** (where a predictor uses them) and **run export scripts** to generate per-gene outputs. Nothing under `data/` is required to browse docs or run unit tests that use small fixtures.

### Bulk downloads (one-time per machine)

Re-run with `--force` to replace existing files.

| Script | Writes to | Required for | Notes |
|--------|-----------|--------------|-------|
| `scripts/download_targetscan.R` | `data/raw/targetscan/` | TargetScan exports | Human (`vert_80`) + mouse (`mmu_80`); ~650 MB total |
| `scripts/download_mirdb.R` | `data/raw/mirdb/` | miRDB exports | Single large gzipped prediction file |
| `scripts/download_diana_bulk.R` | `data/raw/diana/` | DIANA bulk fallback | Optional if the DIANA API is reachable; default fetches `hsa` and `mmu` |

```bash
Rscript scripts/download_targetscan.R
Rscript scripts/download_mirdb.R
Rscript scripts/download_diana_bulk.R

# Optional: only human DIANA bulk
Rscript scripts/download_diana_bulk.R hsa

# Re-download after upstream updates
Rscript scripts/download_targetscan.R --force
```

**mirDIP** has no bulk download script in this repo; exports call the public API and cache responses under `data/cache/mirdip/`.

### Per-gene prediction exports

Each export script queries one predictor for a **gene symbol** and species (`hsa` default, or `mmu`), then writes:

- `data/cache/<GENE>_<species>_<source>.csv`
- `data/cache/<GENE>_<species>_<source>.rds`

Replace `PTEN` with any validated HUGO symbol. Pass `--help` (or run with no args) on any script for flags.

**TargetScan** (requires TargetScan bulk download):

```bash
Rscript scripts/export_score_targetscan.R PTEN hsa --mirna-prefix=species
```

Common flags: `--predictions=all|default`, `--score-name=cwcs|biochem_occupancy`, `--min-cwcs=<numeric>` (keep rows with CWCS ≤ value; more negative = stronger).

**DIANA-microT-CDS** (API by default; bulk fallback if `download_diana_bulk.R` was run):

```bash
Rscript scripts/export_score_diana.R PTEN hsa
```

Common flags: `--threshold=0.7`, `--prefer-api=true|false`, `--mirna-resource=mirbase|mirgenedb`.

**miRDB** (requires miRDB bulk download):

```bash
Rscript scripts/export_score_mirdb.R PTEN hsa
```

Common flag: `--min-score=60` (default mining cutoff in `R/config.R`).

**mirDIP** (network + `httr2`):

```bash
Rscript scripts/export_score_mirdip.R PTEN hsa
```

Common flags: `--score-class=Very High|High|Medium|Low`, `--min-score=<numeric>`.

**Run all four for one gene:**

```bash
Rscript scripts/export_score_targetscan.R PTEN hsa --mirna-prefix=species
Rscript scripts/export_score_diana.R PTEN hsa
Rscript scripts/export_score_mirdb.R PTEN hsa
Rscript scripts/export_score_mirdip.R PTEN hsa
```

### Integrated consensus export

Combines predictors with `integrate_scores()` (default method: **tiered** — not a raw average of incompatible scores). Writes:

| File | Contents |
|------|----------|
| `data/cache/<GENE>_<species>_integrated.csv` / `.rds` | Ranked integrated table |
| `data/cache/<GENE>_<species>_predictions_combined.csv` | Long-format rows from all sources (unless disabled) |

**Live queries** (hits APIs / bulk as needed; does not require prior per-source exports):

```bash
Rscript scripts/export_score_integrated.R PTEN hsa --method=tiered
```

**From cached per-source files** (run the four `export_score_*.R` scripts first; skips re-querying sources that already have cache files):

```bash
Rscript scripts/export_score_integrated.R PTEN hsa --use-cache=true --method=tiered
```

Other flags: `--sources=targetscan,diana,mirdb,mirdip` (comma-separated subset), `--write-combined=true|false`, `--method=tiered|rank_fusion|percentile`.

### End-to-end example (PTEN, human)

```bash
# 1) One-time bulk data
Rscript scripts/download_targetscan.R
Rscript scripts/download_mirdb.R
Rscript scripts/download_diana_bulk.R

# 2) Per-source caches (optional but useful for inspection and offline integration)
Rscript scripts/export_score_targetscan.R PTEN hsa --mirna-prefix=species
Rscript scripts/export_score_diana.R PTEN hsa
Rscript scripts/export_score_mirdb.R PTEN hsa
Rscript scripts/export_score_mirdip.R PTEN hsa

# 3) Integrated ranking
Rscript scripts/export_score_integrated.R PTEN hsa --use-cache=true --method=tiered
```

### Interactive R (same logic as the CLIs)

```r
source("R/main.R")

find_mirna_regulators("PTEN", sources = "targetscan")
find_mirna_regulators(
  "PTEN",
  integrate = TRUE,
  integrate_method = "tiered",
  use_cached_predictions = TRUE
)
```

TargetScan-only options: `query_targets_targetscan("PTEN", predictions = "all", score_name = "cwcs")`, `min_cwcs = -0.5`, `mirna_prefix = "species"` for human mature miRNAs after family expansion.

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
