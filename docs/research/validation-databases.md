# Validation Databases — miRTarBase and DIANA-TarBase

**Status:** complete

**Important:** These are **experimentally supported interaction databases**, not prediction algorithms.

Use [template-algorithm-profile.md](template-algorithm-profile.md). Deep-dive checklist: [deep-dive-checklists.md](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md#validation-databases).

Research plan: [RESEARCH_PLAN.md](RESEARCH_PLAN.md#5-validation-databases)

| Question type | Prediction algorithms | Validation databases |
|---------------|----------------------|----------------------|
| Framing | Could this miRNA regulate this gene? | Has this interaction been experimentally supported? |
| Output | Numeric score (efficacy, ML, meta-integration) | Evidence type / support level (categorical) |
| Use in project | Rank predictions | Read-only overlay — never averaged into predictor scores |

---

## miRTarBase

### What is it?

**miRTarBase** is a manually curated database of **experimentally validated miRNA–target interactions (MTIs)** maintained by the Warshel Institute for Computational Biology (Chinese University of Hong Kong, Shenzhen). Launched in 2010; current release **miRTarBase 2025** (database version 10.0, September 2024) contains **>3,817,550** curated MTIs from **13,690** articles across **39 species** (Cui et al., 2025).

URL: https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025

### Evidence types supported

miRTarBase records the **experimental validation method** per MTI. Documented categories include (Cui et al., 2025; Hsu et al., 2011):

| Category | Methods (documented) |
|----------|---------------------|
| **Strong evidence** | Luciferase **reporter assay**, **Western blot**, **qPCR** / qRT-PCR |
| **Weaker evidence** | **CLIP-seq**, **PAR-CLIP**, **microarray**, pSILAC / SILAC, and other high-throughput screens |

Cui et al. (2025) explicitly state that CLIP-seq/PAR-CLIP provide **weaker experimental evidence** than reporter assays and Western blotting; users should filter by confidence level. Historical release statistics also tracked Western blot, reporter assay, pSILAC, and microarray separately (Hsu et al., 2011).

Additional curated extensions in v10 (not core MTI validation tiers): miRNA–drug associations, drug-resistance links, biomarker annotations, oxidized miRNA sequences (Cui et al., 2025).

### Support levels / confidence

miRTarBase uses a **two-tier evidence strength** model (Hsu et al., 2011; Cui et al., 2025):

| Level | Definition (documented) |
|-------|------------------------|
| **Strong** | Validated by reporter assay, Western blot, and/or qPCR |
| **Weak / less strong** | Supported primarily by high-throughput methods (microarray, pSILAC, CLIP-seq, PAR-CLIP) where direct causality is harder to establish |

Users can download subsets filtered by evidence strength (Cui et al., 2025; KG-Registry miRTarBase strong-evidence CSV product).

### How to query for gene–miRNA pairs

**Web interface (miRTarBase 2025):**

- Browse/search by **miRNA**, **target gene**, or **PubMed ID**
- Download validated MTIs (including by validated miRNA–target sites)
- Network visualization of miRNA–target relationships (Cui et al., 2025)

**R (recommended for this project):**

```r
library(multiMiR)
get_multimir(org = "hsa", mirna = "hsa-miR-21-5p", target = "PTEN",
             table = "mirtarbase", summary = TRUE)
# Or all validated DBs:
get_multimir(org = "hsa", target = "PTEN", table = "validated", summary = TRUE)
```

multiMiR `validated` table searches **miRecords**, **miRTarBase**, and **TarBase** together (Ru et al., multiMiR vignette). Returns columns including `experiment`, `support_type`, `pubmed_id`.

### Deep-dive checklist

- [x] **What miRTarBase curates and current version** — MTIs with experimental evidence; **miRTarBase 2025 / v10.0** (Cui et al., 2025)
- [x] **Evidence types and support levels** — Strong (reporter, Western, qPCR) vs weak (CLIP, microarray, pSILAC) — Hsu et al., 2011; Cui et al., 2025
- [x] **How to query for gene–miRNA pairs** — Web browse/download; multiMiR `get_multimir()` — Cui et al., 2025; multiMiR docs
- [x] **How it differs from TarBase for validation use** — See cross-cutting section below

---

## DIANA-TarBase

### What is it?

**DIANA-TarBase** is the DIANA Lab reference database of **experimentally supported miRNA targets on protein-coding transcripts**, continuously updated since 2006. Current release: **TarBase v9.0** (Skoufos et al., 2024) at https://dianalab.e-ce.uth.gr/tarbasev9/

- **~6 million** entries; **~2 million** unique miRNA–gene pairs
- **37** experimental protocols (22 low-yield + 15 high-throughput)
- **172** tissues/cell-type contexts; **24** species
- First version with **cell-type-resolved binding locations**, **TDMD** events, and **virally encoded miRNAs** (Skoufos et al., 2024)

### Evidence types supported

TarBase v9.0 content is **exclusively experimentally supported**, from two acquisition tracks (Skoufos et al., 2024; TarBase v9 portal):

| Track | Methods (documented) |
|-------|---------------------|
| **High-throughput (de novo NGS analysis)** | AGO-**HITS-CLIP**, **PAR-CLIP**, **CLASH**, **qCLASH**; miRNA transfection/knockdown **RNA-Seq**, RPF-Seq, RIP-Seq; analyzed with **microCLIP** pipeline |
| **Low-yield (manual literature curation)** | **Luciferase reporter assay**, **Western blot**, **RT-qPCR**, microarray, proteomics, biotin miRNA tagging, and others |

TarBase classifies methods as **direct** vs **indirect** experimental protocols (filterable in web UI). Direct methods include CLIP-family techniques that map binding locations; indirect methods include expression changes after miRNA perturbation (Skoufos et al., 2024).

Unique v9 content: **target-directed miRNA degradation (TDMD)** interactions (~45,000 candidate events) and **viral miRNA**–host mRNA pairs (~34,000) (Skoufos et al., 2024).

### Curation approach

1. **Manual curation** — NLP pipeline (BioBERT sentiment + Stanza dependency parsing + BERN2 NER) extracts miRNA–gene pairs from PubMed; human verification before indexing (Skoufos et al., 2024).
2. **De novo high-throughput analysis** — Uniform reprocessing of CLIP/CLASH/RNA-Seq datasets with **microCLIP** (deep-learning super-learning on CLIP peaks) and standardized RNA-Seq pipelines (Skoufos et al., 2024).
3. **Rich metadata** — Each entry links gene/transcript (Ensembl), miRNA (miRBase), sample treatment, tissue/cell type, binding coordinates (3′ UTR or CDS), publication, and regulation direction (Skoufos et al., 2024).
4. **Harmonization** — Ensembl + miRBase nomenclature; interconnected with PubMed, Ensembl genome browser, RNACentral, DIANA-microT-CDS 2023 (Skoufos et al., 2024).

### How to query for gene–miRNA pairs

**Web interface (TarBase v9):**

- Query by **miRNA name/ID** and/or **gene name/ID** simultaneously
- Filters: species, experimental method, direct/indirect protocol, cell line, tissue, cell type, gene region (3′ UTR / CDS), miRNA confidence, publication year
- Download returned interactions with full metadata (Skoufos et al., 2024)

**REST API:**

```
http://62.217.127.8/DianaTools/tarbaseApi?mirnas=hsa-let-7a-5p
http://62.217.127.8/DianaTools/tarbaseApi?genes=PTEN
```

(DIANA Tools API docs)

**R (this project):**

```r
get_multimir(org = "hsa", target = "PTEN", table = "tarbase", summary = TRUE)
```

multiMiR TarBase snapshot: **Sept 2023** (multiMiR database listing); for latest v9 metadata use TarBase web/API directly.

### Deep-dive checklist

- [x] **What TarBase curates and current version** — Experimentally supported miRNA–gene interactions; **TarBase v9.0** (Skoufos et al., 2024)
- [x] **Evidence types and curation approach** — 37 protocols; manual curation + de novo NGS (microCLIP) — Skoufos et al., 2024
- [x] **How to query for gene–miRNA pairs** — TarBase v9 web + REST API + multiMiR — Skoufos et al., 2024; DIANA API
- [x] **How it differs from miRTarBase for validation use** — See cross-cutting section below

---

## Cross-cutting: prediction vs validation

| Aspect | Predictors (TargetScan, DIANA-microT-CDS, miRDB, mirDIP) | Validation DBs (miRTarBase, TarBase) |
|--------|----------------------------------------------------------|--------------------------------------|
| Question | Could this pair be functional? | Was this pair experimentally supported? |
| Output | Numeric score | Evidence method + metadata (categorical) |
| Genome coverage | Comprehensive predictions | Sparse — only tested pairs |
| Negative evidence | Not provided | Absence of entry ≠ proven non-target |
| Project use | Primary ranking layer | **Read-only overlay** on predictions |

**Never** call validation databases "prediction algorithms" or merge their evidence into consensus predictor scores.

---

## How miRTarBase and TarBase differ for validation use

| Dimension | miRTarBase 2025 | DIANA-TarBase v9.0 |
|-----------|-----------------|---------------------|
| Scale | ~3.8M MTIs | ~6M entries (~2M unique pairs) — Skoufos et al., 2024 cite 3× miRTarBase |
| Evidence tiers | **Strong vs weak** (method-based) | **High-throughput vs low-yield** + direct vs indirect protocol filters |
| Binding site detail | Some target-site annotations; CLIP-derived MTIs | **Millions of binding locations** at cell-type resolution (Skoufos et al., 2024) |
| Context metadata | Disease, drug, biomarker extensions | Tissue, cell type, cell line, experimental condition per entry |
| Unique content | miRNA oxidation, drug resistance, biomarkers | TDMD events, viral miRNAs, microCLIP reanalysis |
| Curation | Manual + LLAMA3 NLP pipeline (Cui et al., 2025) | Manual NLP + uniform de novo CLIP/CLASH/RNA-Seq pipelines |
| Best for overlay | Quick strong-evidence filter (reporter/Western/qPCR) | Context-specific validation (cell type, binding region, method) |

**Practical recommendation:** Use **both** as complementary overlays — miRTarBase for straightforward strong-evidence flags; TarBase when cell-type, binding location, or method specificity matters. multiMiR `table = "validated"` queries both (plus miRecords) in one R call.

---

## 4b. Score inventory (validation)

Validation databases do **not** produce prediction scores. Document evidence categories instead:

| Resource | Metric / category | Level | Evidence type | Range / values | Direction | Notes |
|----------|-------------------|-------|---------------|----------------|----------|-------|
| miRTarBase | Support level | gene–miRNA | experimental support | Strong / Weak (2 tiers) | Strong > Weak | Based on validation method — Hsu et al., 2011; Cui et al., 2025 |
| miRTarBase | Validation method | gene–miRNA | experimental support | Categorical (reporter, Western, qPCR, CLIP, microarray, etc.) | Reporter/Western/qPCR = strong | Per-interaction metadata |
| DIANA-TarBase | Experimental method | gene–miRNA | experimental support | 37 protocols (22 low-yield + 15 high-throughput) | Low-yield generally stronger functional evidence | Skoufos et al., 2024 |
| DIANA-TarBase | Protocol type | gene–miRNA | experimental support | Direct / Indirect | Direct = binding-mapping (CLIP family) | Web filter — Skoufos et al., 2024 |
| DIANA-TarBase | Target region | site / gene–miRNA | experimental support | 3′ UTR / CDS | N/A | Skoufos et al., 2024 |

**Not comparable** to predictor numeric scores — use as binary/tiered overlay only.

---

## Integration notes (for consensus project)

### Read-only overlay design

For each predicted gene–miRNA pair from TargetScan/DIANA/miRDB/mirDIP, attach validation flags **without** modifying predictor ranks:

| Field | Source | Logic |
|-------|--------|-------|
| `experimental_support` | miRTarBase and/or TarBase | `TRUE` if either DB reports the pair |
| `support_level` | miRTarBase | `strong` / `weak` / `NA` |
| `validation_methods` | Both | Concatenate distinct method strings |
| `validation_sources` | Both | e.g. `miRTarBase`, `TarBase`, or both |
| `binding_site_reported` | TarBase | `TRUE` if TarBase provides coordinates |

### R module (`R/module_validation.R`)

- Primary interface: **multiMiR** `get_multimir(table = "validated")` for unified lookup
- Optional: DIANA TarBase REST API for rich metadata (cell type, binding location)
- **Do not** use validation presence as a numeric score in `integrate_scores.R`
- Tiered display: predictions with `support_level == "strong"` highlighted; weak/CLIP-only annotated separately

### Why validation does not resolve predictor disagreement

Experimental databases confirm **that** an interaction was observed under specific conditions — they do not adjudicate **which predictor score** is correct. A pair may be predicted by TargetScan but not miRDB, yet still appear in TarBase from a CLIP experiment in a particular cell line — or absent from both despite high computational scores.

---

## Sources

- Cui S, Yu S, Huang H-Y, Lin Y-C-D, Huang Y, Zhang B, et al.; Huang H-D (2025). miRTarBase 2025: updates to the collection of experimentally validated microRNA–target interactions. *Nucleic Acids Res* 53(D1):D147–D156. doi:10.1093/nar/gkae1072
- Hsu S-D, Lin F-M, Wu W-Y, Liang C, Huang W-C, Chan W-L, et al.; Huang H-D (2011). miRTarBase: a database curates experimentally validated microRNA–target interactions. *Nucleic Acids Res* 39:D163–D169. doi:10.1093/nar/gkq1107
- Skoufos G, Kakoulidis P, Tastsoglou S, Zacharopoulou E, Kotsira V, Miliotis M, et al.; Hatzigeorgiou AG (2024). TarBase-v9.0 extends experimentally supported miRNA–gene interactions to cell-types and virally encoded miRNAs. *Nucleic Acids Res* 52(D1):D368–D377. doi:10.1093/nar/gkad1071
- TarBase v9.0 portal: https://dianalab.e-ce.uth.gr/tarbasev9/
- DIANA Tools TarBase REST API: http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=site%2Fhelp&topic=api
- Ru Y, Mulvahill MJ, Mahaffey S, Kechris K (2025). multiMiR: Integration of multiple microRNA-target databases. Bioconductor package. https://bioconductor.org/packages/multiMiR/

---

## Integration notes (for consensus project)

See overlay design above. Validation is **Phase 4 annotation**, not Phase 4 scoring.
