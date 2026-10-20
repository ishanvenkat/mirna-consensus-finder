# TargetScan

**Status:** complete

---

## 1. What is it?

### Who developed it?

TargetScan was developed by the **Bartel lab** (Whitehead Institute for Biomedical Research / Howard Hughes Medical Institute). Key contributors across releases include David P. Bartel, Benjamin P. Lewis, Alexander Grimson, Robin C. Friedman, Vikram Agarwal, Sean E. McGeary, and colleagues.

### When was it introduced?

- **2005** — seed-based target prediction (Lewis et al., 2005)
- **2007** — context scores and site-type efficacy (Grimson et al., 2007)
- **2009** — preferential conservation / PCT (Friedman et al., 2009)
- **2011** — context+ with TA and SPS (García et al., 2011)
- **2015** — Context++ and TargetScan 7.0 (Agarwal et al., 2015)
- **2020** — TargetScan 8.0 default ranking via biochemical model (McGeary et al., 2019)

### Current / relevant version

**TargetScan 8.0** (January 2020 ranking update). Mammalian portals:

- Human: [TargetScanHuman 8.0](https://www.targetscan.org/vert_80/)
- Mouse: [TargetScanMouse 8.0](https://www.targetscan.org/mmu_80/)

Release 8 uses the **same 3′ UTR profiles and miRNA families as Release 7**; the main mammalian change is the **default ranking method** (McGeary et al., 2019), with Context++/PCT retained as options.

### What problem does it solve?

TargetScan predicts which mRNAs are regulated by which miRNAs by scanning 3′ UTRs (and some ORF regions) for **canonical seed-matched sites**, scoring site efficacy and evolutionary conservation, and ranking gene–miRNA pairs for hypothesis generation. It separates **predicted biochemical repression** from **evolutionary conservation evidence**, allowing users to prioritize either cellular response or biological consequence.

---

## 2. How does it work?

```
miRNA families (shared seed) → scan representative 3′ UTR profile
→ classify site type (8mer / 7mer-m8 / 7mer-A1 / 6mer; optional 3′ compensatory)
→ compute conserved branch length & PCT (conservation track)
→ compute Context++ per site (efficacy track, linear model)
→ aggregate to gene–miRNA (total context++, CWCS, aggregate PCT)
→ rank targets (8.0 default: biochemical predicted occupancy; options: CWCS or aggregate PCT)
```

### Step-by-step process

1. **Curate miRNA families** — group mature miRNAs by identical seed (positions 2–8); classify family conservation (broadly conserved, conserved, poorly conserved).
2. **Define 3′ UTR profiles** — Gencode annotations extended with RefSeq, RNA-seq, and 3P-seq tags; **representative transcript** = variant with most 3P-seq tags (TargetScan 7–8).
3. **Identify candidate sites** — search representative 3′ UTRs for perfect seed matches; orthologous UTRs from UCSC multiz alignments enable conservation analysis.
4. **Score each site** — branch length → conserved/nonconserved; PCT (probability of conserved targeting); Context++ (14-feature efficacy model).
5. **Aggregate per gene–miRNA family** — resolve overlapping sites; compute total context++ score, cumulative weighted context++ score (CWCS), and aggregate PCT.
6. **Rank predictions** — TargetScan 8.0 defaults to **predicted occupancy** from the McGeary biochemical model; users can switch to CWCS or aggregate PCT ranking.

### Candidate site identification

- **Canonical sites:** Watson–Crick pairing to miRNA seed (positions 2–7 or 2–8), classified as 8mer, 7mer-m8, 7mer-A1, or 6mer (see Section 8).
- **3′ compensatory sites:** imperfect seed match compensated by conserved 3′ pairing (Friedman et al., 2009); negative site-type codes in download files.
- **Centered sites:** included in TargetScan 7; **removed in TargetScan 8.0** (McGeary et al., 2019 showed they are not functional for most miRNAs).
- **6mer sites:** included in human TargetScan 7+ predictions; always classified **nonconserved** in TargetScanHuman 7+.
- **Conservation filter:** site conservation uses **conserved branch length** thresholds that differ by site type and release (FAQs).

### Scoring and ranking

TargetScan provides **three independent ranking axes** for mammalian predictions:

| Rank by | Level | Evidence | Default in |
|---------|-------|----------|------------|
| Predicted occupancy (biochemical model) | gene–miRNA | efficacy | **TargetScan 8.0 web** |
| Cumulative weighted context++ score (CWCS) | gene–miRNA | efficacy | TargetScan 7.x web; 8.0 option |
| Aggregate PCT | gene–miRNA | conservation | Optional in 7.x and 8.0 |

TargetScan explicitly recommends choosing thresholds at the **target (gene–miRNA) level** rather than individual sites, because multiple weak sites can sum to strong repression (FAQs).

---

## 3. What biological/computational features does it use?

| Feature | Used? | Source / notes |
|---------|-------|----------------|
| Seed matching | ✅ | Lewis et al., 2005; positions 2–7/8 |
| 8mer / 7mer-m8 / 7mer-A1 / 6mer site types | ✅ | Lewis 2005; Friedman 2009; FAQs site-type column |
| Evolutionary conservation | ✅ | Branch length + PCT (Friedman et al., 2009) |
| Binding / free energy | Partial | SPS (seed-pairing stability) in Context++; McGeary model uses measured/predicted *K*<sub>D</sub>, not ΔG alone |
| RNA / site accessibility | ✅ | SA feature in Context++ (RNAplfold; Agarwal et al., 2015) |
| Local AU content | ✅ | Context++ feature (Grimson et al., 2007) |
| Target-site abundance | ✅ | TA feature (García et al., 2011) |
| Position within transcript | ✅ | Minimum distance from 3′ UTR ends (Context++) |
| 3′ compensatory pairing | ✅ | Supplementary pairing feature; separate 3′ compensatory site class |
| G:U pairing | Partial | sRNA/site position-8 identity features in Context++; not primary seed rule |
| Machine learning | ✅ (8.0 default) | CNN predicts *K*<sub>D</sub> for McGeary biochemical model; Context++ itself is **multiple linear regression**, not ML |
| Other sequence/context features | ✅ | 3′ UTR length, ORF length, ORF 8mer count, offset-6mer count, 3P-seq AIR weighting |

---

## 4. What does its score mean?

| Question | Answer |
|----------|--------|
| Score name(s) | Context++, total context++ score, cumulative weighted context++ score (CWCS), PCT, predicted occupancy (biochemical model), predicted relative *K*<sub>D</sub> |
| Numerical range | Context++: site-type-bounded (8mer ≤ −0.03, 7mer-m8 ≤ −0.02, 7mer-A1 ≤ −0.01, 6mer ≤ 0). CWCS: log<sub>2</sub> predicted fold-change in target expression (typically ≤ 0 in 8.0 downloads; weak near 0). PCT: 0–1. Predicted occupancy: fraction of AGO–miRNA complex on target sites (0–1 range). |
| Higher or lower = stronger? | Context++: **more negative** = stronger predicted site repression. CWCS (bulk **Cumulative weighted context++ score**): **more negative** = stronger predicted gene-level repression (same log-scale convention as Context++). Biochemical predicted occupancy: **higher** = stronger. PCT: **higher** = stronger conservation evidence. |
| Score type | Context++: linear regression contribution (log-scale repression term). CWCS: log<sub>2</sub> cumulative repression. PCT: probability. Biochemical: predicted occupancy fraction. |
| Absolute interpretation across genes? | Context++ and CWCS are calibrated on experimental compendia but **not guaranteed comparable across miRNA families** without percentile ranks. PCT is comparable as probability. Occupancy depends on assumed miRNA concentration scenario. |
| Already normalized? | Context++ features are trimmed-normalized to (0,1) before regression (Agarwal et al., 2015). Percentile columns provided in downloads. PCT is inherently normalized 0–1. |
| Recommended cutoff / threshold | **No universal cutoff** — FAQs defer to user sensitivity/specificity trade-offs; thresholds differ by miRNA family and biological question. |
| Site scores aggregated to gene–miRNA score? | Yes — total context++ (sum of non-overlapping sites), CWCS (AIR-weighted marginal repression), aggregate PCT (gene-level conservation rollup). |

### 4b. Score inventory

| Score name | Level | Evidence type | Range | Direction | Notes |
|------------|-------|---------------|-------|-----------|-------|
| Context++ | site | efficacy | Site-type bounded; 8mer ≤ −0.03 | More negative = stronger | Sum of 14 feature contributions; separate linear model per site type |
| Total context++ score | gene–miRNA | efficacy | Unbounded (sum of site scores) | More negative sum = stronger | Sum over **non-overlapping** sites maximizing total |
| Weighted Context++ (CWCS) | gene–miRNA | efficacy | log<sub>2</sub>; typically ≤ 0 in 8.0 Summary Counts | **More negative** = stronger | `CWCS = log₂(1 − Cₙ)` from AIR-weighted marginal repression |
| PCT | site; aggregate at gene–miRNA | conservation | 0–1 | Higher = stronger | Probability of preferential conservation; **separate ranking axis from Context++** |
| Biochemical model score (predicted occupancy) | gene–miRNA | efficacy | ~0–1 (fraction) | Higher = stronger | **TargetScan 8.0 default ranking**; McGeary CNN + biochemical model; three concentration scenarios in downloads |

### Score interpretation notes

- **Context++** = predicted targeting efficacy at a single site (biochemical response to miRNA).
- **PCT** = probability that a site is conserved **because** it is targeted (evolutionary evidence of biological consequence).
- **PCT is also one of the 14 Context++ features**, but the standalone PCT score and aggregate PCT ranking remain **conservation readouts** — do not treat Context++ and PCT as interchangeable outputs.
- **CWCS vs biochemical occupancy:** both predict efficacy; biochemical model (8.0 default) outperformed retrained Context++ in McGeary et al., 2019 cross-validation, but uses a different mathematical framework (*K*<sub>D</sub> + occupancy → repression).
- **No recommended cutoffs** — use target-level thresholds based on your sensitivity/specificity needs (FAQs).

---

## 5. Advantages

- Extensively calibrated against large experimental compendia (74 transfection datasets for Context++; AGO-RBNS + transfection data for biochemical model).
- Separates **efficacy** (Context++/CWCS/biochemical) from **conservation** (PCT), avoiding conflation of biochemical and evolutionary evidence.
- Transparent **site-level output** with site type, context++ score, branch length, and PCT per site.
- Long publication lineage with downloadable parameters, Perl scripts, and bulk files for reproducible offline analysis.
- Representative-transcript and 3P-seq AIR weighting reduces isoform confounding (Nam et al., 2014; Agarwal et al., 2015).
- TargetScan 8.0 retains legacy Context++/PCT ranking as explicit options while improving default efficacy ranking.

---

## 6. Disadvantages

- **3′ UTR–centric** — standard gene pages emphasize 3′ UTR sites; ORF sites exist but are secondary.
- **Representative transcript choice** affects which sites appear; genes with complex isoforms may be under- or over-represented.
- **Multiple incompatible ranking scores** (occupancy vs CWCS vs aggregate PCT) — 8.0 default differs from 7.x, creating integration hazards across tools and releases.
- **No universal score threshold** — users must define cutoffs per study.
- **6mer-only targets** are scored but not displayed on gene pages; targets with only 6mer sites are not listed (Agarwal et al., 2015).
- **McGeary CNN *K*<sub>D</sub> predictions** extrapolate from six experimentally measured miRNAs — accuracy may drop for miRNAs with atypical seed-binding patterns.
- Poorly conserved miRNA families include misannotated small RNAs (FAQs).

---

## 7. When would you use it?

Use TargetScan when you need **seed-based mammalian target predictions** with explicit separation of predicted repression and evolutionary conservation — for example, ranking PTEN targets of miR-21 by CWCS while flagging high aggregate PCT sites as evolutionarily supported.

- Choose **predicted occupancy** (8.0 default) for best-validated cellular repression ranking.
- Choose **CWCS** for compatibility with TargetScan 7 literature and linear Context++ framework.
- Choose **aggregate PCT** when biological consequence / conserved targeting matters more than biochemical potency (especially for broadly conserved miRNA families).
- Always confirm **miRNA and mRNA co-expression** in your cell type — TargetScan does not filter by expression (FAQs).

---

## 8. Method-specific deep dive

- [x] **Seed matching** — Canonical sites require Watson–Crick pairing to miRNA positions 2–7 (6mer, 7mer-A1) or 2–8 (7mer-m8, 8mer). miRNA families grouped by identical seed. (Lewis et al., 2005; FAQs)
- [x] **8mer sites** — Exact match to positions 2–8 of miRNA **plus** an adenosine at target position 1 (immediately 3′ of the seed match). Strongest canonical site type. Download code: `3`. (FAQs; Lewis et al., 2005; Friedman et al., 2009)
- [x] **7mer-m8 sites** — Exact match to positions 2–8 (seed + position 8). Second-strongest canonical type. Download code: `2`. (FAQs; Friedman et al., 2009)
- [x] **7mer-A1 sites** — Exact match to positions 2–7 plus adenosine opposite miRNA position 1. Download code: `1`. (FAQs; Lewis et al., 2005)
- [x] **6mer sites** — Exact match to positions 2–7 only. Weakest canonical type; always nonconserved in TargetScanHuman 7+. Download code: `4`. (FAQs; Friedman et al., 2009; Agarwal et al., 2015)
- [x] **Relative site efficacy** — Experimental hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer (Friedman et al., 2009; Grimson et al., 2007 site-type weights 0.31 / 0.161 / 0.099 for fly sorting). Context++ uses learned site-type intercepts per type (Agarwal et al., 2015 Figure 4).
- [x] **Context++ score** — For each site, CS = sum of 14 feature contributions from a **site-type-specific multiple linear regression** trained on 74 sRNA transfection datasets. Features: site type, supplementary pairing, local AU, minimum distance, sRNA1A/C/G, sRNA8A/C/G, site8A/C/G, 3′ UTR length, SA (site accessibility), ORF length, ORF 8mer count, 3′ UTR offset 6mer count, TA, SPS, PCT. Predicted site repression = `1 − 2^CS`. More negative CS → stronger repression. (Agarwal et al., 2015; [context_score.html](https://www.targetscan.org/vert_71/docs/context_score.html))
- [x] **Weighted Context++ score (CWCS)** — Iterative marginal repression accounting for 3P-seq AIR (isoform abundance) and downstream-site depletion: `Cᵢ = C(i−1) + (1 − 2^CSᵢ)(AIRᵢ − C(i−1))`, with C₀ = 0; `CWCS = log₂(1 − Cₙ)` where Cₙ is cumulative repression at the most proximal site. **More negative CWCS** = greater predicted repression (log<sub>2</sub> fold-change scale; matches sign behavior of Context++ in 8.0 bulk exports). For rank-by-percentile, use site-level **weighted context++ score percentile** (higher = stronger). (Agarwal et al., 2015; TargetScan 8.0 data download page)
- [x] **Total context++ score** — After removing overlapping sites (keeping set that **maximizes** total context score), sum Context++ scores for non-overlapping sites to the representative miRNA. All overlapping sites still shown on gene detail pages. (context_score_totals.html; Agarwal et al., 2015)
- [x] **PCT** — Probability that site conservation exceeds background **because of miRNA targeting**, estimated relative to 50 control k-mers matched for length and dinucleotide conservation rates, using phylogenetic branch lengths binned by UTR conservation (Friedman et al., 2009; updated in Agarwal et al., 2015). **Aggregate PCT** rolls up to gene–miRNA for conservation-based ranking. Contrast: PCT measures evolutionary evidence; Context++ measures predicted biochemical efficacy. PCT as a Context++ feature adds conservation signal to efficacy prediction but does not replace standalone PCT ranking.
- [x] **Evolutionary conservation model** — Conserved branch length (BLS) thresholds vary by site type and release (e.g., TargetScanHuman 7–8: 8mer ≥ 1.8, 7mer-m8 ≥ 2.6, 7mer-A1 ≥ 3.6). miRNA family conservation: broadly conserved (to zebrafish), conserved (placental mammals), poorly conserved. Site and family conservation are **independent** dimensions. (FAQs Release 8)
- [x] **Ranking — web vs download** — **Web (8.0):** default sort by biochemical predicted occupancy; options for CWCS and aggregate PCT. **Downloads:** `Summary_Counts` files provide Total context++ score, Cumulative weighted context++ score, Aggregate PCT, and Predicted occupancy (low/high/transfected miRNA). Site-level files provide context++ score, percentiles, weighted context++ score, predicted relative *K*<sub>D</sub>. BED files use context++ score percentile as score. (TargetScan 8.0 data download page)
- [x] **Version changes 7.x → 8.0** — (1) Default mammalian ranking → McGeary biochemical model + CNN. (2) Centered sites removed. (3) Same 3′ UTR profiles and miRNA families as Release 7. (4) Non-conserved target tables more accessible. (5) Context++/PCT ranking preserved as options. (TargetScan 8.0 help page; vert_80 home page)
- [x] **Score interpretation guide** — For highly conserved miRNA families, choose aggregate PCT (biological consequence / conserved targeting) **or** CWCS/context scores (biochemical efficacy, works for non-conserved miRNAs too). No single best cutoff. Consider co-expression. (FAQs Q11–12)
- [x] **McGeary 2019 biochemical model** — AGO-RBNS measured relative *K*<sub>D</sub> for six miRNAs across ≤12-nt sites. Biochemical model predicts AGO–miRNA occupancy on each mRNA from site *K*<sub>D</sub> values and miRNA concentration, then predicts repression. CNN extends *K*<sub>D</sub> prediction to all miRNAs from sequence. TargetScan 8.0 uses this as **default ranking** via predicted occupancy columns. Outperformed retrained Context++ in held-out miRNA cross-validation (McGeary et al., 2019).

### Site-level → gene–miRNA aggregation

1. **Filter sites** — For each miRNA family, retain 8mer, 7mer-m8, 7mer-A1, 6mer sites; remove overlaps (distal-to-proximal numbering for CWCS).
2. **Representative miRNA** — Compute CWCS for each family member; choose member with greatest predicted repression (ties: lowest miRNA name/number).
3. **Representative transcript** — UTR profile with most 3P-seq tags represents the gene.
4. **Outputs per gene–miRNA family:**
   - **Total context++ score** — sum of Context++ for non-overlapping sites to representative miRNA.
   - **CWCS** — AIR-weighted cumulative repression (formula above).
   - **Aggregate PCT** — conservation rollup for ranking (Friedman et al., 2009 method).
   - **Predicted occupancy** — biochemical model fraction of AGO–miRNA on target sites (8.0 default rank).

Poorly conserved sites receive the **same weight as conserved sites** in CWCS calculation (context_score_totals.html).

### R module: family-representative → all mature miRNAs (Option A)

**Interpretation:** TargetScan bulk `Summary_Counts` rows are **family-level** scores keyed by seed (`miRNA family` = `Seed+m8` in vertebrate 8.0 files) with one **Representative miRNA** label. Cross-tool integration matches on mature miRNA names (DIANA/miRDB/mirDIP), so exporting only the representative silently drops siblings (e.g. `hsa-miR-26a-5p` when the representative is `hsa-miR-26b-5p`).

**Implementation (`R/module_targetscan.R`):** After filtering to the query species (`Species ID` = NCBI taxon from `SPECIES$*_ncbi_tax_id`) and representative transcript, join `miR_Family_Info.txt` on seed (+ species) or legacy family name, emit one row per mature `MiRBase ID`, replicate CWCS/site counts, dedupe by strongest score per miRNA. Extra columns: `mirna_family`, `targetscan_representative_mirna`.

---

## Sources

- Lewis BP, Burge CB, Bartel DP (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. *Cell* 120:15–20.
- Grimson A, Farh KK-H, Johnston WK, Garrett-Engele P, Lim LP, Bartel DP (2007). MicroRNA targeting specificity in mammals: determinants beyond seed pairing. *Mol Cell* 27:91–105.
- Friedman RC, Farh KK-H, Burge CB, Bartel DP (2009). Most mammalian mRNAs are conserved targets of microRNAs. *Genome Res* 19:92–105.
- García DM, Baek D, Shin C, Bell GW, Grimson A, Bartel DP (2011). Weak seed-pairing stability and high target-site abundance decrease the proficiency of lsy-6 and other miRNAs. *Nat Struct Mol Biol* 18:1139–1146.
- Agarwal V, Bell GW, Nam J, Bartel DP (2015). Predicting effective microRNA target sites in mammalian mRNAs. *eLife* 4:e05005. doi:10.7554/eLife.05005
- Nam J-W, Rissland OS, Koppstein D, Abreu-Goodger C, Jan CH, Agarwal V, Lai EC, Bartel DP (2014). Global analyses of the effect of different cellular contexts on microRNA targeting. *Mol Cell* 53:1031–1043.
- McGeary SE, Lin KS, Shi CY, Pham T, Bisaria N, Kelley GM, Bartel DP (2019). The biochemical basis of microRNA targeting efficacy. *Science* 366(6472):eaav1741.
- TargetScan FAQs Release 7+: https://www.targetscan.org/faqs.Release_7.html
- TargetScan FAQs Release 8: https://www.targetscan.org/faqs.Release_8.html
- TargetScan context++ documentation: https://www.targetscan.org/vert_71/docs/context_score.html , https://www.targetscan.org/vert_71/docs/context_score_totals.html
- TargetScanHuman 8.0: https://www.targetscan.org/vert_80/
- TargetScan 8.0 data download (field descriptions): https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
- TargetScan 8.0 changelog: https://www.targetscan.org/vert_80/docs/help.html

---

## Integration notes (for consensus project)

**Recommended R-tool export (this project):**

| Column | Source field | Role |
|--------|--------------|------|
| `cwcs` | `Cumulative weighted context++ score` | Primary efficacy rank; **more negative = stronger** (`score_direction = lower_is_stronger` in R module) |
| `aggregate_pct` | `Aggregate PCT` | Conservation annotation overlay — **not averaged with efficacy** |
| `biochem_occupancy` | `Predicted occupancy (transfected miRNA)` | Optional 8.0-default efficacy rank; document which concentration scenario was chosen |

**TargetScan 8.0 bulk files to use** ([`scripts/download_targetscan.R`](../../scripts/download_targetscan.R)):

- **Gene–miRNA pairs:** `Summary_Counts.all_predictions.txt.zip` (full) or `Summary_Counts.default_predictions.txt.zip` (conserved sites of conserved families only)
- **Family membership:** `miR_Family_Info.txt.zip` — maps each **seed family** (Summary Counts `miRNA family` column = `Seed+m8` in 8.0 vertebrate downloads) to all mature miRNAs per NCBI taxon. The R module joins this file so every family member inherits the family-level CWCS (not only the **Representative miRNA** row).
- **Site-level Context++:** `Conserved_Site_Context_Scores.txt.zip` + `Nonconserved_Site_Context_Scores.txt.zip`
- **Site positions + PCT:** `Conserved_Family_Info.txt.zip` + `Nonconserved_Family_Info.txt.zip`

**Score selection guidance:**

- Default to **CWCS** for cross-tool comparison with DIANA/miRDB unless the analysis explicitly requires TargetScan 8.0 web-default behavior (then use predicted occupancy).
- Attach **aggregate PCT** as a separate conservation tier, not merged into efficacy.
- **Do not** naively average Context++, PCT, CWCS, biochemical occupancy, or TargetScan scores with DIANA/miRDB/mirDIP scores — different evidence types and scales.
- When mirDIP is used, remember it may already include TargetScan-derived evidence (correlated, not independent).

**PTEN example (interpretation):** The same PTEN–miRNA pair may rank high by CWCS (strong 8mer sites, favorable context) but lower by aggregate PCT (sites not preferentially conserved), explaining disagreement with conservation-heavy tools or with PCT-sorted TargetScan views.
