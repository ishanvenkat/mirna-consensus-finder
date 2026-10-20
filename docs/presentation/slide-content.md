# Slide Content — miRNA Algorithm Comparison



Paste-ready bullets for the presentation deck. Max 5 bullets per slide.



**Visual assets:** [`assets/`](assets/) — SVG diagrams referenced per slide below.



---



## Slide 1: Why miRNA Target Prediction Is Hard



- miRNA–mRNA pairing is **imperfect** — functional repression often requires only seed match (positions 2–7/8), not full complementarity (Lewis et al., 2005)

- Combinatorial scale: thousands of miRNAs × tens of thousands of genes → millions of candidate interactions

- Exhaustive experimental testing of every pair is impractical — reporter assays, CLIP, and perturbation studies are costly and incomplete

- Algorithms must **prioritize** from partial sequence and context evidence (conservation, accessibility, site type, etc.)

- Different evidence choices produce different ranked lists — the core question this deck addresses



---



## Slide 2: miRNAs and Targets



- **miRNAs:** ~22-nt non-coding RNAs that repress translation and/or trigger mRNA degradation

- **Seed region** (nt 2–7/8 from the 5′ end) is the primary determinant of target recognition (Lewis et al., 2005)

- Target sites are most often in the **3′ UTR**; DIANA-microT-CDS also scans **CDS** (Reczko et al., 2012)

- Prediction unit: **miRNA → target site(s) on transcript → regulated gene**

- One miRNA can target many genes; one gene can harbor sites for many miRNAs



---



## Slide 3: Prediction vs Experimental Validation



**Visual:** [`assets/slide03-prediction-vs-validation.svg`](assets/slide03-prediction-vs-validation.svg)



- **Predictors (computational):** TargetScan, DIANA-microT-CDS, miRDB, mirDIP — ask *could* this miRNA regulate this gene?

- **Validation databases (experimental):** miRTarBase, DIANA-TarBase — ask *has* this interaction been supported in the lab?

- Validation DBs are **not prediction algorithms** — never treat them as predictors or mix their evidence into predictor scores

- Evidence types in validation DBs: reporter assays, qPCR/Western blot, CLIP, knockdown/overexpression

- Our pipeline: rank predictions first → optionally overlay experimental support as a read-only reference layer



---



## Slide 4: How Predictors Generally Work



**Visual:** [`assets/slide04-predictor-pipeline.svg`](assets/slide04-predictor-pipeline.svg)



- **Input:** miRNA sequence + mRNA/transcript sequence (3′ UTR; DIANA also CDS)

- **Candidate discovery:** scan for seed-matched miRNA recognition elements (MREs) — rules differ by tool

- **Site scoring:** apply tool-specific features — thermodynamics, conservation, accessibility, ML, etc.

- **Aggregation:** combine site scores → single gene–miRNA score; apply tool-specific threshold

- **Output:** ranked miRNA list per gene — score **semantics differ** (efficacy, conservation, ML confidence, meta-integration)



---



## Slide 5: TargetScan Overview



- Bartel-lab seed-based predictor: scans 3′ UTRs for canonical 8mer, 7mer, and 6mer sites (Lewis et al., 2005)

- miRNAs grouped into families by shared seed; targets ranked per gene–miRNA family pair

- Uses representative transcript (most 3P-seq tags) with isoform-aware 3′ UTR profiles (Nam et al., 2014)

- Reports conserved and non-conserved sites; conservation and efficacy scored independently

- Current release: TargetScan 8.0 — same site inventory as 7.x, new default ranking method



---



## Slide 6: TargetScan Scoring System



**Visual:** [`assets/slide06-targetscan-scoring-axes.svg`](assets/slide06-targetscan-scoring-axes.svg)



- **Context++** (site): sum of 14 features predicting site-level repression — more negative = stronger (Agarwal et al., 2015)

- **CWCS** (gene–miRNA): cumulative AIR-weighted repression across non-overlapping sites — **more negative = stronger** (log₂ fold-change; same convention as Context++)

- **PCT** (site / aggregate): probability of conserved targeting — separate conservation axis, not efficacy (Friedman et al., 2009)

- **TargetScan 8.0 default:** biochemical predicted occupancy from CNN + *K*<sub>D</sub> model (McGeary et al., 2019)

- Do not average Context++, CWCS, PCT, or occupancy — different scales, directions, and evidence types



---



## Slide 7: DIANA-microT-CDS Overview



- DIANA Lab algorithm scanning **both 3′ UTR and CDS** — only major online tool designed for dual-region targeting (Reczko et al., 2012)

- Identifies canonical and **non-canonical** MREs (G:U wobble, bulge, 3′ compensatory pairing) via extended-seed alignment

- Separate **GLM models** for CDS vs 3′ UTR, trained on PAR-CLIP binding sites (Reczko et al., 2012)

- Region scores merged by a second GLM trained on 13 microarray miRNA perturbation datasets

- Higher sensitivity than TargetScan 5.0 on independent proteomics benchmark at matched specificity



---



## Slide 8: DIANA Scoring System



- **MRE score** (site): GLM output per recognition element; summed within each region (Reczko et al., 2012)

- **microT score** (gene–miRNA): final GLM combining CDS + 3′ UTR region scores — **higher = stronger**

- Features: binding category, conservation (30 spp CDS / 16 spp UTR), RNAhybrid ΔG, Sfold accessibility (3′ UTR), flanking AU

- User threshold **0.3–1**; human defaults **0.7** (sensitive) / **0.8** (default) / **0.9** (stringent) — DIANA API

- Efficacy score only — conservation is a GLM feature, not a separate axis (contrast TargetScan PCT)



---



## Slide 9: miRDB



- **MirTarget v4.0**: SVM (RBF kernel) with **96 features** trained on CLIP ligation + RNA-seq downregulation data (Liu & Wang, 2019)

- Scans **3′ UTR only** for seed8A1, seed7b, seed7A1 — conservation **not required** for candidacy (Chen & Wang, 2020)

- **MirTarget score** (50–100): S = 100 × (1 − ∏ site probabilities) — higher = stronger ML confidence

- Practical cutoffs: **> 80** high confidence; **< 60** cautious — miRDB FAQ

- ML contrast to rule-based TargetScan/DIANA: black-box but captures non-conserved functional targets



---



## Slide 10: mirDIP



- **Integrative meta-resource**, not a single algorithm — aggregates **24 prediction databases** (Hauschild et al., 2023)

- Includes TargetScan v7.2, DIANA, miRDB v6, PITA, RNAhybrid, miRanda, and 18 others — mirDIP statistics

- **Integrated score**: noisy-or of precision-calibrated confidence scores — proportional to interaction probability (Tokar et al., 2018)

- **Score classes** (per miRNA): Very High (top 1%), High (top 5%), Medium (top 1/3), Low (bottom 2/3)

- Do not equal-weight mirDIP with its constituent tools — correlated evidence / double-counting risk



---



## Slide 11: Validation Databases



- **miRTarBase 2025** (~3.8M MTIs): manually curated experimental interactions — **not a predictor** (Cui et al., 2025)

- **Strong evidence:** reporter assay, Western blot, qPCR; **weaker:** CLIP-seq, PAR-CLIP, microarray (Cui et al., 2025)

- **DIANA-TarBase v9.0** (~6M entries, ~2M unique pairs): largest validated collection with cell-type binding sites (Skoufos et al., 2024)

- TarBase: 37 protocols (low-yield + high-throughput); filters for tissue, cell type, direct/indirect method, 3′ UTR vs CDS

- Project use: read-only **`experimental_support`** overlay via multiMiR — never merge into predictor consensus scores



---



## Slide 12: Side-by-Side Comparison



**Visual:** [`assets/slide12-evidence-types.svg`](assets/slide12-evidence-types.svg)



- **TargetScan:** CWCS / occupancy (efficacy) + PCT (conservation, separate axis) — 3′ UTR canonical sites; no universal cutoff

- **DIANA-microT-CDS:** microT GLM score (efficacy) — **3′ UTR + CDS**, non-canonical sites; threshold 0.7–0.9

- **miRDB:** MirTarget 50–100 (ML confidence) — 3′ UTR only; >80 high confidence, <60 cautious

- **mirDIP:** Integrated 0–1 (meta-integration) — aggregates 24 predictors; Score Class per miRNA — not independent

- **Validation overlay:** miRTarBase/TarBase provide **experimental support** tiers — categorical, read-only; never averaged with predictor scores



---



## Slide 13: Why Algorithms Disagree (PTEN Example)



**Visual:** [`assets/slide13-pten-disagreement.svg`](assets/slide13-pten-disagreement.svg)



- **PTEN** is a well-studied tumor suppressor — dozens of miRNAs reported as regulators; tools return **different top-ranked lists** for the same gene

- **miR-21-5p → PTEN:** luciferase reporter validated; ≥52 literature references; TargetScan-predicted 3′ UTR sites (Meng et al. 2007; Chen et al. 2014)

- **miR-214-3p → PTEN:** reporter validated (Yang et al. 2008); present in TarBase / miRTarBase — may rank lower in conservation-heavy predictors

- Methodological drivers: **CDS vs 3′ UTR scanning**, canonical vs non-canonical seeds, conservation vs ML weighting, mirDIP double-counting constituent tools

- Disagreement reflects **different evidence and assumptions**, not necessarily that one tool is simply wrong



---



## Slide 14: Advantages and Disadvantages

- **TargetScan:** interpretable site types + separate PCT/conservation axis; limited by 3′ UTR focus and multiple incompatible scores (CWCS vs occupancy vs PCT)
- **DIANA-microT-CDS:** dual-region (CDS + UTR) + non-canonical sites + high proteomics sensitivity; limited by opaque GLM and accessibility modeled in UTR only
- **miRDB:** integrative CLIP + RNA-seq ML training + non-conserved targets; limited by black-box SVM, 3′ UTR-only, three seed types
- **mirDIP:** widest coverage + calibrated integrative score outperforming individual tools; limited by correlated constituent evidence — not an independent fourth vote
- **Validation DBs (miRTarBase/TarBase):** ground truth for tested pairs with strong/weak tiers; limited by sparse coverage — absence ≠ disproven interaction



---



## Slide 15: Implications for Consensus Scoring



**Visual:** [`assets/slide15-consensus-scoring-flow.svg`](assets/slide15-consensus-scoring-flow.svg)



- **Cannot naively average** scores — mix efficacy, conservation, ML confidence, and meta-integration on incompatible scales (see score matrix)

- **Per-source percentile normalization** or **rank fusion** (e.g., Borda count) before combining across tools

- Export **one primary rank per predictor** (CWCS, microT, MirTarget); attach PCT and validation flags as **separate tiers**

- Treat **mirDIP as meta-layer** — do not add TargetScan + DIANA + miRDB + mirDIP as four equal votes

- **Validation overlay** (`experimental_support`, strong/weak) via multiMiR — annotates predictions, never replaces them



---



## Slide 16: References



Key sources — full bibliography in [`sources/bibliography.md`](sources/bibliography.md).



- Lewis BP et al. (2005). *Cell* — seed matching

- Agarwal V et al. (2015). *eLife* — Context++

- Friedman RC et al. (2009). *Genome Res* — PCT

- McGeary SE et al. (2019). *Science* — TargetScan 8.0 biochemical model

- Reczko M et al. (2012). *Bioinformatics* — DIANA-microT-CDS

- Liu W, Wang X (2019). *Genome Biology* — MirTarget / miRDB

- Tokar T et al. (2018). *NAR* — mirDIP integrative score

- Hauschild AC et al. (2023). *NAR* — mirDIP 5.2

- Cui S et al. (2025). *NAR* — miRTarBase 2025

- Skoufos G et al. (2024). *NAR* — TarBase v9.0

- Yang H et al. (2008). *Cancer Res* — miR-214 → PTEN (PTEN example)

- Meng F et al. (2007). *Exp Biol Med* — miR-21 → PTEN (PTEN example)



---



## Asset index



| Slide | Asset |

|-------|-------|

| 3 | `assets/slide03-prediction-vs-validation.svg` |

| 4 | `assets/slide04-predictor-pipeline.svg` |

| 6 | `assets/slide06-targetscan-scoring-axes.svg` |

| 12 | `assets/slide12-evidence-types.svg` |

| 13 | `assets/slide13-pten-disagreement.svg` |

| 15 | `assets/slide15-consensus-scoring-flow.svg` |



Speaker notes: [`speaker-notes.md`](speaker-notes.md)


