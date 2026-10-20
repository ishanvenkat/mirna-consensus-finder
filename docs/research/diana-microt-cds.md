# DIANA-microT-CDS

**Status:** complete

Use [template-algorithm-profile.md](template-algorithm-profile.md). Deep-dive checklist: [deep-dive-checklists.md](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md#diana-microt-cds).

Research plan: [RESEARCH_PLAN.md](RESEARCH_PLAN.md#2-diana-microt-cds)

---

## 1. What is it?

### Who developed it?

**DIANA-microT-CDS** was developed by the **DIANA Lab** (Artemis G. Hatzigeorgiou group; Fleming Institute / University of Thessaly ecosystem). Key authors include Martin Reczko, Manolis Maragkakis, Panagiotis Alexiou, Maria D. Paraskevopoulou, and colleagues.

### When was it introduced?

- **2012** — algorithm publication: Reczko et al., *Bioinformatics* (microT-CDS model for 3′ UTR + CDS targeting)
- **2013** — integrated into **DIANA-microT web server v5.0** (Paraskevopoulou et al., *Nucleic Acids Res.*)

### Current / relevant version

**DIANA-microT-CDS v5** via the [DIANA Tools](https://dianalab.e-ce.uth.gr/html/dianauniverse/) web server and REST API (`microT_CDSApi`). Supports *Homo sapiens*, *Mus musculus*, *Drosophila melanogaster*, and *C. elegans* (Paraskevopoulou et al., 2013; DIANA API docs).

### What problem does it solve?

DIANA-microT-CDS predicts **functional miRNA–gene interactions** by scanning both **3′ untranslated regions (3′ UTRs)** and **coding sequences (CDS)** — the first widely available online tool explicitly designed for dual-region targeting (Paraskevopoulou et al., 2013). It addresses cases where canonical 3′ UTR seed matches are absent but CDS sites or combined CDS+UTR targeting drive repression (Reczko et al., 2012).

---

## 2. How does it work?

```
miRNA + transcript (3′ UTR + CDS)
→ scan for candidate MREs (extended seed alignment, 6–9 nt)
→ score each MRE with region-specific GLM (CDS model vs 3′ UTR model)
→ sum MRE scores → region scores (CDS score, 3′ UTR score)
→ combine region scores with final GLM → gene–miRNA microT score
→ filter by user threshold (default 0.7 human)
```

### Step-by-step process

1. **Identify candidate miRNA recognition elements (MREs)** — dynamic programming aligns miRNA positions 1–9 (extended seed) to 9-nt windows in 3′ UTR or CDS; allows G:U wobbles, bulges, mismatches compensated by 3′ pairing (Reczko et al., 2012).
2. **Extract features per MRE** — binding category weight, conservation, free energy (RNAhybrid), accessibility (Sfold, 3′ UTR only), flanking AU content, inter-site distance, position within region, and binding-pattern descriptors (Reczko et al., 2012; Alexiou et al., 2014 review).
3. **Score each MRE** — separate **generalized linear models (GLMs)** for CDS and 3′ UTR, trained on PAR-CLIP true/false MRE sets (Reczko et al., 2012).
4. **Aggregate within region** — **sum** all MRE scores in CDS → CDS region score; sum all MRE scores in 3′ UTR → 3′ UTR region score (Reczko et al., 2012).
5. **Combine regions** — a second GLM trained on 13 microarray datasets (miRNA transfection/knockout expression changes) merges CDS and 3′ UTR region scores into the final **gene–miRNA microT score** (Reczko et al., 2012).
6. **Rank and filter** — results ordered by microT score; user applies threshold (API: 0.3–1).

### Candidate site identification

- Scans **both CDS and 3′ UTR** on the longest annotated transcript (Ensembl; Reczko et al., 2012).
- Supports **canonical and non-canonical** sites: 6mer–9mer extended-seed matches with G:U wobble, single bulge/mismatch (when ≥8 WC pairs), and additional 3′ compensatory pairing (Reczko et al., 2012).
- Positive training MREs derived from **PAR-CLIP** (Hafner et al., 2010); negative set = aligned locations not overlapping CLIP peaks (Reczko et al., 2012).

### Scoring and ranking

- Primary output: **microT score** (gene–miRNA level), higher = stronger predicted interaction (DIANA API docs; Reczko et al., 2012 precision–recall curves use score cutoffs).
- Web results also report **binding type**, **conservation**, and cross-references to TarBase / TargetScan / miRanda overlaps (Alexiou et al., 2014).

---

## 3. What biological/computational features does it use?

| Feature | Used? | Source / notes |
|---------|-------|----------------|
| Seed matching | ✅ | Extended seed (nt 1–9); canonical + non-canonical rules — Reczko et al., 2012 |
| 8mer / 7mer-m8 / 7mer-A1 / 6mer site types | Partial | Uses **64 binding categories** derived from extended-seed alignment patterns, not TargetScan site-type taxonomy — Reczko et al., 2012 |
| Evolutionary conservation | ✅ | Separate models: CDS (30-species excess conservation above amino-acid constraint) vs 3′ UTR (16-species ratio) — Reczko et al., 2012 |
| Binding / free energy | ✅ | RNAhybrid minimum free energy — Reczko et al., 2012 |
| RNA / site accessibility | ✅ (3′ UTR) | Sfold accessibility sum at positions −1, 1, 2 relative to MRE start — Reczko et al., 2012; **not documented for CDS** |
| Local AU content | ✅ | MRE flanking AU content (Grimson et al., 2007 feature) — Reczko et al., 2012 |
| Target-site abundance | Partial | Adjacent MRE distance models multi-site spacing — Reczko et al., 2012 |
| Position within transcript | ✅ | Distance to closest 3′ UTR end; region-specific models — Reczko et al., 2012 |
| 3′ compensatory pairing | ✅ | Non-canonical sites supported by additional 3′ pairing — Reczko et al., 2012 |
| G:U pairing | ✅ | Allowed in extended-seed alignment rules — Reczko et al., 2012 |
| Machine learning | ✅ | **GLMs** (not deep learning): per-MRE GLMs + region-combination GLM — Reczko et al., 2012 |
| Other sequence/context features | ✅ | Binding-pattern descriptors (bnt1–bnt29), binding category weight from logistic regression on PAR-CLIP — Reczko et al., 2012 |

---

## 4. What does its score mean?

| Question | Answer |
|----------|--------|
| Score name(s) | **microT score** (gene–miRNA); **MRE score** (site); **CDS score** / **3′ UTR score** (region sums) |
| Numerical range | **0.3–1** used as operational cutoff range (DIANA API); full GLM output range not explicitly bounded in primary paper — treat as unbounded GLM output filtered by threshold |
| Higher or lower = stronger? | **Higher = stronger** predicted interaction |
| Score type | GLM-derived interaction strength / probability-like ranking metric |
| Absolute interpretation across genes? | Partial — score reflects predicted repression potential for a specific miRNA–gene pair; cross-gene comparison is ranking-oriented |
| Already normalized? | Features normalized within model; final score is model output, not percentile-normalized across genome |
| Recommended cutoff / threshold | Human: **0.7** (sensitive), **0.8** (default), **0.9** (stringent); Mouse: **0.6/0.7/0.8** — DIANA REST API docs |
| Site scores aggregated to gene–miRNA score? | **Yes** — MRE scores summed per region; CDS + 3′ UTR region scores combined by final GLM (Reczko et al., 2012) |

### Score interpretation notes

- The microT score integrates **binding evidence** (PAR-CLIP-trained MRE features) with **functional repression evidence** (microarray-trained region combination). It is an **efficacy-oriented** score, not a standalone conservation metric.
- Conservation enters as a **feature within the GLM**, unlike TargetScan where PCT is a separate conservation axis.
- For integration planning: microT scores are **not comparable** to TargetScan Context++/CWCS, miRDB MirTarget scores, or mirDIP integrated scores without rank fusion or percentile normalization.

### 4b. Score inventory

| Score name | Level | Evidence type | Range | Direction | Notes |
|------------|-------|---------------|-------|-----------|-------|
| MRE score | site (MRE) | efficacy | Not documented | Higher = stronger | Per-site GLM output; summed within region — Reczko et al., 2012 |
| CDS region score | region (CDS) | efficacy | Not documented | Higher = stronger | Sum of CDS MRE scores — Reczko et al., 2012 |
| 3′ UTR region score | region (3′ UTR) | efficacy | Not documented | Higher = stronger | Sum of 3′ UTR MRE scores — Reczko et al., 2012 |
| microT score | gene–miRNA | efficacy | Cutoff 0.3–1 documented | Higher = stronger | Final GLM combining region scores — Reczko et al., 2012; DIANA API |

---

## 5. Advantages

- **Dual-region scanning (CDS + 3′ UTR)** captures targets missed by 3′ UTR–only tools; Reczko et al. (2012) report sensitivity gains on Selbach proteomics data (52% → 65% vs 3′ UTR–only model).
- **Non-canonical site support** (wobble, bulge, 3′ pairing) expands candidate pool beyond strict canonical seeds.
- **Region-specific models** respect different biology of CDS vs 3′ UTR binding (distinct feature sets and conservation logic).
- **Higher sensitivity** than TargetScan 5.0, PicTar, RNA22, miRanda, and DIANA-microT v3.0 on independent proteomics benchmark at matched specificity (Reczko et al., 2012).
- **REST API and bulk downloads** support programmatic pipelines (DIANA API docs; Paraskevopoulou et al., 2013).

---

## 6. Disadvantages

- **Training data age** — core GLM trained on PAR-CLIP (2010) and microarray sets from late 2000s; model may not reflect current transcript annotations or miRNA catalogs.
- **GLM complexity** — less transparent than rule-based site-type models (TargetScan); harder to decompose score into interpretable components.
- **CDS accessibility not modeled** — Sfold accessibility applied to 3′ UTR only (Reczko et al., 2012).
- **Score scale undocumented beyond cutoff range** — full numerical range and cross-tool comparability not specified in primary sources.
- **Species limited** — four species on web server vs miRDB's five or TargetScan's broad mammalian coverage.

---

## 7. When would you use it?

Use DIANA-microT-CDS when you suspect **CDS-mediated targeting** or when 3′ UTR–only tools (e.g., TargetScan) return no hits for a biologically plausible miRNA–gene pair — e.g., genes with **short 3′ UTRs** where CDS sites are enriched (Reczko et al., 2012). Also useful when you want **non-canonical site** sensitivity and a single efficacy score combining both regions. For PTEN disagreement analysis: DIANA may report interactions driven by CDS sites or non-canonical binding that TargetScan's canonical 3′ UTR seed taxonomy misses.

---

## 8. Method-specific deep dive

- [x] **What microT-CDS stands for and what it predicts** — **microT** (DIANA microRNA target algorithm) extended to **CDS**; predicts functional miRNA–gene interactions in 3′ UTR and coding sequences (Reczko et al., 2012; Paraskevopoulou et al., 2013).
- [x] **How candidate sites are identified** — Extended-seed (nt 1–9) dynamic programming scan in both CDS and 3′ UTR with non-canonical alignment rules (Reczko et al., 2012).
- [x] **How seed matching is used** — 64 binding categories from extended-seed alignment patterns; binding category weight feature from logistic regression on PAR-CLIP enrichment (Reczko et al., 2012).
- [x] **How conservation is used** — CDS: excess nucleotide conservation above amino-acid constraint (30 species); 3′ UTR: conserved seed positions / max conservation ratio (16 species) — Reczko et al., 2012.
- [x] **How free energy is used** — RNAhybrid binding free energy as GLM feature (Reczko et al., 2012).
- [x] **How site accessibility is modeled** — Sfold-derived accessibility at positions −1, 1, 2 for **3′ UTR MREs only** (Reczko et al., 2012).
- [x] **How CDS targeting differs from 3′ UTR targeting** — Separate feature sets, separate GLMs, separate conservation models; combined by final GLM trained on expression data (Reczko et al., 2012).
- [x] **What machine-learning component is used** — **Generalized linear models (GLMs)** at MRE, region-combination, and binding-category levels; compared against SVM/neural nets/random forests — GLMs selected by cross-validation (Reczko et al., 2012).
- [x] **What its prediction score represents** — Predicted functional miRNA–gene repression potential integrating binding and expression evidence (Reczko et al., 2012).
- [x] **Its score range and cutoff** — User threshold **0.3–1**; human defaults **0.7/0.8/0.9** (DIANA REST API); full score range not explicitly documented.
- [x] **Advantages/disadvantages compared with TargetScan** — See subsection below.

### Comparison with TargetScan

| Aspect | DIANA-microT-CDS | TargetScan |
|--------|------------------|------------|
| Regions scanned | 3′ UTR + **CDS** | Primarily 3′ UTR (representative transcript); ORF sites in downloads |
| Site definition | 64 binding categories; non-canonical + 3′ pairing | Canonical 8mer/7mer/6mer + 3′ compensatory taxonomy |
| Model type | GLM (PAR-CLIP + microarray trained) | Context++ linear model + separate PCT conservation |
| Conservation | Feature inside efficacy GLM | **Separate PCT axis** independent of Context++ |
| Multi-site aggregation | Sum MRE → region; GLM merge regions | Sum Context++ / CWCS / aggregate PCT (distinct formulas) |
| Sensitivity (Selbach 2008) | Highest at any specificity vs TargetScan 5.0 and others (Reczko 2012) | Strong on canonical conserved 3′ UTR sites; misses CDS-only targets |
| Interpretability | Lower — composite GLM score | Higher — site types, Context++ features, PCT separable |

**Why they disagree for PTEN:** TargetScan ranks by biochemical occupancy / Context++ on canonical 3′ UTR sites; DIANA may additionally score **CDS MREs** and **non-canonical** alignments with different conservation and accessibility weighting — yielding different gene–miRNA pairs and rankings.

### Site-level → gene–miRNA aggregation

1. **MRE → region:** sum all MRE GLM scores within CDS; separately sum within 3′ UTR (Reczko et al., 2012).
2. **Region → gene–miRNA:** final GLM combines CDS score + 3′ UTR score, trained on microarray top/bottom 100 genes per miRNA perturbation experiment (Reczko et al., 2012).

---

## Sources

- Reczko M, Maragkakis M, Alexiou P, Grosse I, Hatzigeorgiou AG (2012). Functional microRNA targets in protein coding sequences. *Bioinformatics* 28(6):771–776. doi:10.1093/bioinformatics/bts043
- Paraskevopoulou MD, Georgiakakis G, Koskinas T, Reczko M, Maragkakis M, et al. (2013). DIANA-microT web server v5.0: service integration into miRNA functional analysis workflows. *Nucleic Acids Res* 41(W1):W169–W173. doi:10.1093/nar/gkt393
- Alexiou P, Maragkakis M, Papadopoulos GL, Reczko M, Hatzigeorgiou AG (2014). Common features of microRNA target prediction tools. *Front Genet* 5:23. doi:10.3389/fgene.2014.00023
- DIANA Tools REST API help — microT-CDS threshold and endpoints: http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=site%2Fhelp&topic=api
- Hafner M, Landthaler M, Burger L, Khorshid M, et al. (2010). PAR-CLIP training data. *Cell* 141(1):129–141.

---

## Integration notes (for consensus project)

- Export **microT score** as DIANA's primary **efficacy** rank; do not average with TargetScan CWCS or miRDB MirTarget score.
- DIANA may add **CDS-only** evidence absent from TargetScan — treat as complementary tier, not duplicate confirmation.
- Recommended R module threshold: **0.7** (sensitive) or **0.8** (default) per DIANA API; expose as user parameter in `R/module_diana.R`.
