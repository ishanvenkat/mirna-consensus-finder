# mirDIP

**Status:** complete

**Important:** mirDIP is an **integrative database/resource**, not a single prediction algorithm like TargetScan.

Use [template-algorithm-profile.md](template-algorithm-profile.md). Deep-dive checklist: [deep-dive-checklists.md](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md#mirdip).

Research plan: [RESEARCH_PLAN.md](RESEARCH_PLAN.md#3-mirdip)

Conceptual flow: TargetScan + DIANA + miRDB + PITA + RNAhybrid + 22 others → mirDIP → integrated score.

---

## 1. What is it?

### Who developed it?

**mirDIP** (microRNA Data Integration Portal) is developed by the **Jurisica Lab** (Krembil Research Institute, University Health Network, Toronto). Key contributors include Igor Jurisica, Tomáš Tokár, Chiara Pastrello, Anne-Christin Hauschild, and colleagues.

### When was it introduced?

- **2011** — NAViGaTOR-based multi-database integration concept (Shirdel et al., *PLoS ONE*)
- **2018** — mirDIP v4.1 with integrative score (Tokar et al., *Nucleic Acids Res.*)
- **2023** — mirDIP v5.2 with tissue context and mirGeneDB support (Hauschild et al., *Nucleic Acids Res.*)

### Current / relevant version

- **Portal:** v5.3.0.3 (April 2026, per live site)
- **Database:** v5.2.2.1 — **46,364,047** integrated predictions for **27,936** genes and **2,734** miRNAs (Hauschild et al., 2023; mirDIP statistics page)
- URL: https://ophid.utoronto.ca/mirDIP

### What problem does it solve?

Individual miRNA target predictors disagree and each covers only part of the interaction space. mirDIP **aggregates human miRNA–gene predictions from 24 computational resources**, calibrates their confidence scores, and assigns a unified **integrative score** plus **confidence class** per interaction — increasing coverage while providing a single ranking metric (Tokar et al., 2018; Hauschild et al., 2023).

---

## 2. How does it work?

```
24 prediction resources (downloaded or re-run)
→ standardize gene symbols (HGNC) and miRNA names (miRBase v22)
→ rank-normalize each resource's scores to [0,1]
→ benchmark each resource → map ranks to precision-based confidence scores
→ noisy-or combine confidence scores per miRNA–gene pair
→ assign confidence class (Very High / High / Medium / Low) per miRNA
→ store + serve via web/API
```

### Step-by-step process

1. **Collect predictions** — Download or execute 24 resources; mirDIP 5.2 also **re-ran** five algorithms (miRanda, BiTargeting, PITA, RNAhybrid, MirMAP) on miRBase v22 and mirGeneDB 2.0 3′ UTRs (Hauschild et al., 2023).
2. **Standardize identifiers** — HGNC gene symbols; miRBase v22 miRNA names (Tokar et al., 2018; Hauschild et al., 2023).
3. **Rank-normalize within each resource** — Assign rank *r* ∈ [0,1] where 0 = most confident prediction in that resource (Tokar et al., 2018).
4. **Calibrate to confidence scores** — For each resource, fit rank→precision curve against experimental validation sets (TarBase, NPInter); interpolate precision for all predictions (Tokar et al., 2018).
5. **Integrate with noisy-or** — For each unique miRNA–gene pair: **S = 1 − ∏(1 − sᵢ)** across all resources predicting that pair (Tokar et al., 2018).
6. **Assign confidence class** — Per miRNA percentile bins (mirDIP 5.2): top 1% = Very High; top 5% = High; top 1/3 = Medium; bottom 2/3 = Low (Hauschild et al., 2023; mirDIP portal).

### Which prediction resources are integrated

**24 resources in mirDIP 5.2** (mirDIP statistics page; Hauschild et al., 2023):

| Resource | Notes |
|----------|-------|
| BiTargeting | Re-run May 2021 |
| BCmicrO | Downloaded |
| CoMeTa | Downloaded |
| Cupid | Downloaded |
| **DIANA** | microT-CDS predictions |
| miRanda | Re-run May 2021 |
| mirbase | Downloaded |
| mirCoX | Downloaded |
| MirMAP | Re-run May 2021 |
| Mirza-G | Downloaded |
| miRcode | Downloaded |
| **miRDB v6** | MirTarget predictions |
| miRTar2GO | Downloaded |
| MBStar | Downloaded |
| MirAncesTar | Downloaded |
| MirSNPInTarget | Downloaded |
| MirTar2 | Downloaded |
| **miRNATIP** | Added in 5.2 |
| MultiMiTar | Downloaded |
| PACCMIT | Downloaded |
| **PITA** | Re-run May 2021 |
| **RNAhybrid** | Re-run May 2021 (−22 kcal/mol filter) |
| RNA22 | Downloaded |
| **TargetScan v7.2** | Downloaded |

Eight outdated resources removed in 5.2 (not updated >10 years; Hauschild et al., 2023).

### How individual predictions are combined

- **Not** a simple vote count or arithmetic mean.
- Each resource's raw score is **rank-normalized**, then **calibrated to a precision estimate** (confidence score) using experimental validation benchmarks (Tokar et al., 2018).
- Confidence scores from all resources predicting the same pair are merged via the **noisy-or** formula: **Sⱼ = 1 − ∏ᵢ(1 − sᵢⱼ)** (Tokar et al., 2018).
- **Interpretation (documented):** integrative score is **proportional to the probability of the interaction occurring** (Hauschild et al., 2023; mirDIP portal).

### Scoring and ranking

- Results sorted by **Integrated Score** (descending).
- **Number of Sources** column shows how many of the 24 resources support the pair.
- Users filter by **Score Class** (Very High / High / Medium / Low) on the web portal.

---

## 3. What biological/computational features does it use?

mirDIP **does not perform its own site finding**. It inherits features indirectly from constituent predictors.

| Aspect | Documented? | Source / notes |
|--------|-------------|----------------|
| Aggregates TargetScan predictions | ✅ | TargetScan v7.2 in 24 resources — mirDIP statistics |
| Aggregates DIANA predictions | ✅ | DIANA-microT-CDS — mirDIP statistics |
| Aggregates miRDB predictions | ✅ | miRDB v6 — mirDIP statistics |
| Aggregates PITA / RNAhybrid / others | ✅ | PITA, RNAhybrid, RNA22, miRanda, etc. — mirDIP statistics |
| Correlated-evidence risk | ✅ | TargetScan, DIANA, and miRDB overlap partially; Jaccard overlap between individual resources is low on average (mean 0.06) but mirDIP **includes all three** — Tokar et al., 2018 |

---

## 4. What does its score mean?

| Question | Answer |
|----------|--------|
| Score name(s) | **Integrated score**; **Score class** (categorical) |
| Numerical range | **0–1** (probability-like); mean **0.05**, median **0.02** across predictions (Tokar et al., 2018) |
| Higher or lower = stronger? | **Higher = stronger** |
| Score type | Noisy-or of calibrated precision/confidence scores |
| Absolute interpretation across genes? | Partial — score approximates interaction probability but distribution is heavily right-skewed (most scores low) |
| Already normalized? | Yes — per-resource rank normalization + precision calibration before integration |
| Recommended cutoff / threshold | Use **Score Class** filters: Very High (top 1% per miRNA), High (top 5%), Medium (top 1/3), Low (bottom 2/3) — mirDIP portal; Hauschild et al., 2023 |
| Confidence categories | Very High / High / Medium / Low — assigned **per miRNA** in v5.2 (not global percentiles) |

### 4b. Score inventory

| Score name | Level | Evidence type | Range | Direction | Notes |
|------------|-------|---------------|-------|-----------|-------|
| Integrated score | gene–miRNA | meta-integration | 0–1 (typical median 0.02) | Higher = stronger | Noisy-or of resource confidence scores — Tokar et al., 2018 |
| Score class | gene–miRNA | meta-integration | Categorical (4 levels) | Very High > High > Medium > Low | Per-miRNA percentile bins since v5.2 — Hauschild et al., 2023 |
| Number of Sources | gene–miRNA | meta-integration | 1–24 (integer) | Higher = more resources agree | Count only; not a calibrated probability |

### Score interpretation notes

- mirDIP is **already meta-integrated** — using mirDIP alongside TargetScan, DIANA, and miRDB in equal-weight consensus **double-counts** the same underlying predictions.
- Integrated score measures **confidence that an interaction exists** across methods, not biochemical repression magnitude (contrast TargetScan Context++/CWCS).
- Low median (0.02) means absolute score values are small; **rank and Score Class** are more useful than raw thresholds.

---

## 5. Advantages

- **Broader coverage** — 46M+ integrated predictions vs any single tool (Hauschild et al., 2023).
- **Unified ranking** — single integrative score replaces comparing incompatible per-tool scales (Tokar et al., 2018).
- **Calibrated confidence** — precision benchmarking against TarBase/NPInter experimental evidence (Tokar et al., 2018).
- **Integrative score outperforms individual resources** on F-score vs experimental validation (Tokar et al., 2018, Figure 1).
- **No bias accumulation** — integrating across resources does not cumulate GO/pathway prediction bias (Tokar et al., 2018).
- **Tissue/disease context** in v5.2 — 330 contexts for filtering (Hauschild et al., 2023).

---

## 6. Disadvantages

- **Correlated evidence** — includes TargetScan, DIANA, and miRDB simultaneously; consensus with those tools risks **double-counting** (project rule).
- **Version heterogeneity** — integrated TargetScan v7.2 while current TargetScan is 8.0; resource versions differ in annotation and scoring (mirDIP statistics).
- **Not a primary algorithm** — cannot explain *why* a specific site is predicted; must trace to constituent resource.
- **Human-only** — unlike TargetScan/miRDB multi-species support.
- **Skewed score distribution** — most integrated scores near zero; raw numeric thresholds less intuitive than Score Class.
- **Binding-energy-only tools** (RNAhybrid, RNA22, BiTargeting) contribute lower-confidence predictions per mirDIP's own analysis (Tokar et al., 2018, Supplementary Figure 2).

---

## 7. When would you use it?

Use mirDIP for a **first-pass meta-ranking** of miRNA–gene interactions when you want maximum recall and a calibrated confidence tier — e.g., "which miRNAs target PTEN with Very High integrative support?" Do **not** use mirDIP as an independent equal-weight vote alongside TargetScan/DIANA/miRDB in naive consensus scoring; instead, treat mirDIP as a **meta-layer** that already subsumes those sources. Filter by **Score Class ≥ High** for hypothesis prioritization; use **Number of Sources** as a secondary support indicator.

---

## 8. Method-specific deep dive

- [x] **Which prediction resources mirDIP integrates (verify current version)** — **24 resources** listed above; portal DB v5.2.2.1 (mirDIP statistics; Hauschild et al., 2023).
- [x] **How individual predictions are combined** — Rank normalization → precision calibration → **noisy-or** (Tokar et al., 2018).
- [x] **How the mirDIP integrated score is calculated** — **S = 1 − ∏(1 − sᵢ)** where sᵢ are per-resource confidence (precision) scores (Tokar et al., 2018).
- [x] **What the integrated score means** — Proportional to **probability of interaction occurring** (Hauschild et al., 2023; mirDIP portal).
- [x] **Its score range and confidence categories** — Range 0–1; mean 0.05 / median 0.02 (Tokar et al., 2018); classes Very High (top 1%), High (top 5%), Medium (top 1/3), Low (bottom 2/3) **per miRNA** (Hauschild et al., 2023).
- [x] **How mirDIP evaluates predictions against experimentally supported interactions** — Benchmark against TarBase v7 + NPInter v3; F-score on top-*h* predictions; precision–rank curves per resource (Tokar et al., 2018); v5.2 uses miRTarBase v8 + NPInter v4 (Hauschild et al., 2023).
- [x] **Advantages of integrating multiple algorithms** — Higher F-score, broader coverage, reduced single-method bias (Tokar et al., 2018).
- [x] **Disadvantages of relying on an integrated resource** — Loss of algorithm-specific interpretability; version staleness; cannot decompose site-level mechanism.
- [x] **Whether combining methods creates correlated evidence or other biases** — Individual resource overlap is low (mean Jaccard 0.06), but mirDIP **includes overlapping tools** (TargetScan + DIANA + miRDB); integrative score bias is moderate vs individual high-bias tools (Tokar et al., 2018).

---

## Sources

- Tokar T, Pastrello C, Rossos AEM, Abovsky M, Hauschild AC, Tsay MJ, Lu R, Jurisica I (2018). mirDIP 4.1 — integrative database of human microRNA target predictions. *Nucleic Acids Res* 46(D1):D360–D370. doi:10.1093/nar/gkx1144
- Hauschild AC, Pastrello C, Ekaputeri GKA, Bethune-Waddell D, Abovsky M, Ahmed Z, Kotlyar M, Lu R, Jurisica I (2023). mirDIP 5.2: tissue context annotation and novel microRNA curation. *Nucleic Acids Res* 51(D1):D217–D225. doi:10.1093/nar/gkac1070
- Shirdel EA, Xie W, Mak TW, Jurisica I (2011). NAViGaTing the Micronome — using multiple miRNA prediction databases. *PLoS ONE* 6(2):e17429. doi:10.1371/journal.pone.0017429
- mirDIP portal and statistics (version 5.3.0.3 / DB 5.2.2.1): https://ophid.utoronto.ca/mirDIP
- GitHub: https://github.com/ijlab/mirdip

---

## Integration notes (for consensus project)

- Treat mirDIP as **meta-layer only** — do not add TargetScan + DIANA + miRDB + mirDIP as four independent votes.
- Preferred use: **Score Class ≥ High** as a tier-1 filter; export Integrated Score for rank fusion **instead of** constituent tool scores when mirDIP is selected.
- **Number of Sources** can annotate support breadth but is not calibrated for averaging.
- For PTEN: mirDIP may rank a pair Very High because TargetScan, DIANA, and miRDB all predict it — this is **correlated evidence**, not four independent confirmations.
