# miRDB

**Status:** complete

Use [template-algorithm-profile.md](template-algorithm-profile.md). Deep-dive checklist: [deep-dive-checklists.md](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md#mirdb).

Research plan: [RESEARCH_PLAN.md](RESEARCH_PLAN.md#4-mirdb)

Provides ML-based contrast to rule/physics methods (TargetScan, DIANA).

---

## 1. What is it?

### Who developed it?

**miRDB** and the **MirTarget** prediction algorithm were developed by **Xiaowei Wang's lab** (Department of Pharmacology, University of Illinois at Chicago). Key publications by Weijun Liu, Yuhao Chen, and Xiaowei Wang.

### When was it introduced?

- **2008** — miRDB database (original release)
- **2015** — MirTarget v3 with SVM framework (Wong & Wang, *Nucleic Acids Res.*)
- **2019** — **MirTarget v4.0** integrative CLIP + RNA-seq model (Liu & Wang, *Genome Biology*)
- **2020** — **miRDB v6.0** database update (Chen & Wang, *Nucleic Acids Res.*)

### Current / relevant version

**miRDB v6.0** — ~3.5 million predicted targets for ~7,000 miRNAs across **human, mouse, rat, dog, and chicken** (Chen & Wang, 2020). URL: https://mirdb.org/

### What problem does it solve?

miRDB provides **genome-wide miRNA target predictions** using a **machine-learning** approach (MirTarget) trained on high-throughput **miRNA binding** (CLIP ligation) and **target downregulation** (RNA-seq) data — aiming to predict **functional** targets, not merely binding sites (Liu & Wang, 2019; Chen & Wang, 2020).

---

## 2. How does it work?

```
miRNA + RefSeq 3′ UTR
→ scan for seed8A1, seed7b, seed7A1 sites (canonical)
→ compute 96 features per candidate site
→ SVM (RBF kernel) → site probability Pᵢ ∈ [0,1]
→ combine sites: S = 100 × (1 − ∏Pᵢ) → gene–miRNA MirTarget score
→ retain targets with score ≥ 50
```

### Step-by-step process

1. **Input sequences** — miRNAs from miRBase v22; 3′ UTRs from NCBI RefSeq (Chen & Wang, 2020).
2. **Candidate site identification** — Search for **seed8A1, seed7b, and seed7A1** canonical sites in 3′ UTRs (Liu & Wang, 2019). **Conservation is not required** for candidacy (Chen & Wang, 2020).
3. **Feature extraction** — **96 features** per site, including seed conservation (phyloP and cross-species word search), UTR length, site position, GC content, seed binding stability, and others (Liu & Wang, 2019, Additional file 3: Table S2).
4. **Site scoring** — **Support vector machine (SVM)** with radial basis function (RBF) kernel; libsvm outputs site probability **Pᵢ ∈ [0, 1]** (Liu & Wang, 2019).
5. **Gene-level aggregation** — **S = 100 × (1 − ∏ᵢ Pᵢ)** across all candidate sites in the 3′ UTR (Liu & Wang, 2019).
6. **Filter and rank** — Present gene–miRNA pairs with **score ≥ 50**, ordered by score descending (Chen & Wang, 2020; miRDB FAQ).

### Candidate site identification

- **3′ UTR only** (RefSeq); no CDS scanning (Chen & Wang, 2020).
- Three seed types: **seed8A1, seed7b, seed7A1** — selected as strongest enrichments in RNA-seq downregulation data (Liu & Wang, 2019, Table 2).
- Both **conserved and non-conserved** sites eligible (Chen & Wang, 2020).

### Scoring and ranking

- **MirTarget score** (0–100) at gene–miRNA level; higher = more confident (miRDB FAQ).
- Search results ordered by prediction score (miRDB FAQ).

---

## 3. What biological/computational features does it use?

| Feature | Used? | Source / notes |
|---------|-------|----------------|
| Seed matching | ✅ | seed8A1, seed7b, seed7A1 canonical types — Liu & Wang, 2019 |
| 8mer / 7mer-m8 / 7mer-A1 / 6mer site types | Partial | Uses seed8A1 ≈ 8mer-A1, seed7b ≈ 7mer-m8, seed7A1 ≈ 7mer-A1; **seed6 excluded** from modeling — Liu & Wang, 2019 |
| Evolutionary conservation | ✅ | phyloP 100-way conservation + cross-species seed match; **top RFE feature: conserved seed8A1** — Liu & Wang, 2019 |
| Binding / free energy | Partial | Seed binding stability (ΔG); **not significant when CLIP + RNA-seq combined** (p = 0.26) — Liu & Wang, 2019 |
| RNA / site accessibility | Partial | GC content / G nucleotide depletion; not explicit structure model — Liu & Wang, 2019 |
| Local AU content | Not documented | — |
| Target-site abundance | ✅ | Multi-site aggregation formula combines all sites in 3′ UTR — Liu & Wang, 2019 |
| Position within transcript | ✅ | Distance to 3′ UTR ends and transcript center — Liu & Wang, 2019 |
| 3′ compensatory pairing | Not documented | Focus on canonical seed types only |
| G:U pairing | Not documented | — |
| Machine learning | ✅ | **SVM (RBF kernel)**, 96 features, MirTarget v4.0 — Liu & Wang, 2019 |
| Other sequence/context features | ✅ | 3′ UTR length; seed-level GC; CLIP- vs overexpression-specific G bias — Liu & Wang, 2019 |

---

## 4. What does its score mean?

| Question | Answer |
|----------|--------|
| Score name(s) | **MirTarget score** (a.k.a. target prediction score in miRDB) |
| Numerical range | **50–100** displayed (scores ≥ 50 retained); underlying site probabilities 0–1 — Chen & Wang, 2020; miRDB FAQ |
| Higher or lower = stronger? | **Higher = stronger** |
| Score type | SVM-derived **ML confidence** score (probability combination) |
| Absolute interpretation across genes? | Partial — designed for ranking; score reflects statistical confidence of functional targeting |
| Already normalized? | Site probabilities from libsvm; gene score scaled to 0–100 via combination formula |
| Recommended cutoff / threshold | **≥ 50** minimum displayed; **> 80** high confidence; **< 60** cautious — miRDB FAQ |
| Site scores aggregated to gene–miRNA score? | **Yes** — **S = 100 × (1 − ∏Pᵢ)** (Liu & Wang, 2019) |

### Score interpretation notes

- MirTarget score reflects **confidence that the gene is a functional miRNA target** (binding + downregulation evidence in training), not predicted repression magnitude.
- **ML confidence** evidence type — distinct from TargetScan efficacy (Context++) and DIANA GLM efficacy.
- Non-conserved targets included — miRDB reports **more targets per miRNA** than prior versions (606 human targets/miRNA in v6.0; Chen & Wang, 2020).

### 4b. Score inventory

| Score name | Level | Evidence type | Range | Direction | Notes |
|------------|-------|---------------|-------|-----------|-------|
| Site probability (Pᵢ) | site | ML confidence | 0–1 | Higher = stronger | libsvm output per candidate site — Liu & Wang, 2019 |
| MirTarget score | gene–miRNA | ML confidence | 50–100 (displayed) | Higher = stronger | Combined via noisy-or-like formula × 100 — Liu & Wang, 2019 |

---

## 5. Advantages

- **Integrative training** — combines CLIP binding + RNA-seq downregulation, targeting **functional** repression not just binding (Liu & Wang, 2019).
- **Strong independent validation** — best ROC AUC (0.78) on HITS-CLIP test vs TargetScan, DIANA-MicroT, miRanda, PITA (Liu & Wang, 2019).
- **Non-conserved target discovery** — conservation not required for candidacy (Chen & Wang, 2020).
- **Multi-species** — human, mouse, rat, dog, chicken (Chen & Wang, 2020).
- **Cell-line expression filter** — v6.0 adds expression profiles to restrict targets to expressed genes (Chen & Wang, 2020).
- **Custom prediction** — users can submit sequences for MirTarget scoring (miRDB website).

---

## 6. Disadvantages

- **Black-box ML** — 96-feature SVM less interpretable than TargetScan site-type taxonomy or DIANA binding categories.
- **Training data dependence** — 25 miRNAs in RNA-seq arm; CLIP datasets from specific cell types (Liu & Wang, 2019).
- **3′ UTR only** — no CDS targeting (contrast DIANA-microT-CDS).
- **Limited seed types** — only three canonical seeds modeled; seed6 and non-canonical sites excluded (Liu & Wang, 2019).
- **Overexpression concern** — RNA-seq training uses transfected miRNAs; may not reflect physiological targeting (Liu & Wang, 2019 Discussion).
- **Score not comparable** to Context++, microT, or mirDIP integrated score without normalization.

---

## 7. When would you use it?

Use miRDB when you want a **machine-learning consensus** on functional 3′ UTR targeting with explicit confidence tiers — especially to capture **non-conserved** or **ML-weighted** predictions that rule-based tools may deprioritize. For PTEN: miRDB may report targets based on feature patterns (conservation, UTR context, GC) weighted differently from TargetScan's Context++ or DIANA's CDS-aware GLM — explaining score and membership disagreements.

---

## 8. Method-specific deep dive

- [x] **How the MirTarget machine-learning model works** — SVM (RBF kernel, libsvm); 96 features; 5-fold cross-validation for hyperparameters; MirTarget v4.0 (Liu & Wang, 2019).
- [x] **What training data it uses** — **CLIP ligation** (Helwak + Grosswendt studies) + **RNA-seq** overexpression of 25 miRNAs in HeLa (1.5B reads, 52 samples); positives = ≥40% downregulation with seed site; negatives = stable expression (Liu & Wang, 2019).
- [x] **What features it considers** — 96 features including conserved seed8A1 (top RFE), phyloP conservation, cross-species seed match, 3′ UTR length, site position, GC/G content; full list in Additional file 3: Table S2 (Liu & Wang, 2019).
- [x] **What its score means** — ML confidence that gene is a **functional miRNA target**; higher = more likely real (miRDB FAQ; Liu & Wang, 2019).
- [x] **Its score range and cutoff** — Displayed **50–100**; **> 80** high confidence, **< 60** cautious (miRDB FAQ); minimum display cutoff **50** (Chen & Wang, 2020).
- [x] **Strengths and weaknesses as ML contrast to TargetScan/DIANA** — Strengths: integrative functional training, non-conserved targets, strong CLIP/knockdown validation. Weaknesses: opaque SVM, 3′ UTR-only, limited seed set, overexpression training bias (Liu & Wang, 2019).

### Site-level → gene–miRNA aggregation

**S = 100 × (1 − ∏ᵢ₌₁ⁿ Pᵢ)** where *n* = number of candidate sites in the 3′ UTR and *Pᵢ* = SVM site probability (Liu & Wang, 2019). For single-site targets, gene score = 100 × P₁.

---

## Sources

- Liu W, Wang X (2019). Prediction of functional microRNA targets by integrative modeling of microRNA binding and target expression data. *Genome Biology* 20:18. doi:10.1186/s13059-019-1629-z
- Chen Y, Wang X (2020). miRDB: an online database for prediction of functional microRNA targets. *Nucleic Acids Res* 48(D1):D127–D131. doi:10.1093/nar/gkz757
- miRDB FAQ (score interpretation and thresholds): https://mirdb.org/mirdb/faq.html
- miRDB website: https://mirdb.org/

---

## Integration notes (for consensus project)

- Export **MirTarget score** as **ML confidence** rank; threshold **≥ 60** (default mining) or **> 80** (high confidence) per miRDB FAQ.
- Do not average MirTarget score with TargetScan CWCS or DIANA microT — use **percentile normalization within source** or **rank fusion**.
- miRDB is also a **constituent of mirDIP** — if mirDIP is used, do not double-count miRDB separately at equal weight.
