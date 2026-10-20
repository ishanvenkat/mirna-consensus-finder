# Score Comparison Matrix

Fill one row per score type when research is executed. Do not invent values — mark `TBD` until verified.

| Resource | Score name | Level (site / gene–miRNA) | Evidence type | Range | Direction | Score type | Absolute across genes? | Normalized? | Cutoff | Aggregation method | Comparable to others? |
|----------|------------|---------------------------|---------------|-------|-----------|------------|------------------------|-------------|--------|-------------------|---------------------|
| TargetScan | Context++ | site | efficacy | Site-type bounded (8mer ≤ −0.03, 7mer-m8 ≤ −0.02, 7mer-A1 ≤ −0.01, 6mer ≤ 0) | More negative = stronger | Linear regression contribution (log repression) | Partial — calibrated but miRNA-family context varies | Features trimmed-normalized; percentiles in downloads | None recommended | N/A (site-level) | No — TargetScan-specific scale |
| TargetScan | Total context++ score | gene–miRNA | efficacy | Unbounded (sum) | More negative sum = stronger | Sum of site Context++ scores | Partial | Derived from normalized features | None recommended | Sum of non-overlapping sites (maximize total) | No |
| TargetScan | Weighted Context++ (CWCS) | gene–miRNA | efficacy | log₂ scale; typically ≤ 0 in 8.0 Summary Counts | **More negative** = stronger | log₂ predicted fold-change in target expression | Partial | AIR-weighted aggregation | None recommended | AIR-weighted marginal repression formula (Agarwal 2015) | No |
| TargetScan | PCT / Aggregate PCT | site; aggregate at gene–miRNA | conservation | 0–1 | Higher = stronger | Probability | Yes (probability scale) | Inherently 0–1 | None recommended | Aggregate PCT rollup for ranking | No — measures conservation, not efficacy |
| TargetScan | Biochemical predicted occupancy | gene–miRNA | efficacy | ~0–1 (fraction) | Higher = stronger | Predicted AGO–miRNA occupancy → repression | Depends on miRNA concentration scenario | Three scenarios in downloads (low/high/transfected) | None recommended | Biochemical model from site *K*<sub>D</sub> values (McGeary 2019) | No — different framework from Context++ despite same evidence type |
| DIANA-microT-CDS | microT score | gene–miRNA | efficacy | Cutoff 0.3–1 (full range not documented) | Higher = stronger | GLM-derived interaction score | Partial — ranking-oriented | Per-resource GLM output | Human: 0.7/0.8/0.9 (DIANA API) | MRE sum → region sum → final GLM | No |
| DIANA-microT-CDS | MRE score | site | efficacy | Not documented | Higher = stronger | GLM site contribution | Partial | GLM feature-normalized | N/A (site-level) | Summed within CDS or 3′ UTR region | No |
| miRDB | MirTarget score | gene–miRNA | ML confidence | 50–100 (displayed) | Higher = stronger | SVM confidence (combined site probabilities) | Partial — ranking-oriented | Scaled from libsvm probabilities | ≥50 display; >80 high confidence; <60 cautious (FAQ) | S = 100 × (1 − ∏Pᵢ) across 3′ UTR sites | No |
| miRDB | Site probability (Pᵢ) | site | ML confidence | 0–1 | Higher = stronger | libsvm probability | Partial | SVM output | N/A (site-level) | Combined to gene score via product formula | No |
| mirDIP | Integrated score | gene–miRNA | meta-integration | 0–1 (median 0.02) | Higher = stronger | Noisy-or of calibrated precision scores | Partial — probability-like | Rank-normalized + precision-calibrated per resource | Use Score Class (Very High / High / Medium / Low) | Noisy-or across 24 resources | No — already meta-integrated |
| mirDIP | Score class | gene–miRNA | meta-integration | 4 categories | Very High > High > Medium > Low | Per-miRNA percentile tier | Relative within miRNA | Assigned from integrated score rank | Very High = top 1% per miRNA (v5.2) | Derived from integrated score | No |
| miRTarBase | Support level | gene–miRNA | experimental support | Strong / Weak | Strong > Weak | Categorical evidence tier | N/A | Method-based assignment | Filter to Strong for high-confidence overlay | N/A — not a score | No — overlay only |
| miRTarBase | Validation method | gene–miRNA | experimental support | Categorical (reporter, Western, qPCR, CLIP, microarray, etc.) | Reporter/Western/qPCR = strong tier | Experimental technique label | N/A | Per-interaction curation | User selects method filter | N/A | No |
| DIANA-TarBase | Experimental method | gene–miRNA | experimental support | 37 protocols (22 low-yield + 15 HT) | Context-dependent | Protocol name | N/A | Curated per entry | Filter by method in TarBase v9 UI | N/A — not a score | No — overlay only |
| DIANA-TarBase | Protocol type | gene–miRNA | experimental support | Direct / Indirect | N/A | CLIP-family = direct binding evidence | N/A | Assigned per interaction | Web filter — Skoufos et al., 2024 | N/A | No |

**Evidence types:** `efficacy` | `conservation` | `ML confidence` | `meta-integration` | `experimental support`

Context++ and PCT measure different things — do not treat as interchangeable.

## Why Naive Averaging Fails

**Exemplar: TargetScan scores for the same gene–miRNA pair**

Consider PTEN targeted by a miRNA family with one strong 8mer site:

| Score | Value (illustrative) | What it measures |
|-------|---------------------|------------------|
| Context++ (site) | −0.025 | Predicted site-level repression (negative = strong) |
| CWCS | −1.8 | Predicted cumulative gene-level repression (log₂; more negative = stronger) |
| Aggregate PCT | 0.92 | Probability of conserved targeting |
| Predicted occupancy | 0.15 | Fraction of AGO–miRNA on PTEN sites |

Naive averaging fails because:

1. **Different evidence types** — CWCS/occupancy measure predicted **biochemical repression**; PCT measures **evolutionary conservation**. Averaging them treats "is this site conserved because it is targeted?" the same as "how much will this miRNA repress PTEN?" — biologically distinct questions.
2. **Incompatible scales and directions** — Context++ and CWCS share a log-scale repression convention (more negative = stronger) but differ in level and aggregation; PCT is 0–1 probability; occupancy is 0–1 fraction with higher = stronger. A mean of (−0.025 + 0.92 + (−1.8) + 0.15) / 4 is meaningless.
3. **PCT is embedded in Context++** — PCT contributes as one of 14 Context++ features, so averaging standalone PCT with Context++ partially double-counts conservation signal.
4. **Release/version hazard** — TargetScan 8.0 defaults to biochemical occupancy; TargetScan 7.x defaulted to CWCS. Averaging across releases or ranking options mixes models trained on different frameworks.
5. **Concentration dependence** — Biochemical occupancy varies by assumed miRNA concentration (low/high/transfected scenarios); PCT and Context++ do not share this parameterization.

**Better approach for this project:** use CWCS (or occupancy, consistently) as the efficacy rank; attach aggregate PCT as a separate conservation tier; do not merge with DIANA/miRDB/mirDIP scores without rank fusion or percentile normalization within each source.

## Integration Options (for Phase 4)

1. Per-source rank fusion (e.g., Borda count)
2. Percentile normalization within source
3. Tiered evidence (mirDIP as meta-layer)
4. Validation overlay (miRTarBase/TarBase) — separate from prediction scores

## Related docs

- Research roadmap: [RESEARCH_PLAN.md](RESEARCH_PLAN.md)
- Deep-dive checklists: [deep-dive-checklists.md](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md)
- TargetScan profile: [targetscan.md](targetscan.md)
- DIANA-microT-CDS profile: [diana-microt-cds.md](diana-microt-cds.md)
- miRDB profile: [mirdb.md](mirdb.md)
- mirDIP profile: [mirdip.md](mirdip.md)
- Validation databases: [validation-databases.md](validation-databases.md)
