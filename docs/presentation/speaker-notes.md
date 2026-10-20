# Speaker Notes � miRNA Algorithm Comparison

Brief delivery notes for each slide. Pair with [`slide-content.md`](slide-content.md) and [`assets/`](assets/).

---

## Slide 1: Why miRNA Target Prediction Is Hard

Open with the central tension: miRNAs bind imperfectly, yet regulate many genes. The combinatorial explosion makes exhaustive validation impossible � every tool is a prioritization heuristic. Frame the deck's question: why do TargetScan, DIANA, miRDB, and mirDIP disagree for the same gene?

---

## Slide 2: miRNAs and Targets

Define seed matching before naming any tool. Emphasize that most predictors scan 3? UTRs; DIANA is the outlier with CDS. One gene (PTEN) can have sites for many miRNAs � ranking is the real task.

---

## Slide 3: Prediction vs Experimental Validation

**Show diagram:** `slide03-prediction-vs-validation.svg`

Critical distinction for slides 11 and 15. Predictors produce scores; validation DBs produce evidence categories. Our R tool ranks first, then flags `experimental_support` � never averages TarBase into TargetScan scores.

---

## Slide 4: How Predictors Generally Work

**Show diagram:** `slide04-predictor-pipeline.svg`

Walk left-to-right once, then note that each tool customizes every box. This sets up slides 5�10 as variations on the same pipeline.

---

## Slide 5: TargetScan Overview

TargetScan is the reference seed-based method. Mention miRNA families (shared seed) � users query one miRNA but TargetScan may report family-level predictions. TargetScan 8.0 changed default ranking, not site inventory.

---

## Slide 6: TargetScan Scoring System

**Show diagram:** `slide06-targetscan-scoring-axes.svg`

Spend time here � this is the most common integration mistake. CWCS and PCT answer different biological questions. If audience asks "which score to use": CWCS for cross-tool efficacy comparison; PCT as conservation annotation only.

---

## Slide 7�8: DIANA-microT-CDS

Highlight CDS scanning and non-canonical sites as the main differentiators from TargetScan. microT is a single efficacy GLM score � conservation is embedded, not separate like PCT. Threshold 0.8 is the practical default for human.

---

## Slide 9: miRDB

Contrast ML (black box) with rule-based methods. MirTarget was trained on both binding (CLIP) and downregulation (RNA-seq) � closer to "functional target" than CLIP-only training. Score >80 is the practical high-confidence cutoff.

---

## Slide 10: mirDIP

Stress: **integrative, not independent.** mirDIP's value is coverage and calibrated Score Class � not a fifth equal vote alongside TargetScan. Mention that mirDIP includes TargetScan v7.2, DIANA, and miRDB internally.

---

## Slide 11: Validation Databases

Neither miRTarBase nor TarBase predicts new interactions. Strong evidence (reporter, Western, qPCR) vs weak (CLIP, microarray) matters for how much trust to place on the overlay. TarBase v9 adds cell-type resolution � useful when tissue context matters.

---

## Slide 12: Side-by-Side Comparison

**Show diagram:** `slide12-evidence-types.svg`

Quick recap table � point to the "Comparable? No" column. Every score row in the research matrix is documented in `docs/research/score-comparison-matrix.md`.

---

## Slide 13: Why Algorithms Disagree (PTEN Example)

**Show diagram:** `slide13-pten-disagreement.svg`

Centerpiece slide. PTEN is ideal because miR-21 and miR-214 are **experimentally validated** regulators � yet tools still rank different miRNAs at the top and assign incompatible numeric scores. Disagreement is methodological, not necessarily biological error. Avoid claiming one tool's top hit is "correct" unless validation supports that specific pair.

---

## Slide 14: Advantages and Disadvantages

Walk the five bullets tool-by-tool. No tool wins on every dimension � choice depends on whether you need CDS coverage, interpretability, ML recall, or meta-integration breadth.

---

## Slide 15: Implications for Consensus Scoring

**Show diagram:** `slide15-consensus-scoring-flow.svg`

Close with the architecture our R modules will implement: per-source ranks ? fusion ? validation overlay. Explicitly reject naive averaging. mirDIP sits above or beside fusion, not inside it as a duplicate vote.

---

## Slide 16: References

Point audience to bibliography file and research docs under `docs/research/` for full methodology profiles.

---

## Timing guide (approximate)

| Slides | Topic | Minutes |
|--------|-------|---------|
| 1�4 | Background | 5�6 |
| 5�6 | TargetScan | 4�5 |
| 7�8 | DIANA | 3�4 |
| 9�10 | miRDB + mirDIP | 4 |
| 11�12 | Validation + comparison | 3�4 |
| 13 | PTEN centerpiece | 4�5 |
| 14�15 | Pros/cons + consensus | 4�5 |
| 16 | References | 1 |

**Total:** ~30�35 minutes with discussion
