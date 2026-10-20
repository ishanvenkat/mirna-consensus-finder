# Slide Outline — 15 Content Slides

Paste-ready bullets go in `docs/presentation/slide-content.md`.

## Slide 1: Why miRNA Target Prediction Is Hard

- Imperfect miRNA–mRNA complementarity
- Thousands of miRNAs × tens of thousands of genes
- Experimental validation of all pairs is infeasible
- Algorithms must prioritize using partial sequence/context evidence

## Slide 2: miRNAs and Targets

- miRNAs: small non-coding RNAs regulating gene expression
- Seed region (nt 2–7/8) critical for target recognition
- Target sites often in 3′ UTR; some tools also use CDS
- Relationship: miRNA → target site → mRNA/gene

## Slide 3: Prediction vs Experimental Validation

- **Predictors:** TargetScan, DIANA-microT-CDS, miRDB, mirDIP
- **Validation DBs:** miRTarBase, DIANA-TarBase
- Predictors ask: could this miRNA regulate this gene?
- Validation DBs ask: has this interaction been experimentally supported?

## Slide 4: How Predictors Generally Work

- Input: miRNA + mRNA/transcript sequences
- Find candidate target sites (seed matching)
- Score using biological/computational features
- Aggregate site scores → gene–miRNA score
- Rank and output predictions

## Slide 5: TargetScan Overview

- Major seed-based predictor; conservation-aware
- Canonical site types: 8mer, 7mer-m8, 7mer-A1, 6mer
- Current version and key developers *(fill from research)*

## Slide 6: TargetScan Scoring System

- Context++: predicted targeting efficacy
- PCT: probability of conserved targeting
- Weighted Context++ / aggregate gene–miRNA scores
- **These are different evidence types — not interchangeable**

## Slide 7: DIANA-microT-CDS Overview

- Predicts sites in 3′ UTR and CDS
- Combines thermodynamic, conservation, accessibility, ML features
- microT-CDS naming and version *(fill from research)*

## Slide 8: DIANA Scoring System

- microT score: meaning, range, direction
- Recommended threshold *(fill from research)*
- How multiple sites aggregate to gene–miRNA score

## Slide 9: miRDB

- MirTarget machine-learning predictor
- Training data and features *(fill from research)*
- ML-based contrast to rule/physics methods

## Slide 10: mirDIP

- Integrative database — combines many prediction resources
- Not a single algorithm like TargetScan
- Integrated score + per-source confidence
- Sources integrated *(fill from research)*

## Slide 11: Validation Databases

- miRTarBase: experimentally supported interactions
- DIANA-TarBase: curated experimental evidence
- Evidence types: reporter assays, qPCR, CLIP, knockdown/overexpression
- Used as reference overlay, not for prediction

## Slide 12: Side-by-Side Comparison

- Feature matrix from `docs/research/score-comparison-matrix.md`
- Compare seed rules, conservation, ML, CDS support, score semantics

## Slide 13: Why Algorithms Disagree (PTEN Example)

- Same gene, different top-miRNA lists across tools
- Differences in: candidate-site selection, seed rules, conservation models
- Thermodynamics, accessibility, ML training data, aggregation methods
- Disagreement ≠ one tool is simply wrong

## Slide 14: Advantages and Disadvantages

- Source-backed pros/cons per method
- Avoid generic limitations without citations

## Slide 15: Implications for Consensus Scoring

- Scores are not directly comparable — cannot naively average
- Options: rank fusion, percentile normalization, tiered evidence
- mirDIP already meta-integrates — avoid double-counting
- Validation overlay separate from prediction scores
