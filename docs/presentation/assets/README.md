# Presentation Assets

SVG diagrams for the miRNA algorithm comparison deck. Import into PowerPoint, Google Slides, or Keynote (Insert ? Image).

All assets use a 960px-wide canvas suitable for 16:9 slides.

## Asset catalog

| File | Slide | Description |
|------|-------|-------------|
| [`slide03-prediction-vs-validation.svg`](slide03-prediction-vs-validation.svg) | 3 | Two-column: predictors vs validation DBs; annotation arrow |
| [`slide04-predictor-pipeline.svg`](slide04-predictor-pipeline.svg) | 4 | Six-step predictor pipeline with per-tool variation notes |
| [`slide06-targetscan-scoring-axes.svg`](slide06-targetscan-scoring-axes.svg) | 6 | Efficacy axis (Context++/CWCS/occupancy) vs PCT conservation axis |
| [`slide12-evidence-types.svg`](slide12-evidence-types.svg) | 12 | Table: resource, score, evidence type, region, comparability |
| [`slide13-pten-disagreement.svg`](slide13-pten-disagreement.svg) | 13 | PTEN hub with four tool branches + validated miR-21/miR-214 callouts |
| [`slide15-consensus-scoring-flow.svg`](slide15-consensus-scoring-flow.svg) | 15 | Per-source rank ? fusion ? consensus + validation overlay |

## Design notes

- Colors encode tool families: blue (TargetScan), cyan (DIANA), purple (miRDB), orange (mirDIP), green (validation)
- Diagrams contain only source-backed labels � no invented numeric scores from live API queries
- PTEN slide cites validated interactions (miR-21, miR-214) from published literature; see bibliography EX-001, EX-002

## Related files

- Paste-ready bullets: [`../slide-content.md`](../slide-content.md)
- Speaker notes: [`../speaker-notes.md`](../speaker-notes.md)
- Score matrix (source for slide 12): [`../../research/score-comparison-matrix.md`](../../research/score-comparison-matrix.md)
- Bibliography: [`../sources/bibliography.md`](../sources/bibliography.md)

