---
name: mirna-presentation
description: >-
  Draft the miRNA algorithm comparison slideshow answering why predictors
  disagree. Use when editing docs/presentation/, slide-content.md, speaker
  notes, or presentation assets for the TargetScan/DIANA/miRDB/mirDIP deck.
---

# miRNA Presentation

## When to Use

Apply when drafting slides, speaker notes, or visuals for the algorithm comparison presentation.

## Central Question

**Why do TargetScan, DIANA, miRDB, and mirDIP give different predictions for the same gene?**

Every slide should support answering this with source-backed differences in methodology and scoring.

## Slide Arc (15 content slides)

See [slide-outline.md](slide-outline.md) for titles and core messages.

| Slide | Focus |
|-------|-------|
| 1–4 | Background: difficulty, miRNA biology, prediction vs validation, general pipeline |
| 5–6 | TargetScan overview + scoring (Context++ vs PCT) |
| 7–8 | DIANA overview + scoring |
| 9 | miRDB (ML contrast) |
| 10 | mirDIP (integrative, not single algorithm) |
| 11 | Validation DBs (miRTarBase, TarBase) |
| 12 | Side-by-side comparison matrix |
| 13 | **Centerpiece:** why algorithms disagree (PTEN example) |
| 14 | Pros/cons by method |
| 15 | Implications for consensus scoring |

## Content Rules

- Max **5 bullets per slide**
- Paste-ready content → [`docs/presentation/slide-content.md`](../../../docs/presentation/slide-content.md)
- One `## Slide N` heading per slide in slide-content.md
- Detailed score tables on slides 6 and 8 only
- Slide 12 pulls from [`docs/research/score-comparison-matrix.md`](../../../docs/research/score-comparison-matrix.md) — do not invent ranges
- Slide 15: explain why naive averaging fails; mention rank fusion, percentile norm, tiered evidence
- Every non-obvious claim → citation in [`docs/presentation/sources/bibliography.md`](../../../docs/presentation/sources/bibliography.md)

## PTEN Running Example

Use PTEN on slide 13 to show:

- Different tools may rank different miRNAs at the top
- Disagreement reflects different evidence, not necessarily error
- Tie explanations to specific methodological differences (seed rules, conservation, ML, integration)

## Visuals

Place diagrams and comparison tables in `docs/presentation/assets/`.

## GHLC Extension

If presenting for GHLC Student Speaker Series, add template/rubric files under `docs/presentation/sources/` and extend this skill with time limits and bookend slides (Title → Introduction → content → References → Acknowledgements → Thank You → Contact).

## Additional Resources

- Slide outline: [slide-outline.md](slide-outline.md)
- Project PDF: [`miRNA project.pdf`](../../../miRNA%20project.pdf)
