# Presentation � miRNA Algorithm Comparison

Paste-ready slide deck content answering: **Why do TargetScan, DIANA, miRDB, and mirDIP disagree for the same gene?**

## Files

| File | Purpose |
|------|---------|
| [`slide-content.md`](slide-content.md) | Bullets for slides 1�16 (15 content + references) |
| [`speaker-notes.md`](speaker-notes.md) | Delivery notes and timing guide |
| [`assets/`](assets/) | SVG diagrams for key slides |
| [`sources/bibliography.md`](sources/bibliography.md) | Citation index |

## Research backing

All score ranges and methodological claims are sourced from completed research profiles:

- [`../research/targetscan.md`](../research/targetscan.md)
- [`../research/diana-microt-cds.md`](../research/diana-microt-cds.md)
- [`../research/mirdb.md`](../research/mirdb.md)
- [`../research/mirdip.md`](../research/mirdip.md)
- [`../research/validation-databases.md`](../research/validation-databases.md)
- [`../research/score-comparison-matrix.md`](../research/score-comparison-matrix.md)

## Building the deck

1. Create one slide per `## Slide N` heading in `slide-content.md`
2. Insert SVG assets where marked (**Visual:** lines)
3. Slide 14 uses five pros/cons bullets (one per tool family) � no table needed
4. Use slide 16 as the References slide

## Running example

**PTEN** � slide 13 and project smoke tests (`find_mirna_regulators("PTEN")`).
