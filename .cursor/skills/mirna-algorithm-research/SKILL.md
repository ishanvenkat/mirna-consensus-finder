---
name: mirna-algorithm-research
description: >-
  Research miRNA target-prediction algorithms and validation databases using
  the 7-question template and PDF deep-dive checklists. Use when writing
  docs/research/*.md, comparing prediction scores, explaining why
  TargetScan/DIANA/miRDB/mirDIP disagree, or updating the score comparison matrix.
---

# miRNA Algorithm Research

## When to Use

Apply when researching TargetScan, DIANA-microT-CDS, miRDB, mirDIP, miRTarBase, or TarBase; drafting algorithm profiles; or preparing score-integration rationale.

## Workflow

1. Read [`docs/research/RESEARCH_PLAN.md`](../../../docs/research/RESEARCH_PLAN.md) for task order and sources
2. Read [`docs/research/template-algorithm-profile.md`](../../../docs/research/template-algorithm-profile.md)
3. Complete the resource checklist in [`deep-dive-checklists.md`](deep-dive-checklists.md)
4. Fill the corresponding stub under `docs/research/`
5. Consult [`reference-sources.md`](reference-sources.md) for canonical URLs and papers
6. Add bibliography entries to [`docs/presentation/sources/bibliography.md`](../../../docs/presentation/sources/bibliography.md)
7. Update [`docs/research/score-comparison-matrix.md`](../../../docs/research/score-comparison-matrix.md)

**Before marking a doc `complete`:** all Section 8 deep-dive items checked or explicitly marked not documented.

## Core Distinctions

- **Prediction:** TargetScan, DIANA-microT-CDS, miRDB, mirDIP
- **Validation (not prediction):** miRTarBase, DIANA-TarBase
- **mirDIP:** integrative meta-resource — not a single algorithm like TargetScan
- **TargetScan:** Context++ (targeting efficacy) and PCT (conservation evidence) are **separate score types**

## Verification Rules

- Tools4miRs is an **index only** — verify every feature in official docs or original papers
- Do not assume feature usage without a source
- Do not invent score ranges or cutoffs

## Score Documentation Checklist

For each score (Section 4b in template):

- Name and level (site vs gene–miRNA pair)
- Evidence type: efficacy | conservation | ML confidence | meta-integration
- Numerical range and direction
- Type (probability, confidence, rank, energy, etc.)
- Absolute vs relative interpretation across genes
- Pre-existing normalization
- Recommended cutoff/threshold
- Site-to-gene aggregation method
- Comparability to other tools' scores

## TargetScan Priority Topics

See [deep-dive-checklists.md#targetscan](deep-dive-checklists.md#targetscan). Emphasize:

1. Context++ mathematical logic (site and gene level)
2. Site-level → gene–miRNA aggregation
3. PCT — independent from Context++
4. Seed site types and ranking logic

## Output Quality

- Tag interpretations vs documented facts clearly
- End each research doc with a **Sources** section
- Use PTEN as the running example when illustrating disagreements

## Additional Resources

- Deep-dive checklists: [deep-dive-checklists.md](deep-dive-checklists.md)
- Research roadmap: [RESEARCH_PLAN.md](../../../docs/research/RESEARCH_PLAN.md)
- Canonical URLs: [reference-sources.md](reference-sources.md)
- Project context: [`miRNA project.pdf`](../../../miRNA%20project.pdf)
