# PDF-Mandated Deep-Dive Checklists

Single source of truth for method-specific research topics from [`miRNA project.pdf`](../../../miRNA%20project.pdf). Each item must be verified from primary sources before marking complete in the corresponding stub under `docs/research/`.

Execution roadmap: [RESEARCH_PLAN.md](../../../docs/research/RESEARCH_PLAN.md)

---

## TargetScan

Stub: [targetscan.md](../../../docs/research/targetscan.md)

- [x] Seed matching — Lewis et al., 2005; FAQs
- [x] 8mer sites — definition, seed-match rules, relative efficacy — FAQs; Friedman et al., 2009
- [x] 7mer-m8 sites — FAQs; Agarwal et al., 2015
- [x] 7mer-A1 sites — FAQs; Lewis et al., 2005
- [x] 6mer sites — FAQs; always nonconserved in TargetScanHuman 7+
- [x] Context++ score — mathematical logic (context features, how combined) — Agarwal et al., 2015; context_score.html
- [x] Weighted Context++ score — context_score_totals.html; Agarwal et al., 2015
- [x] Aggregate/cumulative scores (site → gene–miRNA rollup) — context_score_totals.html
- [x] PCT (Probability of Conserved Targeting) — separate from Context++ — Friedman et al., 2009
- [x] Evolutionary conservation model — FAQs branch-length thresholds
- [x] How TargetScan ranks predictions (gene pages vs download files) — vert_80 data download
- [x] Version changes (current vs prior releases) — TargetScan 8.0 help page
- [x] How different TargetScan scores should be interpreted — FAQs Q11–12
- [x] McGeary 2019 biochemical model (8.0 default ranking) — McGeary et al., 2019; vert_80 home page

**Critical:** Context++ = targeting efficacy; PCT = conservation evidence. Do not conflate.

**Critical:** Document site-level → gene–miRNA aggregation (Agarwal et al. 2015, TargetScan FAQs).

---

## DIANA-microT-CDS

Stub: [diana-microt-cds.md](../../../docs/research/diana-microt-cds.md)

- [x] What microT-CDS stands for and what it predicts — Reczko et al., 2012; Paraskevopoulou et al., 2013
- [x] How candidate sites are identified — Reczko et al., 2012
- [x] How seed matching is used — Reczko et al., 2012
- [x] How conservation is used — Reczko et al., 2012
- [x] How free energy is used — Reczko et al., 2012 (RNAhybrid)
- [x] How site accessibility is modeled — Reczko et al., 2012 (Sfold, 3′ UTR only)
- [x] How CDS targeting differs from 3′ UTR targeting — Reczko et al., 2012
- [x] What machine-learning component is used — Reczko et al., 2012 (GLMs)
- [x] What its prediction score represents — Reczko et al., 2012
- [x] Its score range and cutoff — DIANA REST API (threshold 0.3–1; full range not documented)
- [x] Advantages/disadvantages compared with TargetScan — diana-microt-cds.md §8

---

## mirDIP

Stub: [mirdip.md](../../../docs/research/mirdip.md)

**Not a single algorithm like TargetScan** — integrative meta-resource.

- [x] Which prediction resources mirDIP integrates (verify current version) — 24 resources; mirDIP statistics; Hauschild et al., 2023
- [x] How individual predictions are combined — Tokar et al., 2018
- [x] How the mirDIP integrated score is calculated — Tokar et al., 2018 (noisy-or)
- [x] What the integrated score means — Hauschild et al., 2023; mirDIP portal
- [x] Its score range and confidence categories — Tokar et al., 2018; Hauschild et al., 2023
- [x] How mirDIP evaluates predictions against experimentally supported interactions — Tokar et al., 2018
- [x] Advantages of integrating multiple algorithms — Tokar et al., 2018
- [x] Disadvantages of relying on an integrated resource — mirdip.md §6
- [x] Whether combining methods creates correlated evidence or other biases — Tokar et al., 2018; mirdip.md §6

---

## miRDB

Stub: [mirdb.md](../../../docs/research/mirdb.md)

- [x] How the MirTarget machine-learning model works — Liu & Wang, 2019
- [x] What training data it uses — Liu & Wang, 2019
- [x] What features it considers — Liu & Wang, 2019 (96 features)
- [x] What its score means — Liu & Wang, 2019; miRDB FAQ
- [x] Its score range and cutoff — Chen & Wang, 2020; miRDB FAQ
- [x] Strengths and weaknesses as ML contrast to TargetScan/DIANA — mirdb.md §5–6

---

## Validation databases

Stub: [validation-databases.md](../../../docs/research/validation-databases.md)

**Not prediction algorithms.**

### miRTarBase

- [x] What it curates; current version — miRTarBase 2025 / v10.0 (Cui et al., 2025)
- [x] Evidence types (reporter assays, qPCR, CLIP, knockdown/overexpression, etc.) — Cui et al., 2025; Hsu et al., 2011
- [x] Support levels / confidence categories — Strong vs Weak (Cui et al., 2025)
- [x] How to use as reference overlay in our project — validation-databases.md integration section

### DIANA-TarBase

- [x] What it curates; current version — TarBase v9.0 (Skoufos et al., 2024)
- [x] Evidence types and curation approach — Skoufos et al., 2024
- [x] How to use as reference overlay in our project — validation-databases.md integration section

### Cross-cutting

- [x] Explicit prediction vs validation distinction documented — validation-databases.md
- [x] How miRTarBase and TarBase differ for validation use — validation-databases.md cross-cutting section

---

## Completion rule

A research doc may be marked `complete` only when:

1. All seven template sections are filled with source-backed content
2. Every applicable deep-dive item is checked OR marked "not documented" with source search noted
3. Section 4b score inventory rows are copied to [score-comparison-matrix.md](../../../docs/research/score-comparison-matrix.md)
4. Bibliography entries added to [bibliography.md](../../../docs/presentation/sources/bibliography.md)
