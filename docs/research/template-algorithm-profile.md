# Algorithm / Resource Profile Template

Copy this template for each prediction algorithm or validation database. Replace `{RESOURCE}` with the name (e.g., TargetScan, mirDIP).

**Status:** draft | in-review | complete

---

## 1. What is it?

### Who developed it?

<!-- TODO: names, institutions -->

### When was it introduced?

<!-- TODO: year, key publication -->

### Current / relevant version

<!-- TODO: version number, release date -->

### What problem does it solve?

<!-- TODO: 2–3 sentences -->

---

## 2. How does it work?

Describe the pipeline:

```
miRNA/mRNA input → candidate target sites → features → model → score → ranking/output
```

### Step-by-step process

<!-- TODO -->

### Candidate site identification

<!-- TODO -->

### Scoring and ranking

<!-- TODO -->

---

## 3. What biological/computational features does it use?

Verify each feature from primary sources. Mark ✅ only if documented.

| Feature | Used? | Source / notes |
|---------|-------|----------------|
| Seed matching | | |
| 8mer / 7mer-m8 / 7mer-A1 / 6mer site types | | |
| Evolutionary conservation | | |
| Binding / free energy | | |
| RNA / site accessibility | | |
| Local AU content | | |
| Target-site abundance | | |
| Position within transcript | | |
| 3′ compensatory pairing | | |
| G:U pairing | | |
| Machine learning | | |
| Other sequence/context features | | |

---

## 4. What does its score mean?

| Question | Answer |
|----------|--------|
| Score name(s) | |
| Numerical range | |
| Higher or lower = stronger? | |
| Score type (probability, confidence, rank, energy, etc.) | |
| Absolute interpretation across genes? | |
| Already normalized? | |
| Recommended cutoff / threshold | |
| Site scores aggregated to gene–miRNA score? | How? |

### Score interpretation notes

<!-- TODO: especially important for integration planning -->

### 4b. Score inventory

List every named score for this resource. Copy rows to [`score-comparison-matrix.md`](score-comparison-matrix.md) when verified.

| Score name | Level (site / gene–miRNA) | Evidence type | Range | Direction | Notes |
|------------|---------------------------|---------------|-------|-----------|-------|
| TBD | TBD | TBD | TBD | TBD | TBD |

Evidence types: `efficacy` | `conservation` | `ML confidence` | `meta-integration` | `experimental support`

---

## 5. Advantages

<!-- TODO: source-backed strengths only -->

---

## 6. Disadvantages

<!-- TODO: source-backed limitations only -->

---

## 7. When would you use it?

<!-- TODO: realistic research scenario -->

---

## 8. Method-specific deep dive

Complete the resource-specific checklist from [`deep-dive-checklists.md`](../../.cursor/skills/mirna-algorithm-research/deep-dive-checklists.md) or the stub file for `{RESOURCE}`. Mark each item `- [ ]` until verified; use `- [x]` only with a source citation.

<!-- Replace with PDF-mandated topics for this resource -->

- [ ] TBD — see deep-dive-checklists.md for {RESOURCE}

### Site-level to gene–miRNA aggregation (if applicable)

<!-- TODO: document how site scores roll up — especially critical for TargetScan -->

TBD

---

## Sources

<!-- Primary references only — one line per source -->

- 

---

## Integration notes (for consensus project)

<!-- How does this score compare to others? Can it be naively averaged? -->
