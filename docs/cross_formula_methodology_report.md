# Cross-formula methodology and convergence report

This report documents the second test case added to turn the single-formula
network-pharmacology study into a reusable, de-circularized methodological
framework. The disease axis (liver inflammation to pain sensitization) is held
fixed and only the formula input changes, so the case tests whether the pipeline
transfers rather than whether a new disease is explained.

## 1. Why a second formula

The first case uses Chaihu Shugan San, with vinegar-processed Cyperus as its
qi-regulating component. A reviewer's concern was that a framework demonstrated
on one formula and one disease axis could be a tuned, single-use analysis. The
second case applies the identical acquisition, filtering, mapping, and
benchmarking steps to a different formula indicated for the same liver-related
pain phenotype.

## 2. Formula selection

Two candidates were processed.

### Jinlingzi San (selected second case)

Jinlingzi San is a two-herb formula from *Suwen Bingji Qiyi Baoming Ji*,
composed of Toosendan Fructus and Corydalis Rhizoma, used for liver qi
transforming into fire and for chest and hypochondriac pain. It shares no herbs
with Chaihu Shugan San, which is the strongest possible input-level separation.

### Xiaoyao San (exploratory control, not the formal case)

Xiaoyao San was processed first but shares the monarch herb Chaihu and several
other herbs with the first formula. Its filtered gene set overlaps the first
formula at Jaccard 0.932 (246 shared genes), so it cannot serve as a separate
case. It is retained as an exploratory control that illustrates how shared
herbs drive shared targets.

## 3. Acquisition and processing

The two Jinlingzi herbs were fetched from TCMSP using the standard search form
and the English registered names "Toosendan Fructus" and "Corydalis Rhizoma".
The embedded kendoGrid JSON arrays were parsed and the same filter applied to
both formulas: oral bioavailability at least 30 and drug-likeness at least 0.18.

Counts:

| Quantity | Chaihu Shugan San | Jinlingzi San |
|---|---|---|
| Herbs | 7 | 2 |
| Unique compounds | 938 | 109 |
| Compounds passing filter | 139 | 57 |
| Filtered human genes | 254 | 204 |
| DEG intersection | 51 | 42 |
| External pain overlap | 6 | 5 |

Even with zero herb overlap, the filtered gene sets overlap at Jaccard 0.703
(189 shared genes). The all-mapped, unfiltered background gives a lower Jaccard
of 0.359. This is an honest finding: the TCMSP human target network is densely
connected, so distinct formulas converge at the gene level. The second case
therefore cannot be described as an independent gene set.

## 4. Symmetric de-circularized benchmark

The disease DEG rule and the external pain reference are identical for both
formulas. For each formula, the selected set is the DEG intersection with its
own universe and the positive labels are the literature-curated pain genes
detectable in that universe. A 10,000-draw permutation within each universe
gives the null, with a one-sided Fisher exact test as a closed-form check.

| Metric | Chaihu Shugan San | Jinlingzi San |
|---|---|---|
| Universe size | 254 | 204 |
| External pain genes in universe | 9 | 8 |
| Selected (DEG intersect) | 51 | 42 |
| Observed pain overlap | 6 | 5 |
| Observed precision | 0.1176 | 0.1190 |
| Chance expected overlap | 1.81 | 1.65 |
| Enrichment vs chance | 3.32 | 3.04 |
| Null precision mean | 0.0353 | 0.0388 |
| Null precision 95th percentile | 0.0784 | 0.0952 |
| Fisher p | 0.002593 | 0.010265 |
| Empirical p | 0.0020 | 0.0092 |

The two formulas reach almost identical precision, about 0.118, each roughly
three times the permutation expectation, and both are significant. The two
precision values should not be counted as two independent confirmations, because
the universes overlap and the disease and reference inputs are shared.

## 5. Core convergence and the PTGES difference

The strict six-gene core is PTGS2, IL6, IL1B, IL1A, PTGES, and PTGER3.
Jinlingzi San recovers five of the six and lacks PTGES, because none of its
compounds hits prostaglandin E synthase. The shared five are IL1A, IL1B, IL6,
PTGER3, and PTGS2. The imperfect match is more credible than a perfect copy:
PTGES is a genuine first-formula-specific node rather than a forced convergence.

## 6. Engineering decision on downstream layers

The single-cell and spatial layers were not rerun for the overlapping genes,
because those analyses operate on the same genes and would produce the same
outputs. The second case adds an upstream portability and convergence test, not
new downstream evidence. This is a deliberate, recorded engineering choice, and
the limitation is stated in the manuscript rather than hidden.

## 7. Reproducibility

The benchmark is fully reproducible from committed files. The self-contained
script `scripts/13_cross_formula/01_cross_formula_benchmark.py` reads only the
small data in `data/processed_small` and reproduces every metric in this report
without network access. The fetch and parse scripts in the same folder document
the end-to-end acquisition and mapping steps.
