# Prediction vs measured evidence

This table checks whether the study's computational predictions agree with
independently published experimental measurements. It does not replace those
experiments; it shows the docking and network predictions are not isolated
computational artifacts.

## Concordance summary

| Item | Our prediction | Published measurement | Assay |
|---|---|---|---|
| Luteolin-PTGS2 | DUD-E score -9.63; SwissTarget 0.941 | COX-2 IC50 ~0.65 ug/mL | cell-free enzyme |
| Luteolin-COX2 protein | predicted engagement | COX-2 protein suppressed dose-dependently | RAW264.7 + LPS |
| alpha-Cyperone-COX2/PGE2 | docking predicts COX-2 | inhibits COX-2 and PGE2 via NF-kB; lowers IL-6 | RAW264.7 + LPS |
| alpha-Cyperone-pain | predicted nociceptive role | reduces paclitaxel neuropathic pain | animal pain model |
| IL1B->PTGS2 | CellOracle/network direction | established experimental route | Samad 2001, Piston 2007 |
| Luteolin L1000 | Down signature includes PTGER3/PTGS2 | independent perturbation signature | Enrichr |

## Boundaries (must be stated)

- Every measured result concerns a **single compound in an isolated system**, not
  the whole seven-herb formula and not a liver-disease dose context.
- The alpha-cyperone pain result is animal evidence and is external to this
  computational study; it is used only to show the predicted nociceptive role is
  biologically plausible.
- Agreement with measured data supports the reliability of the **protocol and
  ranking**, not direct binding of the formula, which still requires a dedicated
  assay (cell-free COX-2 kit, CETSA/DARTS, or a cell model).

## Sources

- Luteolin COX-2 IC50: https://www.tandfonline.com/doi/full/10.1080/13880209.2023.2200841
- Luteolin RAW264.7: https://pmc.ncbi.nlm.nih.gov/articles/PMC3865263/
- alpha-Cyperone: https://pmc.ncbi.nlm.nih.gov/articles/PMC11679560/
- Enrichr: https://maayanlab.cloud/Enrichr/
