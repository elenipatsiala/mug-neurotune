# Preprint outline (bioRxiv, "Bioinformatics" or "Neuroscience" category)

**Working title:** Maturation or stress? A validated transcriptomic scorecard for human midbrain
organoids, applied to dopaminergic organoids cultured on the International Space Station

**Target venues:** bioRxiv preprint → *npj Microgravity* (Brief Communication) or
*Journal of Emerging Investigators* (if submitting as a high-school author with a mentor).

## Abstract (fill numbers after the real run; ≤ 250 words)
Problem (organoid immaturity limits PD modelling; ISS organoids reported "faster maturation") →
gap (maturation-like transcription can come from stress/hypoxia; no validated way to tell) →
what we built (scorecard validated on known-age organoids, anchored to fetal midbrain) →
what we found on GSE259421 (classification + hypoxia-adjusted result + fetal-age estimate) →
O₂ model result → design implication (donors needed; O₂ must be controlled) → availability.

## Figures
1. **Concept + validation.** (a) Maturation vs fidelity vs stress axes; (b) Step A on Fiorenzano
   time-course; (c) Step B fetal-age estimates with published 11/15 PCW cross-check.
2. **ISS result.** Maturation–Fidelity map (LEO vs ground, fill = hypoxia); effect table with
   CIs/q; hypoxia-adjusted bar (raw vs adjusted g).
3. **Oxygen physics.** Vial schematic; paired Monte Carlo (fig_o2_montecarlo); regime map;
   conditional analysis given normoxic ground organoids.
4. **Next mission.** Power grid; recommended design (donors × organoids) + O₂-controlled vessel.
Supplementary: gene sets + sources, coverage, leave-one-gene-out, per-line results, z-score
sensitivity, classifier calibration, cortical-organoid secondary check.

## Methods sections
Data sources · scoring (background-referenced ranks) · validation criteria · anchoring ·
stratified permutation/bootstrap · classification rules · hypoxia adjustment · power simulation ·
O₂ model and verification · pre-registration (OSF link) · code and data availability (GitHub + Zenodo DOI).

## Author checklist before posting
- [ ] Pre-registration posted on OSF *before* real-data run (timestamp)
- [ ] Real-data results reproduced from a clean clone with one command
- [ ] At least 2 domain experts read it (organoid biologist + space biologist); feedback logged
- [ ] Contacted Marotta / Fossati / Loring with results before posting (courtesy + error check)
- [ ] Every claim mapped to a figure; non-claims section present
- [ ] Mentor/affiliation listed; bioRxiv requires an institutional or independent affiliation
