# Pre-registration: µG-NeuroTune analysis of ISS dopaminergic organoids

Frozen on 2026-11-17, **before** any real ISS, time-course or fetal data were scored.
Gene sets v0.1.0, fingerprint printed by `python -c "from mug_neurotune import genesets; print(genesets.fingerprint())"`.
Post this file to OSF (https://osf.io) with a timestamp before running the real analysis.

## Questions
1. **A (instrument validity).** Does the Maturation score rise with days in culture in
   Earth-grown midbrain organoids of known age?
2. **B (developmental anchoring).** What fetal age (PCW) do ground and LEO organoids most resemble?
3. **C (main question).** In LEO vs ground dopaminergic organoids, does maturation rise while
   dopaminergic identity/fidelity is preserved and stress/hypoxia programs are not raised?
4. **D (oxygen hypothesis).** Is any LEO maturation effect partly explained by a hypoxia
   program (hypoxia-adjusted test)?
5. **E (design).** How many donors and organoids would a future mission need?

## Analysis plan (fixed)
* **Step 0:** check gene symbols against HGNC; report coverage; programs with < 60% genes found
  are flagged unreliable. Any gene-set edit after this point bumps the version and is logged.
* **Scores:** background-referenced rank score (primary); mean z-score (sensitivity).
* **A pass criterion:** Spearman ρ(maturation, day) ≥ 0.6, bootstrap 95% CI lower bound > 0,
  and progenitor/cycle ρ < 0. **If A fails, the maturation axis is not used for any claim in C.**
* **B:** age-informative genes = top 300 by |ρ with PCW| in the fetal reference;
  softmax-weighted age estimate. Cross-check: Fiorenzano day 40–70 organoids should land at
  roughly 11–15 PCW. Report reference similarity; < 0.5 = "does not resemble fetal midbrain".
* **C:** DA organoids only. Strata = cell line × microglia. Effect = size-weighted
  within-stratum Hedges' g; p from 5,000 within-stratum permutations; CI from 2,000 stratified
  bootstraps of the raw difference; Benjamini–Hochberg across programs.
  A program is "changed" only if |g| ≥ 0.5 AND CI excludes 0 AND q < 0.05 AND the direction agrees in both cell lines.
* **Classification (exact rules):** *insufficient evidence* if maturation is not up;
  *dysregulated* if maturation up and DA identity or fidelity down; *mixed* if maturation up and
  any dysregulation program (glycolysis, ER stress, hypoxia, p53, reactive glia) up;
  otherwise *promising maturation*.
* **D:** within-stratum regression of maturation on hypoxia; share of the LEO effect explained
  = 1 − g_adjusted / g_raw. > 0.5 → "hypoxia explains most of the effect".
  Interpreted together with the O₂ transport model (o2model/), never alone.
* **Robustness (all reported, none selected):** z-score method; leave-one-gene-out for
  maturation; each line separately; ± microglia separately; same analysis on cortical organoids
  as a secondary, non-confirmatory check.
* **E:** variance components from the DA data; simulation over 2–12 donors × 2–12 organoids
  per donor per arm, donor×LEO SD = 0, 0.5, 1 × organoid SD; target power 0.8.

## Things this analysis cannot show (non-claims)
* that microgravity (rather than radiation, launch, return at 4 °C, or the CubeLab) caused a change;
* functional maturity (electrophysiology, dopamine release): transcripts only;
* that LEO organoids are better Parkinson's disease models;
* anything that generalises beyond 2 donors: with 2 lines, donor-level inference has 1 df;
* actual oxygen levels in the flight vials (not measured; the model gives plausibility only).

## Calibration of the classifier
`examples/calibration.py` simulates CRS-19-sized missions with known truth (no effect /
true maturation / hypoxia-driven maturation) and reports how often each call is right.
Results (`results/calibration.json`, 200 simulated missions per truth, q < 0.05 rule):
| Truth | Correct call | Accuracy (95% CI) |
|---|---|---|
| no effect | insufficient evidence | 99.5% (98.5–100) |
| true maturation | promising maturation | 95.5% (92.6–98.4) |
| hypoxia-driven maturation | mixed state | 92.5% (88.8–96.2) |
| maturation + identity loss | dysregulated state | 100% |
The q < 0.05 cut-off was chosen over q < 0.10 by a rule fixed in advance (best worst-case
accuracy: 92.5% vs 92.0%), on simulation only.
