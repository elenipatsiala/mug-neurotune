# µG-NeuroTune

**Is your midbrain organoid really maturing, or just stressed?**

µG-NeuroTune is an open, validated scorecard that profiles human midbrain (dopaminergic)
organoids on three axes from RNA-seq: **Neuronal age (validated clock)**, **Dopaminergic fidelity** and
**Stress/Hypoxia**. It ships with a spaceflight module, built around the NASA ISS organoid
experiment (Marotta et al. 2024; GEO GSE259421), and an **oxygen-transport model** of the
flight culture vials.

> Why: organoids are too immature to model Parkinson's disease well, and "maturation-like"
> gene expression can also come from cellular stress or low oxygen (Bhaduri et al. 2020;
> Vértesy et al. 2022). The ISS study reported faster maturation in orbit. This tool asks
> whether that signal is *faithful* maturation, and what a future mission must control.

## What's inside
| Step | What it does | Module / command |
|---|---|---|
| A · validate | Checks the maturation score rises with known culture age | `mugnt validate` |
| B · anchor | Estimates the fetal age (PCW) each organoid resembles | `mugnt anchor` |
| C · compare | LEO vs ground within cell line × microglia strata; permutation p, bootstrap CI, FDR; pre-registered verdict | `mugnt compare` |
| D · hypoxia | Hypoxia-adjusted maturation test + O₂ physics model of the vial | `mugnt compare`, `o2model/` |
| E · design | Power for the next mission (donors × organoids) | `mugnt power` |
| F · app | Browser scorecard: upload counts, get the profile (data stays local) | `webapp/index.html` |

## Quick start
```bash
pip install -e .            # numpy, scipy, pandas, matplotlib
mugnt demo --out results/demo_synthetic      # end-to-end on SYNTHETIC data with known truth
pytest -q                                    # 17 tests (~30 s)
python o2model/run_o2.py results/o2 3000     # O2 Monte Carlo + figures
```
Real data: follow **`docs/DATA_GUIDE.md`** (exact files and commands).

## Current status (v0.2.0, 2026-11-17): run on real data
* ✅ **Validated maturation axis.** A neuronal age clock learned from human fetal midbrain neurons
  (GSE192405) orders an *independent* organoid time-course (GSE168323, days 15–120) with
  **ρ = 0.91** (95% CI 0.71–0.96). The v0.1 marker set failed the same test (ρ = −0.16) and is now
  descriptive only.
* ✅ **Applied to NASA ISS organoids (GSE259421).** No change on the validated clock (g = 0.00,
  q = 0.995). Glycolysis up (g = +1.41, q = 0.017) and p53/DNA-damage down (g = −3.76, q = 0.011),
  same direction in both cell lines and in cortical organoids from two other donors.
* ✅ **Mission design:** 2 donors gives ≤ 28% power for a 2-week maturation effect; 4 donors × 4
  organoids per arm gives 88%.
* ✅ Classifier calibrated on 800 simulated missions (92.5–100% correct). Browser and Python give
  identical scores on the real data. 17 tests pass.
* Full write-up: **`docs/RESULTS.md`** · figure: `results/real/fig_main_results.png`.

## Design choices that make it trustworthy
* **Pre-registered** gene sets and decision rules (`docs/PREREGISTRATION.md`), fingerprinted
  in every output.
* **Non-overlapping programs** and **background-referenced ranks**, so one program changing
  cannot fake changes in another (a test enforces this exactly).
* **Within-stratum statistics**, so a donor difference can never look like a LEO effect.
* **Every robustness check reported**, none selected: z-score method, leave-one-gene-out,
  per line, per microglia status.
* **Explicit non-claims:** transcripts ≠ function; LEO ≠ microgravity only; 2 donors ≠ general.

## Repository map
```
mug_neurotune/   genesets · io · scoring · validate · anchor · compare · power · plots · cli · synthetic
o2model/         vial_o2.py (model) · run_o2.py (Monte Carlo, figures)
webapp/          index.html (built) · template.html · core.js (scoring, identical to Python)
tests/           ground-truth tests
examples/        calibration.py
docs/            DATA_GUIDE · PREREGISTRATION · O2_MODEL · PREPRINT_OUTLINE
results/         demo_synthetic/ (SYNTHETIC) · o2/ · calibration.json
```

## Roadmap
* **v0.3:** more fetal ages (2nd trimester) for the clock, day-90 organoids, expert review of gene sets, preprint.
* **Phase 2 Create:** O₂-controlled culture for microgravity. Use the O₂ model to size a perfused
  or gas-permeable vessel that keeps the organoid normoxic in µg, plus a ground test plan.

## Cite
See `CITATION.cff`. Primary data: Marotta D. et al. (2024) *Stem Cells Transl Med* 13(12):1186–1197,
doi:10.1093/stcltm/szae070.

License: MIT.
