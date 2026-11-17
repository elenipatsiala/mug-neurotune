# µG-NeuroTune

**Is your midbrain organoid really maturing, or just stressed?**

µG-NeuroTune is an open, validated scorecard that profiles human midbrain (dopaminergic)
organoids on three axes from RNA-seq: **Maturation**, **Dopaminergic fidelity** and
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
pytest -q                                    # 12 ground-truth tests (~30 s)
python o2model/run_o2.py results/o2 3000     # O2 Monte Carlo + figures
```
Real data: follow **`docs/DATA_GUIDE.md`** (exact files and commands).

## Current status (v0.1.0, 2026-11-17)
* ✅ Pipeline built and tested on synthetic data with planted effects: it recovers true
  maturation, rejects null data, flags hypoxia-driven "maturation" as *mixed state*, and
  recovers a planted fetal age (12 PCW → 12.3–12.5).
* ✅ Classifier calibrated on **800 simulated CRS-19-sized missions** (200 per truth):
  no effect → correct 99.5% · true maturation → 95.5% · hypoxia-driven → 92.5% ·
  identity loss → 100% (`results/calibration.json`).
* ✅ O₂ model run (3,000 paired simulations; solver verified analytically). See `docs/O2_MODEL.md`.
* ⏳ **Real-data runs pending**: ISS (GSE259421), Fiorenzano time-course, fetal reference.
  No biological conclusion is claimed until those runs are done.

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
* **v0.2:** real-data results (A–E) + preprint.
* **Phase 2 Create:** O₂-controlled culture for microgravity. Use the O₂ model to size a perfused
  or gas-permeable vessel that keeps the organoid normoxic in µg, plus a ground test plan.

## Cite
See `CITATION.cff`. Primary data: Marotta D. et al. (2024) *Stem Cells Transl Med* 13(12):1186–1197,
doi:10.1093/stcltm/szae070.

License: MIT.
