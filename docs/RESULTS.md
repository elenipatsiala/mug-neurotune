# Results on real data (v0.2.0, 2026-11-17)

Main figure: `results/real/fig_main_results.png` (script: `examples/figure_real_results.py`).

## Data
| Dataset | What | Used for |
|---|---|---|
| GSE259421 (Marotta et al. 2024) | Bulk RNA-seq, one library per organoid. NCBI counts contain **18 of 23 dopaminergic** organoids (UEC741 control, HDF410 PD; ± microglia) and **18 of 19 cortical** organoids (2 donors). GEO labels the flight as "38 days microgravity exposure". | Steps C, D, E |
| GSE168323 (Fiorenzano et al. 2021) | scRNA-seq of midbrain organoids. **16 organoids** intact: day 15 (2), 30 (5), 60 (5), 120 (4). Day 90 and one day-120 file were incomplete in the download. | Step A (validation) |
| GSE192405 (Birtele et al. 2022) | scRNA-seq of human fetal ventral midbrain tissue at 6, 8, 10.5, 11.5 PCW. | Training the age clock; Step B |

## Step A: does the maturation score measure maturation? (pre-registered test)
* **v0.1 marker set: FAILED.** Spearman ρ with culture day = −0.16 (CI −0.66 to 0.52) on whole organoids;
  0.41 (CI −0.19 to 0.81) on neurons only. Cause 1: **cell composition**. Neurons fall from 35%
  of cells (day 60) to 4% (day 120) while glia rise to 48%, so a bulk marker score tracks neuron
  *content*. Cause 2: **ceiling**. Within neurons the synaptic genes are already near-maximal at day 15.
  Per the pre-registration, the v0.1 score is **not used for any maturation claim**.
* **v0.2 neuronal age clock: PASSED.** Learned on fetal neurons only (50 genes up with age, dominated by
  mitochondrial respiration; 50 down, dominated by immature-neuron cytoskeleton), frozen, then tested
  blind on the organoids: **ρ = 0.91 (CI 0.71–0.96) whole organoid; ρ = 0.84 (CI 0.60–0.93) neurons only.**
  Monotonic across days 15 → 30 → 60 → 120.
* Honest caveats: the clock was designed *after* v0.1 failed. It was trained only on fetal data, which
  are independent of the organoids, so the organoid test is out-of-sample. Each fetal age is one
  donor, so age and batch are confounded in training; dissociation-stress and glial genes were removed
  by fixed published lists, not by tuning.

## Step B: which fetal age do the organoids resemble? (approximate)
Linear calibration on 4 fetal tissue samples (r = 0.78; valid only within 6–11.5 PCW).
Organoids: day 15 ≈ 9.7, day 60 ≈ 11.2, day 120 ≈ 12.2 PCW (extrapolated). The published atlas places
day 70–120 organoids near 15 PCW; our reference stops at 11.5 PCW and cannot reach that far.
**ISS and ground dopaminergic organoids both map to ≈ 10.8 PCW, with no difference.**

## Step C: ISS vs ground, dopaminergic organoids (n = 18; strata = line × microglia)
| Program | Hedges' g | q (BH) | Same direction in both lines | Call |
|---|---|---|---|---|
| **Age clock (validated)** | −0.00 | 0.995 | yes | **no change** |
| Synaptic markers (v0.1, descriptive) | +0.83 | 0.143 | yes | n.s. |
| DA identity | −0.26 | 0.816 | no | n.s. |
| Progenitor / cell cycle | −1.04 | 0.143 | yes | n.s. |
| **Glycolysis** | **+1.41** | **0.017** | yes | **up** |
| Hypoxia (HIF targets) | +0.91 | 0.143 | yes | n.s. (trend) |
| **p53 / DNA damage** | **−3.76** | **0.011** | yes | **down** |

**Pre-registered verdict: "insufficient evidence" for faster maturation.**
* All 9 glycolysis genes rise (LDHA ≈ 2.1×, HK2 ≈ 1.8×). Several HIF targets rise too (PDK1, BNIP3,
  NDRG1, P4HA1), others do not (VEGFA, DDIT4).
* p53 targets fall (GADD45A ≈ 0.54×, CDKN1A ≈ 0.84×), consistent with the original paper.
* **Cortical organoids (2 other donors, secondary check)** show the same direction: glycolysis
  g = +1.80 (q = 0.11), hypoxia +1.23, p53 −3.31 (q = 0.017), and no clock change (g = −0.49, n.s.).
* Robustness: z-score method gives the same signs for 82% of programs (glycolysis g = +1.55, q = 0.10).

## Step D: is the "maturation" signal a metabolic signal?
* Within strata, the v0.1 synaptic-marker score correlates with glycolysis (r = 0.68). Adjusting for
  glycolysis removes 79% of its apparent ISS effect. **Exploratory**: the pre-registered covariate was the
  hypoxia program, which removes 24%.
* On the validated clock there is no effect left to explain.
* Interpretation (hypothesis, not proven): ISS culture shifted organoid metabolism toward glycolysis,
  which is the opposite of the glycolysis-to-respiration switch of maturing neurons. Marker-based
  "maturation" read-outs pick up this shift. Fits the O₂ model (`docs/O2_MODEL.md`): both culture
  conditions sit near a diffusion limit, and µg can tip the organoid into O₂ limitation.

## Step E: what should the next mission fly? (variance from these data)
* Organoid-to-organoid SD on the clock ≈ 6.6 culture-days of maturation.
* Power to detect **2 weeks** of faster maturation: **2 donors (CRS-19): ≤ 28%** even with 12 organoids
  per donor per arm; **4 donors × 4 organoids: 88%; 6 donors × 2: 96%** (donor × LEO SD = 0.5 × organoid SD).
* Confirming the glycolysis shift needs about **4 donors × 2 organoids per arm** (power ≈ 0.81–0.88).
* Rule: **add donors before adding organoids.**

## What this does and does not show
Shows: a validated, reusable maturation clock for midbrain organoids. On it, ISS dopaminergic organoids
were not detectably more mature after the mission, while glycolysis rose and p53 targets fell
reproducibly in both organoid types.
Does not show: that microgravity (rather than launch, radiation, return transport or culture hardware)
caused the shift; any functional property; anything general from 2 donors per organoid type; actual O₂
levels in the flight vials. The CI on the clock excludes large effects, but not small ones (< ~1 week).

## Process notes (transparency)
* Gene sets and decision rules were frozen in git (commit 076fcad, 2026-11-17) **before** any real data
  were scored, but were **not** yet posted on OSF; do that now, citing this history.
* v0.2 changes after seeing real data: the clock was added (because v0.1 failed Step A), and the
  "share explained" metric is now not reported when there is no raw effect. Logged in `docs/CHANGELOG.md`.
