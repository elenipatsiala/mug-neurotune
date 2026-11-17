# Oxygen-transport model of the CRS-19 organoid vials (hypothesis D)

**Question.** In microgravity there is no buoyancy-driven convection and organoids do not
sediment. Could the ISS organoids have been more O₂-limited than ground controls *in identical
sealed vials*, so that part of the "faster maturation" signal is really a hypoxia response?

**Code:** `o2model/vial_o2.py` (model), `o2model/run_o2.py` (Monte Carlo + figures).
**Outputs:** `results/o2/` (`o2_summary.json`, `o2_conditional.json`, `fig_o2_montecarlo.png`,
`fig_o2_regime_map.png`, `sensitivity_delta_hypoxic.csv`).

## Facts from the paper vs assumptions
| Item | Value | Status |
|---|---|---|
| Vessel | Nunc cryovial #374088, 1 organoid/vial, 1 mL medium, sealed, 37 °C, 30 d, no medium change | paper |
| Organoid size | 200–500 µm diameter, formed from 1×10⁵ cells | paper |
| Hypoxia | "no signs of hypoxia in either group" (method not stated) | paper |
| Headspace volume, wall O₂ permeability | not reported | **assumed** (sealed, impermeable = worst case) |
| O₂ use per cell | 5–100 amol/cell/s (log-uniform) | assumed range |
| Convection in ground vial | 0.5–50 µm/s | assumed range |
| Residual mixing on ISS | 0.001–0.5 µm/s | assumed range |

## Model (three levels)
* **L0 budget.** Headspace + dissolved O₂ vs 30-day demand.
* **L1 vial transport.** Source (gas) → organoid. 1 g: organoid at the bottom, gas at the top,
  thermal convection adds mixing. µg: bubble and organoid at random positions, no convection.
* **L2 organoid.** External boundary layer (Sherwood number; 2 for pure diffusion) in series
  with internal spherical diffusion + Michaelis–Menten uptake, solved by Newton finite volumes.
  **Verified** against the exact zero-order solution (max error 3×10⁻⁶) and mass balance (exact).

## Results (3,000 paired Monte Carlo runs, same parameters in both gravities)
* Most organoid-scale transport resistance sits in the **boundary layer around the organoid**
  (vial-scale resistance is only ~16–20% of it). At the organoid's scale, gentle convection
  barely helps (Péclet number ≲ 6), so **the two conditions are mostly similar**:
  * similar within 5 percentage points of hypoxic volume: **65%** of runs
  * ground more hypoxic: **22%** · ISS more hypoxic: **13%**
* What decides the direction (Spearman with the µg − 1 g difference): **convection speed in the
  ground vial** (ρ = 0.57), organoid–bubble distance on ISS (0.27), and how much of the
  organoid surface touches the vial wall on the ground (0.24).
* **Tension with the paper.** Across the full ranges, most runs predict substantial hypoxia
  in *both* groups. Only ~6% of runs match the reported normoxic ground organoids, and those
  need low O₂ use (median ~6 amol/cell/s).
  **Conditional on normoxic ground organoids**, the ISS organoid becomes partly hypoxic in
  ~15% of runs (up to 37% of its volume). Near the threshold it is close to a coin flip
  (41% ISS worse vs 36% ground worse).
* A sealed, gas-impermeable vial would run out of O₂ before day 30 in ~44% of runs. Whether the
  vial wall lets O₂ through is the biggest unknown and should be asked of the authors / Space Tango.

## What this means
1. Physics says hypothesis D is **plausible but not guaranteed**. It is not a reason to assume
   the ISS maturation signal is hypoxia. It is a reason to **test** it with the hypoxia program
   and the hypoxia-adjusted analysis (`mugnt compare`).
2. Both culture conditions sit **near a diffusion limit**. So oxygen was an uncontrolled
   variable in the experiment, and that alone justifies the Phase-2 Create: a vessel that
   controls O₂ delivery in µg (perfused or gas-permeable culture), sized with this model.

## Honest limitations
Steady state (diffusion time over the vial ≈ 15 h ≪ 30 d, so reasonable); spherical organoid;
uniform cell density; no necrotic-core feedback on uptake; convection represented by
mixing-length and creeping-flow correlations rather than CFD; priors are wide ranges, not
measurements. All of this is stated in every figure caption and must be in any write-up.
