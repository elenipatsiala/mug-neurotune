# Clock v0.3 frozen before the held-out test (2026-11-17)

`mug_neurotune/data/neuronal_age_clock_v03.json` is frozen in this commit, **before** it is applied to
the held-out organoid dataset (GSE277032 midbrain-like organoids). The test and its pass rule are fixed now:

* Held-out data: GSE277032 organoid single-cell samples (time points from the series matrix).
* Pseudobulk per sample (whole organoid; neurons only with the same gate as before).
* **Pass:** Spearman ρ(clock, culture day) ≥ 0.6 and bootstrap 95% CI lower bound > 0, on whole organoids.
* The result is reported whatever it is, and the clock is not changed in response.

## Why v0.3 exists (history, for transparency)
* v0.2 (trained on one fetal dataset) passed the first organoid test (Lund, ρ = 0.91) but **failed** the
  second (BrainSTEM GSE281535, ρ = 0.49). Its "up" genes did not replicate in an independent fetal
  dataset (mean ρ = −0.01 in La Manno 2016), so they were dataset-specific.
* v0.3 rule (decided before scoring any organoids with it): keep only genes whose age trend has the
  **same sign in both** fetal datasets (Birtele 2022, La Manno 2016); rank by the weaker of the two.
* v0.3 was then scored on both development organoid sets: Lund ρ = 0.74 (whole) / 0.95 (neurons);
  BrainSTEM ρ = 0.99 / 0.99. Because v0.3 was designed after seeing v0.2 fail, these are development
  results, not proof; the held-out test above is the proof.

## Outcome of the held-out test (run 2026-11-17, after this file and the clock were pushed)
GSE277032 organoid pools, day 40 / 70 / 120 (3 pools each; day-40/70 pools also contain some 2D cultures).
**Clock v0.3, whole organoid: ρ = 0.95 (95% CI 0.77–0.97) → PASS.** Neurons only: ρ = 0.95.
Mean clock by day: 40 → −0.023, 70 → −0.012, 120 → −0.007 (monotonic), although neurons fall from
57% to 23% of cells by day 120. The clock was not changed after this test.
