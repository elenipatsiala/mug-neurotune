# Changelog (all changes made before any real data were scored)

## v0.1.0 — 2026-11-17
* Scoring changed from whole-transcriptome ranks to **background-referenced ranks** after a
  synthetic test showed rank leakage (a rise in one program created fake decreases in others).
  A deterministic test now proves untouched programs are unaffected exactly.
* Bootstrap CIs moved from Hedges' g to the raw score difference (g explodes when 3-organoid
  cells resample to zero variance).
* Added FDR control and a **cross-line consistency** requirement to the classification.
* FDR cut-off set to **q < 0.05** by the pre-specified calibration rule (best worst-case accuracy
  over four simulated truths, 200 missions each).
