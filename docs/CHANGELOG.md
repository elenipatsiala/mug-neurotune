# Changelog

## v0.2.0 — 2026-11-17 (after first real-data run)
* **Step A failed for the v0.1 maturation marker set** on real organoids (ρ = −0.16). Kept as a
  descriptive "synaptic markers" score; never drives the verdict.
* **Added the neuronal age clock** (`mug_neurotune/clock.py`, trained by
  `examples/train_fetal_clock.py` on fetal neurons only, frozen before the organoid test), now the
  maturation axis. Validated out-of-sample: ρ = 0.91.
* Gene-name aliases for renamed HGNC symbols (ATP5MPL→ATP5MJ, ATP5MD→ATP5MK, H3F3A→H3-3A).
* CLI: uses exactly the samples in the metadata file and reports the rest; `pseudobulk` command;
  clear errors for Ensembl/NCBI IDs without an annotation file.
* Fixed: clock silently dropped from `score_all` (argument overwritten); z-score sensitivity
  comparison crashed when programs differed; "share explained" no longer reported without a raw effect.
* Web app: validated clock, real ISS demo data, identical to Python on real data (max diff 7e-16).

# Pre-data changes (v0.1.0)

## v0.1.0 — 2026-11-17
* Scoring changed from whole-transcriptome ranks to **background-referenced ranks** after a
  synthetic test showed rank leakage (a rise in one program created fake decreases in others).
  A deterministic test now proves untouched programs are unaffected exactly.
* Bootstrap CIs moved from Hedges' g to the raw score difference (g explodes when 3-organoid
  cells resample to zero variance).
* Added FDR control and a **cross-line consistency** requirement to the classification.
* FDR cut-off set to **q < 0.05** by the pre-specified calibration rule (best worst-case accuracy
  over four simulated truths, 200 missions each).
