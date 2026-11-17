# Benchmark: which maturity read-outs can be trusted in midbrain organoids?

Figure: `results/real/fig_benchmark.png` · script: `examples/figure_benchmark.py` · tables: `results/real/benchmark/`.

Three midbrain-organoid single-cell time courses from three labs, scored as whole-organoid pseudobulk:
Lund (GSE168323, day 15–120, 16 organoids), BrainSTEM (GSE281535, day 20–60, 48 samples), and a
held-out set (GSE277032, day 40/70/120, 9 pools; tested once after freezing the clock).

| Read-out | Lund | BrainSTEM | Held-out | All 3 ≥ 0.6 |
|---|---|---|---|---|
| **Age clock v0.3** | 0.74 | 0.99 | 0.95 | ✅ |
| **Fetal-correlation age** | 0.93 | 0.84 | 0.84 | ✅ |
| **Progenitor loss (SOX2+NES)** | 0.94 | 0.98 | 0.95 | ✅ |
| SOX2 ↓ alone | 0.87 | 0.95 | 0.95 | ✅ |
| NES ↓ alone | 0.77 | 0.98 | 0.95 | ✅ |
| MKI67 ↓ | 0.73 | 0.94 | 0.26 | ❌ |
| Age clock v0.2 (retired) | 0.91 | 0.49 | 0.95 | ❌ |
| Synaptic set v0.1 (retired) | −0.16 | 0.97 | −0.58 | ❌ |
| NeuN/RBFOX3 ↑ | −0.31 | 0.92 | −0.47 | ❌ |
| SYP ↑ | −0.02 | 0.90 | −0.47 | ❌ |
| TH ↑ | −0.26 | 0.88 | −0.74 | ❌ |
| MAP2 ↑ | 0.12 | 0.16 | −0.95 | ❌ |

## Take-aways
1. **"Maturity" markers that papers commonly show (NeuN, SYP, TH, MAP2) do not track organoid age** whenever
   cell composition changes. They correlate with the fraction of neurons (ρ 0.75–0.98 in Lund and held-out).
   They only "work" in BrainSTEM, where neuron content stays high.
2. **Age is reliably visible as the loss of progenitor/immature genes** (SOX2, NES) and in fetal-referenced
   read-outs (clock v0.3, fetal correlation).
3. **The age clock is not more accurate than a SOX2+NES score.** Its value is a fetal-anchored, multi-gene
   read-out that is independent of the progenitor genes (SOX2 and NES are excluded from the clock), so the
   two give independent evidence. That is why v0.3 calls a maturation change only when the panel agrees.

## Caveats
* The progenitor-loss read-out was chosen after seeing the Lund and BrainSTEM benchmark (it was also
  validated on the held-out set, but the choice was not pre-registered). The clock v0.3 was frozen first.
* Three datasets, all midbrain; cortical organoid time courses are not yet tested.
* Composition is estimated with a simple marker gate (neuron vs glia).
