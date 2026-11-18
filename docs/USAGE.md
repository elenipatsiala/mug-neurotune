# Run µG-NeuroTune on your own organoid data (≈ 10 minutes)

Written for teams with new spaceflight or ground organoid RNA-seq (e.g. later ISS missions), but it works
for any comparison of two conditions in human midbrain or cortical organoids.

## 1. Install
```bash
pip install git+https://github.com/elenipatsiala/mug-neurotune     # Python ≥ 3.10; numpy, scipy, pandas, matplotlib
mugnt --help
```

## 2. What you need
**A raw count table**: genes in rows, one column per organoid (or per pseudobulked sample).
* Row IDs can be HGNC symbols, Ensembl gene IDs or NCBI GeneIDs. For Ensembl or NCBI IDs, add
  `--annotation FILE` with an ID → `Symbol` table (GEO's `Human.GRCh38.p13.annot.tsv.gz`, or a BioMart export).
* featureCounts output works as is. Use **raw counts**, not TPM or normalised values.
* Single-cell data? Make one profile per organoid first:
  `mugnt pseudobulk 10x_folder/ --cells cells.csv --out counts.csv` (`cells.csv`: `barcode,sample`).

**A metadata CSV** with one row per sample:
```
sample,condition,line,microglia
GSM0001,ground,CTRL_A,noMG
GSM0002,LEO,CTRL_A,noMG
...
```
* `condition`: `ground` / `LEO` (the reference and test labels).
* `line`: donor or cell line. Comparisons are made **within** each line × microglia group, so a donor
  difference can never look like a flight effect. With ≥ 2 lines, a change must point the same way in every line.
* `microglia`: any grouping label; use one constant value (e.g. `none`) if not applicable.
* Only samples listed in the metadata are used. Others in the count table are ignored and reported.

## 3. Run
```bash
# condition comparison + verdict (permutation tests, bootstrap CIs, FDR, hypoxia/glycolysis-adjusted test)
mugnt compare counts.tsv.gz --annotation annot.tsv.gz --meta meta.csv --out results/mission2

# per-sample scores only (no metadata needed)
mugnt score counts.tsv.gz --annotation annot.tsv.gz --out results/scores

# how many donors / organoids would the next mission need? (uses the variance in your data)
mugnt power counts.tsv.gz --annotation annot.tsv.gz --meta meta.csv --program age_clock --out results/power
```

## 4. Read the output (`results/mission2/`)
| File | What it tells you |
|---|---|
| `summary.json` | **Verdict** and reason, adjusted tests, robustness checks, gene-set version and fingerprint |
| `leo_vs_ground.csv` | Per program: Hedges' g, raw difference with bootstrap 95% CI, permutation p, BH q, `consistent_lines` |
| `scores_rank.csv` | Per-sample scores (maturity panel, identity, stress, hypoxia, p53, glia, microglia QC) |
| `maturation_fidelity_map.png` | Each organoid on maturity × dopaminergic fidelity, coloured by hypoxia |
| `leo_vs_ground_zscore_sensitivity.csv` | The same comparison with a different scoring method (signs should agree) |
| `maturation_leave_one_gene_out.csv` | Whether any single marker gene drives the synaptic-marker result |

**Verdicts:**
* **promising maturation**: ≥ 2 of the 3 validated maturity read-outs up, none down, identity kept, no stress program up
* **mixed state**: maturation up together with glycolysis / ER stress / hypoxia / p53 / reactive glia
* **dysregulated state**: maturation up but dopaminergic identity or fidelity down
* **inconsistent maturity read-outs**: the validated read-outs point in opposite directions
* **insufficient evidence**: no robust maturation change

A program counts as changed only if |g| ≥ 0.5, the CI excludes 0, q < 0.05, and (with ≥ 2 lines) the
direction agrees in every line. On 800 simulated missions of mission-1 size, these rules gave the right call
92.5–100% of the time.

## 5. The maturity panel (what "maturation" means here)
| Read-out | What it is | Tracks culture age in 3 labs (ρ) |
|---|---|---|
| `age_clock` | 100 genes that change with age in fetal midbrain neurons in two independent fetal datasets | 0.74 / 0.99 / 0.95 |
| `fetal_corr_age` | whole-transcriptome similarity to fetal ventral midbrain at 6–11.5 PCW (estimated PCW) | 0.93 / 0.84 / 0.84 |
| `progenitor_loss` | decline of SOX2 + NES | 0.94 / 0.98 / 0.95 |

Do **not** read maturity from NeuN/RBFOX3, SYP, TH or MAP2 alone. In organoids they track how many
neurons there are, not how old they are (`docs/BENCHMARK.md`).

## 6. Python API
```python
from mug_neurotune import io, scoring, compare, maturity
counts = io.load_counts("counts.tsv.gz", annotation="annot.tsv.gz")
meta = io.load_metadata("meta.csv")
logcpm = io.log_cpm(counts)[meta.index]       # keep the samples listed in the metadata
scores = scoring.score_all(logcpm)            # programs + maturity panel
res = compare.compare_programs(scores, meta)  # stratified effects, CIs, q-values
print(compare.classify(res))
```

## 7. Limits to keep in mind
* Transcription only: no protein, electrophysiology or dopamine release.
* The maturity panel is validated on **midbrain** organoid time courses (day 15–120). For cortical organoids,
  treat it as provisional.
* `fetal_corr_age` is calibrated only within 6–11.5 PCW. Older organoids are extrapolated.
* With fewer than ~4 donors, results describe *these* lines; see `mugnt power` before generalising.

Questions or problems: open an issue at https://github.com/elenipatsiala/mug-neurotune/issues
