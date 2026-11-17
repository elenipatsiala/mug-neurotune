# Data guide: exactly which files to download

The analysis workspace that built this repo had no internet access, so **no real data
has been analysed yet**. Every result in `results/` so far is from synthetic data
(clearly labelled) or from the physics model. Download these four things, put them in
`data/`, and run the commands at the bottom.

## 1. ISS dopaminergic organoids (the target dataset)
* **What:** Marotta et al. 2024, *Stem Cells Transl Med* 13:1186. One bulk RNA-seq library
  per organoid; 20 LEO + 22 ground organoids across cortical and dopaminergic types;
  DA lines HDF410iPS504 (sporadic PD, male) and UEC741iPS517 (control, male);
  ± isogenic microglia.
* **Where:** GEO **GSE259421** → https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259421
* **Easiest download:** on the GEO page, click **"Download RNA-Seq counts"** (NCBI-generated)
  and download:
  * `GSE259421_raw_counts_GRCh38.p13_NCBI.tsv.gz` → `data/iss_counts.tsv.gz`
  * `Human.GRCh38.p13.annot.tsv.gz` → `data/human_annot.tsv.gz` (maps GeneID → Symbol)
  * If the button is missing, use the authors' supplementary count file instead and
    tell me its column layout.
* **Metadata:** open the "Samples (GSM…)" list and build `data/iss_meta.csv` with columns
  `sample,condition,line,microglia,organoid_type` using values `ground/LEO`,
  `PD/CTRL`, `MG/noMG`, `DA/cortical`. Keep only `organoid_type == DA` for the main run.
  Also check NASA OSDR for a matching OSD entry (the paper's dataset DOI 10.26030/c4yh-8e94).

## 2. Known-age midbrain organoids (Step A: validation)
* **What:** Fiorenzano et al. 2021 *Nat Commun* 12:7302 (PMID 34911939). 10x scRNA-seq of
  ventral-midbrain organoids at **days 15, 30, 60, 90, 120**, several organoids per day
  (n = 2, 4–5, 5, 2, 6).
* **Where:** the paper's "Data availability" section gives the GEO accession. Download the
  per-sample 10x folders (matrix.mtx, features.tsv, barcodes.tsv) into
  `data/fiorenzano/<sample>/` and write `data/fiorenzano_age.csv` (`sample,day`).
* `mug_neurotune.io.pseudobulk_10x()` sums cells per sample so they compare with bulk.

## 3. Human fetal ventral midbrain (Step B: anchoring)
* **Option 1:** Birtele et al. 2022 *Development* 149:dev200504 (PMID 36305490). Fetal VM at
  6, 8, 11 post-conception weeks (PCW).
* **Option 2 (wider range):** "An in vivo and in vitro spatiotemporal profile of human midbrain
  development", *Nat Commun* 2025 (doi 10.1038/s41467-025-67779-1). First and second trimester
  fetal midbrain. It also gives our **cross-check**: day-40 organoids ≈ 11 PCW,
  day-70/120 organoids ≈ 15 PCW.
* Pseudobulk per donor/age → `data/fetal_counts.csv` + `data/fetal_age.csv` (`sample,pcw`).

## 4. Commands
```bash
pip install -e .
# A: validation on known-age organoids
mugnt validate data/fiorenzano_pseudobulk.csv --age data/fiorenzano_age.csv --out results/A
# B: fetal anchoring (ISS samples and Fiorenzano samples both)
mugnt anchor data/iss_counts.tsv.gz --annotation data/human_annot.tsv.gz \
     --ref data/fetal_counts.csv --ref-age data/fetal_age.csv --out results/B
# C + D: LEO vs ground, hypoxia-adjusted test, classification
mugnt compare data/iss_counts.tsv.gz --annotation data/human_annot.tsv.gz \
     --meta data/iss_meta_DA.csv --out results/CD
# E: power for the next mission
mugnt power data/iss_counts.tsv.gz --annotation data/human_annot.tsv.gz \
     --meta data/iss_meta_DA.csv --out results/E
```
