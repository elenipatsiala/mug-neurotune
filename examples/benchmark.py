"""Benchmark: how well do common maturity read-outs track organoid age?

Methods (what labs use today vs this tool):
  * single markers (higher = mature): RBFOX3 (NeuN), MAP2, SYP, TH
  * single markers (lower = mature, sign flipped): MKI67, SOX2, NES
  * v0.1 synaptic marker set (rank score)
  * fetal-correlation age: whole-transcriptome Spearman correlation to fetal tissue
    pseudobulk, softmax-weighted age (the reference-mapping idea behind VoxHunt-style tools)
  * neuronal age clock (this tool)
Metrics on each known-age organoid dataset:
  rho_day           Spearman with culture day (higher = better)
  rho_day|neurons   partial Spearman with day, controlling for neuron fraction (robust to composition)
  rho_neuron_frac   Spearman with neuron fraction (closer to 0 = less fooled by composition)
Usage: python examples/benchmark.py DATASET_NAME PSEUDOBULK.csv AGES.csv COMPOSITION.csv OUTDIR
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mug_neurotune import io, scoring, clock, anchor

name, pb_path, age_path, comp_path, out = sys.argv[1:6]
FETAL = os.environ.get("FETAL_PB", "results/real/derived/fetal_tissue_pseudobulk.csv.gz")
FETAL_AGE = os.environ.get("FETAL_AGE", "results/real/derived/fetal_age.csv")
os.makedirs(out, exist_ok=True)
pb = pd.read_csv(pb_path, index_col=0)
age = pd.read_csv(age_path, index_col=0).iloc[:, 0].loc[pb.columns]
comp = pd.read_csv(comp_path, index_col=0).loc[pb.columns]
lc = io.log_cpm(pb)
scores = {}
for g, sign in [("RBFOX3", 1), ("MAP2", 1), ("SYP", 1), ("TH", 1), ("MKI67", -1), ("SOX2", -1), ("NES", -1)]:
    if g in lc.index:
        scores[f"marker {g}" + (" (−)" if sign < 0 else "")] = sign * lc.loc[g]
scores["v0.1 synaptic set"] = scoring.score_all(lc)["maturation"]
fet = io.log_cpm(pd.read_csv(FETAL, index_col=0))
fage = pd.read_csv(FETAL_AGE, index_col=0).iloc[:, 0]
shared = lc.index.intersection(fet.index)
expressed = shared[(lc.loc[shared] > 1).mean(1) > 0.5]
genes = anchor.age_informative_genes(fet.loc[expressed], fage, n_genes=300)
scores["fetal-correlation age"] = anchor.estimate_age(lc, fet, fage, genes=genes)["est_pcw"]
scores["age clock v0.2"] = clock.clock_score(lc)["age_clock"]
import json as _j
_orig = clock.load_clock
clock.load_clock = lambda: _j.load(open(os.path.join(os.path.dirname(clock.__file__), "data", "neuronal_age_clock_v03.json")))
scores["age clock v0.3 (this tool)"] = clock.clock_score(lc)["age_clock"]
clock.load_clock = _orig

def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    res = lambda a, b: a - np.polyval(np.polyfit(b, a, 1), b)
    return stats.pearsonr(res(rx, rz), res(ry, rz))[0]

rows = []
for m, s in scores.items():
    s = s.loc[pb.columns].astype(float)
    rows.append({"method": m, "rho_day": stats.spearmanr(s, age)[0],
                 "rho_day|neurons": partial_spearman(s.values, age.values, comp["neuron_frac"].values),
                 "rho_neuron_frac": stats.spearmanr(s, comp["neuron_frac"])[0]})
res = pd.DataFrame(rows).set_index("method")
res.to_csv(os.path.join(out, f"benchmark_{name}.csv"))
pd.DataFrame(scores).loc[pb.columns].assign(day=age).to_csv(os.path.join(out, f"benchmark_scores_{name}.csv"))
print(f"== {name}: n={len(pb.columns)} organoids, days {sorted(age.unique())}")
print(res.round(2).sort_values("rho_day", ascending=False).to_string())
