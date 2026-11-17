"""Clock v0.3: genes must change with age in the SAME direction in two independent fetal
ventral-midbrain datasets (Birtele 2022 10x, 4 fetuses; La Manno 2016 C1, 35 chips, wk 6-11).
Score per gene = min(|rho_Birtele|, |rho_LaManno|) with matching sign. Organoid data untouched."""
import json, os, sys, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/home/claude/mug-neurotune")
from mug_neurotune import genesets
W = "/home/claude/work"
rng = np.random.default_rng(0)
DISSOC = {"FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","IER2","IER3","DUSP1","ZFP36","NR4A1","SOCS3",
          "HSPA1A","HSPA1B","HSPA8","HSP90AA1","HSPB1","DNAJA1","DNAJB1","UBC","BTG2","KLF6","CEBPD","PPP1R15A"}
GLIAL = {"FABP7","SLC1A3","NFIA","NFIB","VIM","C1orf61","PTPRZ1","HES1","HES5","SOX2","SOX9","TTYH1","GFAP",
         "AQP4","S100B","PAX6","NES","MKI67","TOP2A","CLU","SPARC","SPARCL1"}
PROGRAM = {g for v in genesets.gene_lists().values() for g in v}
BAD = re.compile(r"^(RP[SL]\d|MT-|MTRNR2L|HB[AB]|XIST|RPS4Y|DDX3Y|UTY|KDM5D|EIF1AY|MALAT1|NEAT1|HSP|DNAJ)")

def gene_rho(X, ages, min_frac=0.10):
    keep = (X > 0).mean(axis=1) >= min_frac
    X = X[keep]
    r = stats.spearmanr(X.T.values, ages[:, None]).statistic[-1, :-1]
    return pd.Series(r, index=X.index).dropna()

# --- La Manno neurons (cell level)
X = pd.read_pickle(f"{W}/lamanno_cells.pkl")
meta = pd.read_csv(f"{W}/lamanno_cells_meta.csv.gz", index_col=0)
Xn = X.loc[:, meta.is_neuron.values]
Xn = np.log1p(Xn / Xn.sum(axis=0) * 1e4)
rho_lm = gene_rho(Xn, meta.loc[meta.is_neuron, "week"].values)
print("La Manno neurons:", Xn.shape[1], "genes tested", len(rho_lm), flush=True)

# --- Birtele neurons (cell level, 400/sample, same gate as v0.2)
exec(open(f"{W}/fetal_clock.py").read().split("for name, pcw in TISSUE.items():")[0])
cells, ages = [], []
for name, pcw in TISSUE.items():
    sub, _, _ = two_pass(os.path.join(SRC, name + ".csv.gz"), N_PER)
    cells.append(sub); ages += [pcw] * sub.shape[1]
Xb = pd.concat(cells, axis=1).fillna(0)
rho_b = gene_rho(Xb, np.array(ages))
print("Birtele neurons:", Xb.shape[1], "genes tested", len(rho_b), flush=True)

common = rho_lm.index.intersection(rho_b.index)
common = [g for g in common if not BAD.match(g) and g not in DISSOC | GLIAL | PROGRAM]
df = pd.DataFrame({"rho_birtele": rho_b[common], "rho_lamanno": rho_lm[common]})
same = np.sign(df.rho_birtele) == np.sign(df.rho_lamanno)
df["consensus"] = np.where(same, np.sign(df.rho_birtele) * np.minimum(df.rho_birtele.abs(), df.rho_lamanno.abs()), 0.0)
df.to_csv(f"{W}/clock_v03_gene_table.csv")
up = list(df.consensus.sort_values(ascending=False).index[:50])
down = list(df.consensus.sort_values().index[:50])
print("agreement of signs across datasets (all genes):", round(same.mean(), 3),
      "| corr of rho:", round(stats.spearmanr(df.rho_birtele, df.rho_lamanno)[0], 3))
print("UP:", up[:30]); print("DOWN:", down[:30])
print("consensus |rho| range up", round(df.consensus[up].min(), 3), "down", round(df.consensus[down].max(), 3))
json.dump({"name": "neuronal_age_clock", "version": "0.3.0",
           "trained_on": "consensus of fetal VM neurons: GSE192405 (Birtele 2022; 6/8/10.5/11.5 PCW, 400 neurons each) "
                         "and GSE76381 (La Manno 2016; weeks 6-11, 35 C1 chips, author-annotated neurons)",
           "rule": "same sign in both datasets; score = min(|rho|); top 50 each direction",
           "up": up, "down": down}, open(f"{W}/neuronal_age_clock_v03.json", "w"), indent=1)
