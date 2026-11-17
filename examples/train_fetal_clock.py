"""Learn a neuronal age clock from fetal ventral-midbrain NEURONS (GSE192405 tissue samples).

Training data only: 6, 8, 10.5, 11.5 PCW tissue. Organoid data are NOT touched here.
Per gene: Spearman rho between expression (log CP10k) and PCW across neurons, with
neurons subsampled to the same number per sample so no sample dominates.
Clock = top-N age-UP genes minus top-N age-DOWN genes (background-referenced rank scores).
Genes excluded: ribosomal, mitochondrial, haemoglobin, sex-chromosome, and all
v0.1 program genes (keeps the clock independent of the stress/hypoxia programs).
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/home/claude/mug-neurotune")
from mug_neurotune import genesets

SRC = "/home/claude/data_fetal"
TISSUE = {"GSM5746451_MP06-hVM-6wks": 6.0, "GSM5746448_MP03-hVM-8wks": 8.0,
          "GSM5746447_MP02-hVM-10-5wks": 10.5, "GSM5746446_MP01-hVM-11-5wks": 11.5}
NEU, GLIA = ["SNAP25", "STMN2", "SYT1", "GAP43", "ELAVL4"], ["GFAP", "SLC1A3", "AQP4", "S100B"]
N_PER = 400
N_GENES = 50
rng = np.random.default_rng(0)
cells, ages, comp = [], [], {}
def two_pass(path, n_pick):
    """Pass 1: library size + marker sums per cell. Pass 2: load only picked cells."""
    lib = None; neu = None; gl = None
    for ch in pd.read_csv(path, index_col=0, chunksize=1500, dtype={0: str}):
        ch = ch.astype(np.float32)
        s = ch.sum(axis=0); lib = s if lib is None else lib + s
        n = ch.loc[ch.index.isin(NEU)].sum(axis=0); neu = n if neu is None else neu + n
        g = ch.loc[ch.index.isin(GLIA)].sum(axis=0); gl = g if gl is None else gl + g
    neu, gl = neu / lib * 1e4, gl / lib * 1e4
    is_neu = (neu >= 10) & (gl < 2)
    idx = np.where(is_neu.values)[0]
    pick = lib.index[rng.choice(idx, min(n_pick, len(idx)), replace=False)] if len(idx) else []
    parts = [ch[list(pick)].astype(np.float32) for ch in pd.read_csv(path, index_col=0, chunksize=1500)]
    sub = pd.concat(parts)
    sub = np.log1p(sub.div(lib[list(pick)], axis=1) * 1e4)
    return sub, int(len(lib)), int(is_neu.sum())


for name, pcw in TISSUE.items():
    sub, ncell, nneu = two_pass(os.path.join(SRC, name + ".csv.gz"), N_PER)
    comp[name] = {"pcw": pcw, "cells": ncell, "neurons": nneu}
    cells.append(sub)
    pick = range(sub.shape[1])
    ages += [pcw] * len(pick)
    print(name, comp[name], flush=True)
X = pd.concat(cells, axis=1).fillna(0)
ages = np.array(ages)
program = {g for v in genesets.gene_lists().values() for g in v}
bad = X.index.str.match(r"^(RP[SL]\d|MT-|MTRNR2L|HB[AB]|XIST|RPS4Y|DDX3Y|UTY|KDM5D|EIF1AY|MALAT1|NEAT1|HSP|DNAJ)")
# Dissociation-induced immediate-early/heat-shock genes (van den Brink et al. 2017 Nat Methods core set)
DISSOC = {"FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "IER2", "IER3", "DUSP1", "ZFP36",
          "NR4A1", "SOCS3", "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSPB1", "DNAJA1", "DNAJB1",
          "UBC", "BTG2", "KLF6", "CEBPD", "PPP1R15A"}
# Radial-glia / astroglial genes (keep the clock neuronal even if a few glia pass the neuron gate)
GLIAL = {"FABP7", "SLC1A3", "NFIA", "NFIB", "VIM", "C1orf61", "PTPRZ1", "HES1", "HES5", "SOX2", "SOX9",
         "TTYH1", "GFAP", "AQP4", "S100B", "PAX6", "NES", "MKI67", "TOP2A", "CLU", "SPARC", "SPARCL1"}
bad = bad | X.index.isin(DISSOC | GLIAL)
expressed = (X > 0).mean(axis=1) >= 0.10
X = X[expressed & ~bad & ~X.index.isin(program)]
rho = pd.Series(stats.spearmanr(X.T.values, ages[:, None]).statistic[-1, :-1], index=X.index)
rho = rho.dropna().sort_values()
up = list(rho.index[-N_GENES:][::-1])
down = list(rho.index[:N_GENES])
clock = {"name": "neuronal_age_clock", "version": "0.2.0",
         "trained_on": "GSE192405 fetal VM tissue neurons, 6/8/10.5/11.5 PCW, 400 neurons each",
         "up": up, "down": down,
         "rho_up_range": [float(rho[up].min()), float(rho[up].max())],
         "rho_down_range": [float(rho[down].min()), float(rho[down].max())],
         "n_genes_tested": int(len(rho)), "composition": comp}
json.dump(clock, open("/home/claude/work/neuronal_age_clock.json", "w"), indent=1)
print("UP:", up[:25]); print("DOWN:", down[:25])
print("rho up", clock["rho_up_range"], "down", clock["rho_down_range"])
