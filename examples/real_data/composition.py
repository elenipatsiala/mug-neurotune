"""Per organoid: neuron / glia / progenitor fractions, and pseudobulk of NEURONS ONLY."""
import glob, os, re, sys
import numpy as np, pandas as pd
src, out = sys.argv[1], sys.argv[2]
NEU = ["SNAP25", "STMN2", "SYT1", "GAP43", "ELAVL4"]
GLIA = ["GFAP", "SLC1A3", "AQP4", "S100B"]
CYC = ["MKI67", "TOP2A"]
comp, neu_pb = {}, {}
for f in sorted(glob.glob(os.path.join(src, "GSM*.csv.gz"))):
    name = os.path.basename(f).split(".")[0]
    try:
        df = pd.read_csv(f, index_col=0)
    except Exception:
        continue
    lib = df.sum(0)
    def s(genes):
        g = [x for x in genes if x in df.index]
        return (df.loc[g].sum(0) / lib * 1e4)
    neu, gl, cy = s(NEU), s(GLIA), s(CYC)
    is_neu = (neu >= 10) & (gl < 2)
    is_gl = (gl >= 2) & (neu < 10)
    comp[name] = {"day": int(re.search(r"day(\d+)", name).group(1)), "cells": df.shape[1],
                  "neuron_frac": is_neu.mean(), "glia_frac": is_gl.mean(), "cycling_frac": (cy >= 1).mean()}
    neu_pb[name] = df.loc[:, is_neu.values].sum(1)
    print(name, {k: round(v, 3) if isinstance(v, float) else v for k, v in comp[name].items()}, flush=True)
pd.DataFrame(comp).T.rename_axis("sample").to_csv(os.path.join(out, "fiorenzano_composition.csv"))
pb = pd.DataFrame(neu_pb); pb.index.name = "gene"
pb.to_csv(os.path.join(out, "fiorenzano_neuron_pseudobulk.csv"))
