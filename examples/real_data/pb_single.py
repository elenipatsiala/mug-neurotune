"""Per-sample 10x matrix (genes x cells, no barcodes needed) -> whole + neuron-only pseudobulk."""
import gzip, sys, os, glob
import numpy as np, pandas as pd
src, feat, out = sys.argv[1:4]
f = pd.read_csv(feat, sep="\t", header=None)
sym = f[1].astype(str).values
NEU = ["SNAP25", "STMN2", "SYT1", "GAP43", "ELAVL4"]; GLIA = ["GFAP", "SLC1A3", "AQP4", "S100B"]
nm, gm = np.isin(sym, NEU), np.isin(sym, GLIA)
whole, neur, comp = {}, {}, {}
for path in sorted(glob.glob(os.path.join(src, "GSM*_matrix.mtx.gz"))):
    name = os.path.basename(path).split("_")[0]
    with gzip.open(path, "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        ng, nc, _ = map(int, line.split())
        d = pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.float32})
    g, c, v = d.g.values - 1, d.c.values - 1, d.v.values
    if ng != len(sym) and nc == len(sym):   # cells x genes orientation
        g, c = c, g; ng, nc = nc, ng
    lib = np.bincount(c, weights=v, minlength=nc)
    neu = np.bincount(c[nm[g]], weights=v[nm[g]], minlength=nc) / np.maximum(lib, 1) * 1e4
    gl = np.bincount(c[gm[g]], weights=v[gm[g]], minlength=nc) / np.maximum(lib, 1) * 1e4
    keep = lib >= 500
    is_neu = (neu >= 10) & (gl < 2) & keep
    whole[name] = np.bincount(g[keep[c]], weights=v[keep[c]], minlength=ng)
    neur[name] = np.bincount(g[is_neu[c]], weights=v[is_neu[c]], minlength=ng)
    comp[name] = {"cells": int(keep.sum()), "neuron_frac": float(is_neu.sum() / max(keep.sum(), 1)),
                  "glia_frac": float(((gl >= 2) & (neu < 10) & keep).sum() / max(keep.sum(), 1))}
    print(name, comp[name], flush=True)
for nmx, dct in (("held_pseudobulk", whole), ("held_neuron_pseudobulk", neur)):
    df = pd.DataFrame(dct, index=sym).groupby(level=0).sum(); df.index.name = "gene"
    df.to_csv(os.path.join(out, nmx + ".csv.gz"))
pd.DataFrame(comp).T.rename_axis("sample").to_csv(os.path.join(out, "held_composition_noday.csv"))
