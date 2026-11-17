import glob, os, re, sys
import pandas as pd
src = sys.argv[1]; out = sys.argv[2]
cols, ages, ncells = {}, {}, {}
for f in sorted(glob.glob(os.path.join(src, "GSM*.csv.gz"))):
    name = os.path.basename(f).split(".")[0]
    try:
        tot = None; n = None
        for chunk in pd.read_csv(f, index_col=0, chunksize=2000):
            s = chunk.sum(axis=1)
            tot = s if tot is None else pd.concat([tot, s])
            n = chunk.shape[1]
    except Exception as e:
        print("skip", name, e); continue
    cols[name] = tot; ncells[name] = n
    ages[name] = int(re.search(r"day(\d+)", name).group(1))
    print(name, ages[name], n, "cells", int(tot.sum()), "UMIs", flush=True)
pb = pd.DataFrame(cols).fillna(0)
pb.index.name = "gene"
pb.to_csv(os.path.join(out, "fiorenzano_pseudobulk.csv"))
pd.DataFrame({"day": ages, "n_cells": ncells}).rename_axis("sample").to_csv(os.path.join(out, "fiorenzano_age.csv"))
