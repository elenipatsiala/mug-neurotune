import os, sys, numpy as np, pandas as pd
SRC = "/home/claude/data_fetal"
T = {"GSM5746451_MP06-hVM-6wks": 6.0, "GSM5746448_MP03-hVM-8wks": 8.0,
     "GSM5746447_MP02-hVM-10-5wks": 10.5, "GSM5746446_MP01-hVM-11-5wks": 11.5}
cols = {}
for n in T:
    tot = None
    for ch in pd.read_csv(os.path.join(SRC, n + ".csv.gz"), index_col=0, chunksize=1500):
        s = ch.sum(axis=1); tot = s if tot is None else pd.concat([tot, s])
    cols[n] = tot
pb = pd.DataFrame(cols).fillna(0); pb.index.name = "gene"
pb.to_csv("/home/claude/work/fetal_tissue_pseudobulk.csv")
pd.Series(T, name="pcw").rename_axis("sample").to_csv("/home/claude/work/fetal_age.csv")
print(pb.shape)
