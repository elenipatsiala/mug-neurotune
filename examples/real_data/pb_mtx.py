"""Stream a 10x-style matrix.mtx.gz (genes x cells) into per-sample pseudobulk,
plus neuron-only pseudobulk and composition. Two passes, constant memory."""
import gzip, sys, os
import numpy as np, pandas as pd
mtx, genes_f, bc_f, out = sys.argv[1:5]
genes = pd.read_csv(genes_f, sep="\t", header=None)
sym = genes[1 if genes.shape[1] > 1 else 0].astype(str).values
bcs = pd.read_csv(bc_f, sep="\t", header=None)[0].astype(str)
day = bcs.str.extract(r"^D(\d+)")[0].astype(int).values
samp = (bcs.str.split("__").str[0].str.split("_").str[0] + "_" + bcs.str.split("__").str[1]).values
groups = pd.Index(sorted(set(samp)))
gidx = groups.get_indexer(samp)
NEU = ["SNAP25", "STMN2", "SYT1", "GAP43", "ELAVL4"]; GLIA = ["GFAP", "SLC1A3", "AQP4", "S100B"]
neu_rows = set(np.where(np.isin(sym, NEU))[0]); gl_rows = set(np.where(np.isin(sym, GLIA))[0])
neu_mask_g = np.isin(np.arange(len(sym)), list(neu_rows)); gl_mask_g = np.isin(np.arange(len(sym)), list(gl_rows))
n_cells = len(bcs)

def chunks():
    with gzip.open(mtx, "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        for ch in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=5_000_000,
                              dtype={"g": np.int32, "c": np.int32, "v": np.float32}):
            yield ch.g.values - 1, ch.c.values - 1, ch.v.values

lib = np.zeros(n_cells); neu = np.zeros(n_cells); gl = np.zeros(n_cells)
pb = np.zeros((len(sym), len(groups)))
for g, c, v in chunks():
    np.add.at(lib, c, v)
    m = neu_mask_g[g]; np.add.at(neu, c[m], v[m])
    m = gl_mask_g[g]; np.add.at(gl, c[m], v[m])
    np.add.at(pb, (g, gidx[c]), v)
neu_cp, gl_cp = neu / lib * 1e4, gl / lib * 1e4
is_neu = (neu_cp >= 10) & (gl_cp < 2)
pbn = np.zeros((len(sym), len(groups)))
for g, c, v in chunks():
    m = is_neu[c]
    np.add.at(pbn, (g[m], gidx[c[m]]), v[m])
pbdf = pd.DataFrame(pb, index=sym, columns=groups).groupby(level=0).sum(); pbdf.index.name = "gene"
pbndf = pd.DataFrame(pbn, index=sym, columns=groups).groupby(level=0).sum(); pbndf.index.name = "gene"
pbdf.to_csv(os.path.join(out, "brainstem_pseudobulk.csv.gz"))
pbndf.to_csv(os.path.join(out, "brainstem_neuron_pseudobulk.csv.gz"))
comp = pd.DataFrame({"sample": samp, "day": day, "neuron": is_neu, "glia": (gl_cp >= 2) & (neu_cp < 10)})
c2 = comp.groupby("sample").agg(day=("day", "first"), cells=("day", "size"), neuron_frac=("neuron", "mean"), glia_frac=("glia", "mean"))
c2.to_csv(os.path.join(out, "brainstem_composition.csv"))
c2[["day"]].to_csv(os.path.join(out, "brainstem_age.csv"))
print(c2.groupby("day")[["cells", "neuron_frac", "glia_frac"]].mean().round(3))
