import gzip, sys, re
import numpy as np, pandas as pd
f = sys.argv[1]; out = sys.argv[2]
with gzip.open(f, "rt") as fh:
    lines = fh.read().split("\n")
cid = lines[1].split("\t")[2:]; ctype = lines[2].split("\t")[2:]; tp = lines[3].split("\t")[2:]
n = len(cid)
genes, rows = [], []
for l in lines[5:]:
    if not l.strip():
        continue
    p = l.split("\t")
    genes.append(re.sub(r"_(p\d+|loc\d+)$", "", p[0]))
    rows.append(np.array(p[2:2 + n], dtype=np.float32))
X = pd.DataFrame(np.vstack(rows), index=genes, columns=cid).groupby(level=0).sum()
meta = pd.DataFrame({"cell_type": ctype, "week": [float(t.replace("week_", "")) for t in tp],
                     "chip": [c.rsplit("_", 1)[0] for c in cid]}, index=cid)
NEURON = re.compile(r"^h(DA\d|Nb|OMTN|RN|Sert|Gaba|NbGaba)")
meta["is_neuron"] = meta.cell_type.str.match(NEURON)
print(meta.groupby(["week", "chip"]).agg(cells=("cell_type", "size"), neurons=("is_neuron", "sum")).to_string())
print(meta.cell_type.value_counts().to_string())
# pseudobulk per chip (whole tissue) and per chip (neurons only)
pb = X.T.groupby(meta.chip).sum().T
pbn = X.loc[:, meta.is_neuron].T.groupby(meta.chip[meta.is_neuron]).sum().T
ages = meta.groupby("chip").week.first()
pb.index.name = pbn.index.name = "gene"
pb.to_csv(out + "/lamanno_tissue_pseudobulk.csv.gz"); pbn.to_csv(out + "/lamanno_neuron_pseudobulk.csv.gz")
ages.rename("pcw").to_csv(out + "/lamanno_age.csv")
meta.to_csv(out + "/lamanno_cells_meta.csv.gz")
X.to_pickle(out + "/lamanno_cells.pkl")
