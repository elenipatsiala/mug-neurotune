"""Benchmark figure: which maturity read-outs track organoid age across three labs?"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
B = "results/real/benchmark"
sets = [("lund_whole", "Lund\n(GSE168323)\nday 15–120"), ("brainstem_whole", "BrainSTEM\n(GSE281535)\nday 20–60"),
        ("heldout_whole", "Held-out\n(GSE277032)\nday 40–120")]
tabs = {k: pd.read_csv(f"{B}/benchmark_{k}.csv", index_col=0) for k, _ in sets}
# progenitor-loss panel read-out
from mug_neurotune import io, maturity
from scipy import stats
paths = {"lund_whole": ("results/real/derived/fiorenzano_pseudobulk.csv.gz", "results/real/derived/fiorenzano_age.csv", "results/real/derived/fiorenzano_composition.csv"),
         "brainstem_whole": ("/home/claude/work/brainstem_pseudobulk.csv.gz", "/home/claude/work/brainstem_age.csv", "/home/claude/work/brainstem_composition.csv"),
         "heldout_whole": ("/home/claude/work/held_pseudobulk.csv.gz", "/home/claude/work/held_age.csv", "/home/claude/work/held_composition.csv")}
for k, (p, a, c) in paths.items():
    lc = io.log_cpm(pd.read_csv(p, index_col=0)); age = pd.read_csv(a, index_col=0).iloc[:, 0]
    comp = pd.read_csv(c, index_col=0)
    s = maturity.progenitor_loss(lc)
    tabs[k].loc["progenitor loss (SOX2+NES)"] = [stats.spearmanr(s, age.loc[s.index])[0], np.nan,
                                                 stats.spearmanr(s, comp.loc[s.index, "neuron_frac"])[0]]
order = ["age clock v0.3 (this tool)", "fetal-correlation age", "progenitor loss (SOX2+NES)", "marker SOX2 (−)",
         "marker NES (−)", "marker MKI67 (−)", "age clock v0.2", "v0.1 synaptic set", "marker RBFOX3", "marker SYP",
         "marker TH", "marker MAP2"]
label = {"age clock v0.3 (this tool)": "Age clock v0.3 ★", "fetal-correlation age": "Fetal-correlation age ★",
         "progenitor loss (SOX2+NES)": "Progenitor loss SOX2+NES ★", "marker SOX2 (−)": "SOX2 ↓",
         "marker NES (−)": "NES ↓", "marker MKI67 (−)": "MKI67 ↓", "age clock v0.2": "Age clock v0.2 (retired)",
         "v0.1 synaptic set": "Synaptic set v0.1 (retired)", "marker RBFOX3": "NeuN / RBFOX3 ↑",
         "marker SYP": "SYP ↑", "marker TH": "TH ↑", "marker MAP2": "MAP2 ↑"}
R = pd.DataFrame({k: tabs[k].reindex(order)["rho_day"] for k, _ in sets})
N = pd.DataFrame({k: tabs[k].reindex(order)["rho_neuron_frac"] for k, _ in sets})
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
for ax, M, title, cmap, vmin in ((axes[0], R, "A  Tracks culture age? (Spearman ρ with day)", "RdBu", -1),
                                (axes[1], N, "B  Fooled by neuron content? (ρ with neuron fraction)", "PuOr_r", -1)):
    im = ax.imshow(M.values, cmap=cmap, vmin=vmin, vmax=1, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in sets], fontsize=8.5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([label[o] for o in order] if ax is axes[0] else [])
    for i in range(M.shape[0]):
        for j in range(3):
            v = M.values[i, j]
            if np.isfinite(v):
                ok = (ax is axes[0] and v >= 0.6)
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if abs(v) > 0.6 else "#1c2230", fontweight="bold" if ok else "normal")
    ax.axhline(2.5, color="#1c2230", lw=1.2)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8)
passes = (R >= 0.6).all(axis=1)
for i, o in enumerate(order):
    axes[0].text(3.55, i, "✓ all 3" if passes[o] else "✗", va="center", fontsize=9,
                 color="#15803d" if passes[o] else "#b91c1c", clip_on=False)
fig.text(0.01, -0.02, "★ = validated maturity panel used by µG-NeuroTune v0.3 (a change is called only when its read-outs agree). "
         "Bold = passes ρ ≥ 0.6. Whole-organoid pseudobulk. Held-out = tested once, after freezing the clock (docs/FREEZE_v0.3.md).",
         fontsize=8, color="#64748b")
fig.tight_layout()
fig.savefig("results/real/fig_benchmark.png", dpi=200, bbox_inches="tight")
R.to_csv("results/real/benchmark/summary_rho_day.csv"); N.to_csv("results/real/benchmark/summary_rho_neuronfrac.csv")
print(R.round(2).assign(all3=passes).to_string())
