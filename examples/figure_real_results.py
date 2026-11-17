"""Main results figure (real data). Run after the real-data analyses in results/real/."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREY, ORANGE, LIGHT = "#1d4ed8", "#64748b", "#c2410c", "#cbd5e1"
fig = plt.figure(figsize=(15, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 1.0], wspace=0.38)

# A: validation
ax = fig.add_subplot(gs[0])
ck = pd.read_csv("/home/claude/work/clock_fiorenzano_pseudobulk.csv", index_col=0)
old = pd.read_csv("results/real/A/validation.csv", index_col=0)
from mug_neurotune import io, scoring
pb = pd.read_csv("/home/claude/work/fiorenzano_pseudobulk.csv", index_col=0)
mat = scoring.score_all(io.log_cpm(pb))["maturation"]
z = lambda s: (s - s.mean()) / s.std()
jit = np.random.default_rng(0).uniform(-4, 4, len(ck))
ax.scatter(ck.day + jit, z(ck.age_clock), color=BLUE, s=30, label="age clock (v0.2)", zorder=3)
ax.scatter(ck.day + jit, z(mat.loc[ck.index]), color=LIGHT, edgecolor=GREY, s=30, label="marker set (v0.1)")
for s_, col in ((z(ck.age_clock), BLUE), (z(mat.loc[ck.index]), GREY)):
    m = s_.groupby(ck.day).mean(); ax.plot(m.index, m.values, color=col, lw=1.5)
r_new = stats.spearmanr(ck.day, ck.age_clock)[0]; r_old = stats.spearmanr(ck.day, mat.loc[ck.index])[0]
ax.set_xticks([15, 30, 60, 120]); ax.set_xlabel("days in culture (independent organoids, GSE168323)")
ax.set_ylabel("score (z)")
ax.set_title(f"A  Validation: clock ρ = {r_new:.2f} ✓ · markers ρ = {r_old:.2f} ✗", loc="left", fontweight="bold", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="upper left")

# B: ISS effects
ax = fig.add_subplot(gs[1])
r = pd.read_csv("results/real/CD/DA_with_clock.csv", index_col=0)
rc = pd.read_csv("results/real/CD/cortical_with_clock.csv", index_col=0)
progs = ["age_clock", "maturation", "da_identity", "progenitor_cycle", "stress_glycolysis", "hypoxia", "p53_dna_damage"]
names = ["Age clock ✓", "Synaptic markers", "DA identity", "Progenitor/cycle", "Glycolysis", "Hypoxia", "p53/DNA damage"]
y = np.arange(len(progs))[::-1]
for k, (p, yy) in enumerate(zip(progs, y)):
    g = r.loc[p, "hedges_g"]; sig = r.loc[p, "q_bh"] < 0.05 and r.loc[p, "consistent_lines"]
    ax.barh(yy + 0.18, g, height=0.34, color=(ORANGE if sig else LIGHT), edgecolor=GREY, lw=0.5)
    ax.barh(yy - 0.18, rc.loc[p, "hedges_g"], height=0.34, color="none", edgecolor=GREY, lw=1, hatch="///")
    ax.text(max(g, 0) + 0.15 if g >= 0 else g - 0.15, yy + 0.18, f"q={r.loc[p,'q_bh']:.3f}", va="center",
            ha="left" if g >= 0 else "right", fontsize=7.5, color="#1c2230" if sig else "#94a3b8")
ax.axvline(0, color="k", lw=0.7); ax.set_yticks(y); ax.set_yticklabels(names)
ax.set_xlabel("ISS − ground effect (Hedges' g, within line × microglia)")
ax.set_xlim(-4.8, 3.0)
ax.set_title("B  ISS organoids: metabolic shift, no faster aging", loc="left", fontweight="bold", fontsize=10)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=ORANGE, label="dopaminergic, significant"), Patch(color=LIGHT, label="dopaminergic, n.s."),
                   Patch(facecolor="none", edgecolor=GREY, hatch="///", label="cortical (replication)")],
          frameon=True, framealpha=0.95, edgecolor="none", fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=3)

# C: power
ax = fig.add_subplot(gs[2])
G = pd.read_csv("results/real/E_power_clock.csv")
P = G[(G.effect_days == 14) & (G.dxc_ratio == 0.5)].pivot(index="n_donors", columns="organoids_per_donor_per_arm", values="power")
im = ax.imshow(P.values, cmap="Blues", vmin=0, vmax=1, aspect="auto", origin="lower")
ax.set_xticks(range(P.shape[1])); ax.set_xticklabels(P.columns); ax.set_yticks(range(P.shape[0])); ax.set_yticklabels(P.index)
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        v = P.values[i, j]; ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8, color="white" if v > 0.6 else "#1c2230",
                                    fontweight="bold" if v >= 0.8 else "normal")
ax.add_patch(plt.Rectangle((-0.5, -0.5), P.shape[1], 1, fill=False, edgecolor=ORANGE, lw=2))
ax.text(P.shape[1] - 0.45, 0, " CRS-19", color=ORANGE, va="center", fontsize=8, clip_on=False)
ax.set_xlabel("organoids per donor per arm"); ax.set_ylabel("donors (cell lines)")
ax.set_title("C  Power to detect 2 weeks faster maturation", loc="left", fontweight="bold", fontsize=10)
fig.text(0.01, -0.12, "Data: GSE259421 (ISS, 18 DA + 18 cortical organoids), GSE168323 (organoid time-course, 16 organoids), GSE192405 (fetal ventral midbrain, 6–11.5 PCW). "
         "Clock trained on fetal neurons only, tested blind on organoids. q = BH-adjusted within-stratum permutation test.", fontsize=7.5, color=GREY)
fig.savefig("results/real/fig_main_results.png", dpi=200, bbox_inches="tight")
print("ok")
