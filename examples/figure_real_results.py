"""Main results figure (v0.3, real data)."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mug_neurotune import io, clock, maturity
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREY, ORANGE, LIGHT, GREEN = "#1d4ed8", "#64748b", "#c2410c", "#cbd5e1", "#15803d"
D = "results/real/derived"
SETS = [("Lund (GSE168323)", f"{D}/fiorenzano_pseudobulk.csv.gz", f"{D}/fiorenzano_age.csv", "#1d4ed8", "o"),
        ("BrainSTEM (GSE281535)", f"{D}/brainstem_pseudobulk.csv.gz", f"{D}/brainstem_age.csv", "#0891b2", "s"),
        ("Held-out (GSE277032)", f"{D}/held_pseudobulk.csv.gz", f"{D}/held_age.csv", "#7c3aed", "D")]
fig = plt.figure(figsize=(15.5, 4.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.35, 0.95], wspace=0.6)

# A: clock v0.3 vs day in three labs
ax = fig.add_subplot(gs[0])
for name, p, a, col, mk in SETS:
    lc = io.log_cpm(pd.read_csv(p, index_col=0)); age = pd.read_csv(a, index_col=0).iloc[:, 0]
    s = clock.clock_score(lc)["age_clock"]; age = age.loc[s.index]
    z = (s - s.mean()) / s.std()
    r = stats.spearmanr(s, age)[0]
    ax.scatter(age, z, s=18, color=col, marker=mk, alpha=0.75, label=f"{name}: ρ = {r:.2f}")
    m = z.groupby(age).mean(); ax.plot(m.index, m.values, color=col, lw=1.4)
ax.set_xscale("log"); ax.set_xticks([15, 20, 30, 40, 60, 120]); ax.set_xticklabels([15, 20, 30, 40, 60, 120])
ax.set_ylim(-1.9, 2.6); ax.set_xlabel("days in culture (log scale)"); ax.set_ylabel("age clock v0.3 (z within dataset)")
ax.set_title("A  Age clock v0.3 works in 3 labs", loc="left", fontweight="bold", fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.0, 1.0))

# B: ISS effects, DA (filled) and cortical (hatched)
ax = fig.add_subplot(gs[1])
r = pd.read_csv("results/real/CD_v03_DA/leo_vs_ground.csv", index_col=0)
rc = pd.read_csv("results/real/CD_v03_cortical/leo_vs_ground.csv", index_col=0)
progs = ["age_clock", "fetal_corr_age", "progenitor_loss", "stress_glycolysis", "hypoxia", "p53_dna_damage"]
names = ["Age clock v0.3 ★", "Fetal-corr. age ★", "Progenitor loss ★", "Glycolysis", "Hypoxia", "p53 / DNA damage"]
y = np.arange(len(progs))[::-1]
for p, yy in zip(progs, y):
    for tab, off, hatch in ((r, 0.18, None), (rc, -0.18, "///")):
        g = tab.loc[p, "hedges_g"]; sig = tab.loc[p, "q_bh"] < 0.05 and tab.loc[p, "consistent_lines"]
        ax.barh(yy + off, g, height=0.34, color=(ORANGE if sig else LIGHT) if hatch is None else "none",
                edgecolor=(ORANGE if sig else GREY), lw=1 if hatch else 0.5, hatch=hatch)
        ax.text(g + (0.12 if g >= 0 else -0.12), yy + off, f"q={tab.loc[p, 'q_bh']:.3f}", va="center",
                ha="left" if g >= 0 else "right", fontsize=7, color="#1c2230" if sig else "#94a3b8")
ax.axvline(0, color="k", lw=0.7); ax.axhline(2.5, color="#94a3b8", lw=0.8, ls=":")
ax.set_yticks(y); ax.set_yticklabels(names); ax.set_xlim(-5, 3.3)
ax.text(3.25, 4.0, "maturity panel:\nread-outs disagree\n→ insufficient\n    evidence", ha="right", va="center", fontsize=7.5, color=GREY)
ax.set_xlabel("ISS − ground (Hedges' g, within line × microglia)")
ax.set_title("B  ISS organoids: no consistent faster maturation;\n    glycolysis ↑, p53 ↓", loc="left", fontweight="bold", fontsize=10)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=ORANGE, label="significant (q < 0.05, both lines)"), Patch(color=LIGHT, label="not significant"),
                   Patch(facecolor="none", edgecolor=GREY, hatch="///", label="cortical organoids")],
          frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=3)

# C: power (2-week effect)
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
ax.text(P.shape[1] - 0.45, 0, " mission 1", color=ORANGE, va="center", fontsize=8, clip_on=False)
ax.set_xlabel("organoids per donor per arm"); ax.set_ylabel("donors (cell lines)")
ax.set_title("C  Power to detect 2 weeks\n    faster maturation", loc="left", fontweight="bold", fontsize=10)
fig.text(0.01, -0.12, "★ validated maturity panel (each read-out tracks culture age in 3 organoid datasets from 3 labs; docs/BENCHMARK.md). "
         "Data: GSE259421 (ISS: 18 DA + 18 cortical organoids), GSE168323, GSE281535, GSE277032 (organoid time courses), "
         "GSE192405 + GSE76381 (fetal). Panel C: v0.3 clock variance from the ISS data; q = BH-adjusted within-stratum permutation test.",
         fontsize=7.3, color=GREY, wrap=True)
fig.savefig("results/real/fig_main_results.png", dpi=200, bbox_inches="tight")
print("ok")
