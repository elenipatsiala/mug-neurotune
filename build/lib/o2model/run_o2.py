"""Run the vial O2 Monte Carlo, sensitivity analysis and figures."""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from o2model import vial_o2 as m  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/o2"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
C1, CU = "#64748b", "#1d4ed8"

CSV = os.path.join(OUT, "o2_montecarlo.csv")
if os.path.exists(CSV) and "--rerun" not in sys.argv:
    df = pd.read_csv(CSV)
else:
    df = m.monte_carlo(N, seed=42)
    df.to_csv(CSV, index=False)
LABELS = {"U_1g_um_s": "convection speed in ground vial", "u_pos": "µg: organoid–bubble distance",
          "contact_1g": "1 g: surface not touching wall", "ocr_amol_s": "O₂ use per cell",
          "c_mix": "vial mixing efficiency", "n_cells": "cell number", "C_far_frac": "headspace O₂ level",
          "radius_um": "organoid radius", "km_uM": "Km for O₂", "D_med": "O₂ diffusivity (medium)",
          "D_tissue_frac": "O₂ diffusivity (tissue)", "U_ug_um_s": "µg residual mixing",
          "vial_radius_mm": "vial radius"}

anyhyp_1g = df.hypoxic_frac_1pct_1g > 0.05
anyhyp_ug = df.hypoxic_frac_1pct_ug > 0.05
summary = {
    "n": N,
    "P(ug more hypoxic than 1g, by >=5 pct-points of volume)": float((df.delta_hypoxic > 0.05).mean()),
    "P(1g more hypoxic than ug, by >=5 pct-points)": float((df.delta_hypoxic < -0.05).mean()),
    "P(similar within 5 pct-points)": float((df.delta_hypoxic.abs() <= 0.05).mean()),
    "P(ug surface O2 lower than 1g)": float((df.delta_surface_pct < 0).mean()),
    "median surface O2 % 1g": float(df.surface_pct_o2_1g.median()),
    "median surface O2 % ug": float(df.surface_pct_o2_ug.median()),
    "median hypoxic fraction 1g": float(df.hypoxic_frac_1pct_1g.median()),
    "median hypoxic fraction ug": float(df.hypoxic_frac_1pct_ug.median()),
    "P(>5% of organoid hypoxic) 1g": float(anyhyp_1g.mean()),
    "P(>5% of organoid hypoxic) ug": float(anyhyp_ug.mean()),
    "median R_far/R_near 1g": float((df.R_far_1g / df.R_near_1g).median()),
    "median R_far/R_near ug": float((df.R_far_ug / df.R_near_ug).median()),
    "sealed-vial O2 budget days (median, 5th-95th pct)":
        [float(df.budget_days.median()), float(df.budget_days.quantile(.05)), float(df.budget_days.quantile(.95))],
    "P(sealed impermeable vial runs out of O2 before day 30)": float((df.budget_days < 30).mean()),
}
# global sensitivity: Spearman of each input with delta_hypoxic
inputs = list(m.PRIORS) + ["u_pos"]
sens = pd.Series({k: stats.spearmanr(df[k], df.delta_hypoxic)[0] for k in inputs}).sort_values(key=abs, ascending=False)
sens.to_csv(os.path.join(OUT, "sensitivity_delta_hypoxic.csv"), header=["spearman_rho"])
summary["top_drivers_of_difference"] = sens.head(5).round(3).to_dict()
json.dump(summary, open(os.path.join(OUT, "o2_summary.json"), "w"), indent=2)
print(json.dumps(summary, indent=2))

# --------------------------------------------------------------- figure 1: paired outcome
fig, axes = plt.subplots(1, 3, figsize=(14, 4), gridspec_kw={"width_ratios": [1, 1, 1.25]})
ax = axes[0]
bins = np.linspace(0, max(df.surface_pct_o2_1g.max(), df.surface_pct_o2_ug.max()), 40)
ax.hist(df.surface_pct_o2_1g, bins, color=C1, alpha=.7, label="1 g ground")
ax.hist(df.surface_pct_o2_ug, bins, color=CU, alpha=.6, label="µg ISS")
ax.set_xlabel("O₂ at organoid surface (% O₂ equiv.)")
ax.set_ylabel("Monte Carlo runs")
ax.set_title("A  Surface O₂", loc="left", fontweight="bold")
ax.legend(frameon=False)
ax = axes[1]
ax.scatter(df.hypoxic_frac_1pct_1g, df.hypoxic_frac_1pct_ug, s=5, c=CU, alpha=.35)
ax.plot([0, 1], [0, 1], color="#94a3b8", lw=1, ls="--")
ax.set_xlabel("hypoxic fraction (<1 % O₂), 1 g")
ax.set_ylabel("hypoxic fraction (<1 % O₂), µg")
ax.set_title("B  Same organoid, two gravities", loc="left", fontweight="bold")
ax.text(.03, .93, f"above line: µg worse ({(df.delta_hypoxic > 0.05).mean():.0%})", transform=ax.transAxes, fontsize=8)
ax = axes[2]
top = sens.head(7)[::-1]
ax.barh([LABELS.get(k, k) for k in top.index], top.values, color=[CU if v > 0 else C1 for v in top.values])
ax.axvline(0, color="k", lw=.6)
ax.set_xlabel("ρ with (µg − 1 g) hypoxic fraction")
ax.set_title("C  Drivers of the difference", loc="left", fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_o2_montecarlo.png"), dpi=200)
plt.close(fig)

# --------------------------------------------------------------- figure 2: regime map
base = pd.Series({k: (lo * hi) ** .5 if sc == "log" else (lo + hi) / 2 for k, (lo, hi, sc) in m.PRIORS.items()})
base["u_pos"] = 0.5
radii = np.linspace(100, 250, 31)
ocrs = np.exp(np.linspace(np.log(5), np.log(100), 31))
grid = {"1g": np.zeros((31, 31)), "ug": np.zeros((31, 31))}
for i, oc in enumerate(ocrs):
    for j, ra in enumerate(radii):
        p = base.copy()
        p["ocr_amol_s"], p["radius_um"] = oc, ra
        for g in grid:
            grid[g][i, j] = m.run_condition(p, g)["hypoxic_frac_1pct"]
np.savez(os.path.join(OUT, "regime_grid.npz"), radii=radii, ocrs=ocrs, **grid)
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
for ax, (title, Z, cmap, vmin, vmax) in zip(axes, [
        ("1 g: hypoxic fraction", grid["1g"], "Greys", 0, 1),
        ("µg: hypoxic fraction", grid["ug"], "Blues", 0, 1),
        ("µg − 1 g", grid["ug"] - grid["1g"], "RdBu_r", -.5, .5)]):
    im = ax.pcolormesh(radii, ocrs, Z, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    ax.set_yscale("log")
    ax.set_xlabel("organoid radius (µm)")
    ax.set_title(title, loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax)
axes[0].set_ylabel("O₂ consumption (amol/cell/s)")
fig.suptitle("Regime map, other parameters at median (cell number FIXED at ~1e5, so a larger organoid is less dense and less hypoxic)", x=.01, ha="left", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_o2_regime_map.png"), dpi=200)
plt.close(fig)
print("figures written to", OUT)
