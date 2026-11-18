"""Figures. Colour roles: ground = slate grey, LEO = deep blue; stress = orange ramp."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COND_COLORS = {"ground": "#64748b", "LEO": "#1d4ed8"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "DejaVu Sans"})


def maturation_fidelity_map(scores: pd.DataFrame, meta: pd.DataFrame, path: str,
                            cond_col="condition", stress_col="hypoxia", title=None):
    fig, ax = plt.subplots(figsize=(6.2, 5))
    s = scores.loc[meta.index].copy()
    xcol = "age_clock" if "age_clock" in s.columns else "maturation"
    s["_x"] = s[xcol]
    stress = s[stress_col]
    norm = plt.Normalize(stress.min(), stress.max())
    markers = {}
    for i, line in enumerate(sorted(meta.get("line", pd.Series("all", index=meta.index)).astype(str).unique())):
        markers[line] = "os^Dv"[i % 5]
    for (cond, line), sub in meta.groupby([meta[cond_col].astype(str),
                                          meta.get("line", pd.Series("all", index=meta.index)).astype(str)]):
        ax.scatter(s.loc[sub.index, "_x"], s.loc[sub.index, "fidelity"],
                   c=stress.loc[sub.index], cmap="Oranges", norm=norm, s=70,
                   marker=markers[line], edgecolors=COND_COLORS.get(cond, "k"), linewidths=2,
                   label=f"{cond} · {line}")
    for cond, sub in meta.groupby(meta[cond_col].astype(str)):
        mx, my = s.loc[sub.index, "_x"].mean(), s.loc[sub.index, "fidelity"].mean()
        ax.scatter(mx, my, marker="X", s=180, color=COND_COLORS.get(cond, "k"), zorder=5)
        ax.annotate(f"{cond} mean", (mx, my), xytext=(6, 6), textcoords="offset points",
                    color=COND_COLORS.get(cond, "k"), fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap="Oranges", norm=norm)
    fig.colorbar(sm, ax=ax, label=f"{stress_col} score (fill)")
    ax.set_xlabel("Neuronal age clock v0.3  →  older" if xcol == "age_clock" else "Synaptic marker score (v0.1, not validated)")
    ax.set_ylabel("Dopaminergic fidelity (DA identity − off-target)")
    ax.set_title(title or "µG-NeuroTune Age–Fidelity map", loc="left", fontweight="bold")
    ax.legend(fontsize=7, frameon=False, loc="best", title="edge = condition")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def validation_plot(scores: pd.DataFrame, age: pd.Series, res: pd.DataFrame, path: str):
    progs = [p for p in ("maturation", "da_identity", "progenitor_cycle", "hypoxia") if p in scores]
    fig, axes = plt.subplots(1, len(progs), figsize=(3.1 * len(progs), 3.2), sharex=True)
    for ax, p in zip(np.atleast_1d(axes), progs):
        ax.scatter(age.loc[scores.index], scores[p], color="#1d4ed8", s=28)
        r = res.loc[p]
        ax.set_title(f"{p}\nρ = {r.rho:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]", fontsize=9, loc="left")
        ax.set_xlabel("days in culture")
    np.atleast_1d(axes)[0].set_ylabel("score")
    verdict = "PASS" if res.attrs.get("maturation_validated") else "FAIL"
    fig.suptitle(f"Step A: maturation score validation on known-age organoids — {verdict}",
                 x=0.01, ha="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def power_plot(grid: pd.DataFrame, path: str, title: str):
    ratios = sorted(grid.dxc_ratio.unique())
    fig, axes = plt.subplots(1, len(ratios), figsize=(3.4 * len(ratios), 3.3), sharey=True)
    shades = ["#bfdbfe", "#60a5fa", "#2563eb", "#1e3a8a"]
    for ax, r in zip(np.atleast_1d(axes), ratios):
        sub = grid[grid.dxc_ratio == r]
        for i, (no, g) in enumerate(sub.groupby("organoids_per_donor_per_arm")):
            ax.plot(g.n_donors, g.power, "-o", color=shades[i % 4], label=f"{no} organoids/donor/arm", ms=4)
        ax.axhline(0.8, color="#94a3b8", ls="--", lw=1)
        ax.set_title(f"donor×LEO SD = {r:g} × organoid SD", fontsize=9, loc="left")
        ax.set_xlabel("donors (cell lines)")
    np.atleast_1d(axes)[0].set_ylabel("power (α = 0.05)")
    np.atleast_1d(axes)[-1].legend(fontsize=7, frameon=False)
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
