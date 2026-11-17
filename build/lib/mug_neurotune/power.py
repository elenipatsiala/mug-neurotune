"""Step E - how big must the next flight be?

Model for one program score y:
    y = mu + (beta + b_d) * LEO + u_d + e
      u_d  ~ N(0, sd_donor)        donor (cell-line) baseline
      b_d  ~ N(0, sd_donor_x_cond) donor-specific LEO response  <- limits generalisation
      e    ~ N(0, sd_organoid)     organoid-to-organoid noise
Each donor contributes organoids to BOTH arms (as on CRS-19).

Analysis per simulated mission: per-donor mean difference (LEO - ground), then a
one-sample t-test across donors. This is the honest test if the claim is meant to
generalise to new donors. With 2 donors it has 1 degree of freedom.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def variance_components(scores: pd.Series, meta: pd.DataFrame, cond_col="condition",
                        donor_col="line", ref="ground", test="LEO") -> dict:
    """Estimate sd_organoid, sd_donor and sd_donor_x_cond from an existing dataset."""
    df = pd.DataFrame({"y": scores, "c": meta.loc[scores.index, cond_col].astype(str),
                       "d": meta.loc[scores.index, donor_col].astype(str)})
    cell = df.groupby(["d", "c"])["y"]
    sd_org = float(np.sqrt(cell.var(ddof=1).mean()))
    n_per = cell.size().mean()
    means = cell.mean().unstack()
    diff = (means[test] - means[ref]).dropna()
    beta = float(diff.mean())
    # var(diff_d) = sd_dxc^2 + 2*sd_org^2/n
    var_dxc = max(float(diff.var(ddof=1)) - 2 * sd_org ** 2 / n_per, 0.0) if len(diff) > 1 else np.nan
    base = means[ref].dropna()
    var_d = max(float(base.var(ddof=1)) - sd_org ** 2 / n_per, 0.0) if len(base) > 1 else np.nan
    return {"beta": beta, "sd_organoid": sd_org, "sd_donor": float(np.sqrt(var_d)),
            "sd_donor_x_cond": float(np.sqrt(var_dxc)), "n_donors_used": int(len(diff)),
            "note": "with <=2 donors, donor-level SDs are extremely uncertain; "
                    "power grid therefore sweeps sd_donor_x_cond as a ratio of sd_organoid"}


def simulate_power(beta: float, sd_organoid: float, sd_donor_x_cond: float,
                   n_donors: int, n_org: int, alpha: float = 0.05, n_sim: int = 4000,
                   seed: int = 0) -> float:
    if n_donors < 2:
        return np.nan
    rng = np.random.default_rng(seed)
    b = rng.normal(0, sd_donor_x_cond, (n_sim, n_donors))
    e_leo = rng.normal(0, sd_organoid, (n_sim, n_donors, n_org)).mean(axis=2)
    e_gnd = rng.normal(0, sd_organoid, (n_sim, n_donors, n_org)).mean(axis=2)
    d = beta + b + e_leo - e_gnd
    t = d.mean(axis=1) / (d.std(axis=1, ddof=1) / np.sqrt(n_donors))
    p = 2 * stats.t.sf(np.abs(t), df=n_donors - 1)
    return float(np.mean(p < alpha))


def power_grid(beta: float, sd_organoid: float, ratios=(0.0, 0.5, 1.0),
               donors=(2, 3, 4, 6, 8, 12), organoids=(2, 4, 8, 12)) -> pd.DataFrame:
    rows = []
    for r in ratios:
        for nd in donors:
            for no in organoids:
                rows.append({"dxc_ratio": r, "n_donors": nd, "organoids_per_donor_per_arm": no,
                             "total_organoids_flown": nd * no,
                             "power": simulate_power(beta, sd_organoid, r * sd_organoid, nd, no)})
    return pd.DataFrame(rows)


def minimal_design(grid: pd.DataFrame, target: float = 0.8) -> pd.DataFrame:
    ok = grid[grid.power >= target].sort_values(["dxc_ratio", "total_organoids_flown", "n_donors"])
    return ok.groupby("dxc_ratio").head(1)
