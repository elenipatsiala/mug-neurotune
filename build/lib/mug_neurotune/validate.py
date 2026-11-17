"""Step A - does the Maturation score actually track maturation?

Test on organoids of KNOWN age (e.g. Fiorenzano et al. 2021, days 15-120).
Pre-registered pass criterion (docs/PREREGISTRATION.md):
  Spearman rho(maturation score, culture day) >= 0.6 AND bootstrap 95% CI lower bound > 0
  AND progenitor_cycle score falls with age (rho < 0).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

PASS_RHO = 0.6


def spearman_ci(x, y, n_boot: int = 5000, seed: int = 0) -> tuple[float, float, float, float]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = []
    n = len(x)
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(np.unique(y[i])) < 2 or len(np.unique(x[i])) < 2:
            continue
        boots.append(stats.spearmanr(x[i], y[i])[0])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(rho), float(lo), float(hi), float(p)


def validate_timecourse(scores: pd.DataFrame, age: pd.Series,
                        programs=("maturation", "da_identity", "progenitor_cycle",
                                  "stress_glycolysis", "stress_er", "hypoxia")) -> pd.DataFrame:
    rows = []
    for prog in programs:
        if prog not in scores:
            continue
        rho, lo, hi, p = spearman_ci(scores[prog].values, age.loc[scores.index].values)
        rows.append({"program": prog, "rho": rho, "ci_low": lo, "ci_high": hi, "p": p})
    res = pd.DataFrame(rows).set_index("program")
    mat = res.loc["maturation"] if "maturation" in res.index else None
    prog = res.loc["progenitor_cycle"] if "progenitor_cycle" in res.index else None
    passed = bool(mat is not None and mat.rho >= PASS_RHO and mat.ci_low > 0
                  and (prog is None or prog.rho < 0))
    res.attrs["maturation_validated"] = passed
    return res
