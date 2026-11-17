"""LEO vs ground comparison, hypoxia-adjusted test (step D) and evidence classification.

Design-aware: comparisons are made WITHIN strata (cell line x microglia status), so a
difference between donors can never masquerade as a LEO effect. P-values come from
permuting the condition label within strata.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genesets import DYSREGULATION

# Pre-registered classification thresholds (Hedges' g, pooled across strata).
# v0.1.0 rule set: |g| >= 0.5, bootstrap CI excludes 0, BH q < 0.05, and the same
# direction in every cell line (added 2026-11-17 before any real data, after
# simulation showed chance single-line hits; see docs/CHANGELOG.md).
G_MEANINGFUL = 0.5
Q_MAX = 0.05


def _hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return np.nan
    j = 1 - 3 / (4 * (na + nb) - 9)
    return j * (b.mean() - a.mean()) / sp


def stratified_effect(values: pd.Series, cond: pd.Series, strata: pd.Series,
                      ref: str, test: str) -> float:
    gs, ws = [], []
    for _, idx in strata.groupby(strata).groups.items():
        v, c = values.loc[idx], cond.loc[idx]
        a, b = v[c == ref].values, v[c == test].values
        g = _hedges_g(a, b)
        if not np.isnan(g):
            gs.append(g)
            ws.append(len(a) + len(b))
    return float(np.average(gs, weights=ws)) if gs else np.nan


def stratified_diff(values: pd.Series, cond: pd.Series, strata: pd.Series,
                    ref: str, test: str) -> float:
    """Size-weighted mean of within-stratum (test - ref) differences, in score units."""
    ds, ws = [], []
    for _, idx in strata.groupby(strata).groups.items():
        v, c = values.loc[idx], cond.loc[idx]
        a, b = v[c == ref].values, v[c == test].values
        if len(a) and len(b):
            ds.append(b.mean() - a.mean())
            ws.append(len(a) + len(b))
    return float(np.average(ds, weights=ws)) if ds else np.nan


def _permute_within(cond: pd.Series, strata: pd.Series, rng) -> pd.Series:
    out = cond.copy()
    for _, idx in strata.groupby(strata).groups.items():
        out.loc[idx] = rng.permutation(cond.loc[idx].values)
    return out


def _boot_within(index: pd.Index, cond: pd.Series, strata: pd.Series, rng) -> list:
    pick = []
    for _, idx in strata.groupby(strata).groups.items():
        for c in cond.loc[idx].unique():
            members = [i for i in idx if cond.loc[i] == c]
            pick += list(rng.choice(members, len(members), replace=True))
    return pick


def compare_programs(scores: pd.DataFrame, meta: pd.DataFrame, cond_col: str = "condition",
                     ref: str = "ground", test: str = "LEO", strata_cols=("line", "microglia"),
                     n_perm: int = 5000, n_boot: int = 2000, seed: int = 1) -> pd.DataFrame:
    """Stratified LEO-vs-ground effects for every program (vectorised over permutations)."""
    rng = np.random.default_rng(seed)
    meta = meta.loc[scores.index]
    cond = meta[cond_col].astype(str)
    keep = cond.isin([ref, test]).values
    scores, meta, cond = scores[keep], meta[keep], cond[keep]
    strata = (meta[list(strata_cols)].astype(str).agg("|".join, axis=1) if strata_cols
              else pd.Series("all", index=meta.index))
    lines = meta["line"].astype(str).values if "line" in meta.columns else None
    is_test = (cond.values == test)
    groups = [np.where(strata.values == k)[0] for k in sorted(strata.unique())]
    groups = [g for g in groups if is_test[g].any() and (~is_test[g]).any()]
    n = len(cond)

    def perm_masks(k):
        M = np.tile(is_test, (k, 1))
        for g in groups:
            M[:, g] = rng.permuted(M[:, g], axis=1)
        return M

    def eff(V, M):
        """V: (k, n) values, M: (k, n) bool test mask -> (diff, hedges g) arrays of len k."""
        d_sum = np.zeros(V.shape[0]); w_sum = 0.0
        g_sum = np.zeros(V.shape[0]); wg = 0.0
        for g in groups:
            v, m = V[:, g], M[:, g]
            nb = m.sum(1); na = (~m).sum(1)
            mb = (v * m).sum(1) / nb
            ma = (v * ~m).sum(1) / na
            ntot = len(g)
            d_sum += (mb - ma) * ntot; w_sum += ntot
            if na.min() > 1 and nb.min() > 1:
                vb = ((v - mb[:, None]) ** 2 * m).sum(1) / (nb - 1)
                va = ((v - ma[:, None]) ** 2 * ~m).sum(1) / (na - 1)
                sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (ntot - 2))
                J = 1 - 3 / (4 * ntot - 9)
                with np.errstate(divide="ignore", invalid="ignore"):
                    g_sum += np.where(sp > 0, J * (mb - ma) / sp, np.nan) * ntot
                wg += ntot
        return d_sum / w_sum, (g_sum / wg if wg else np.full(V.shape[0], np.nan))

    Mperm = perm_masks(n_perm)
    # stratified bootstrap: resample within each stratum x condition cell
    boot_idx = np.empty((n_boot, n), dtype=int)
    for g in groups:
        for flag in (True, False):
            members = g[is_test[g] == flag]
            boot_idx[:, members] = rng.choice(members, (n_boot, len(members)), replace=True)
    rows = []
    for prog in scores.columns:
        v = scores[prog].values.astype(float)
        diff, g = (x[0] for x in eff(v[None, :], is_test[None, :]))
        _, g_null = eff(np.tile(v, (n_perm, 1)), Mperm)
        ok = ~np.isnan(g_null)
        p = (np.sum(np.abs(g_null[ok]) >= abs(g) - 1e-12) + 1) / (ok.sum() + 1)
        d_boot, _ = eff(v[boot_idx], np.tile(is_test, (n_boot, 1)))
        lo, hi = np.nanpercentile(d_boot, [2.5, 97.5])
        consistent = True
        if lines is not None:
            for L in np.unique(lines):
                sel = lines == L
                dl = []
                for gi in groups:
                    gg = gi[sel[gi]]
                    if is_test[gg].any() and (~is_test[gg]).any():
                        dl.append(v[gg][is_test[gg]].mean() - v[gg][~is_test[gg]].mean())
                if dl and np.sign(np.mean(dl)) != np.sign(diff):
                    consistent = False
        rows.append({"program": prog, "hedges_g": g, "diff": diff, "ci_low": lo, "ci_high": hi,
                     "p_perm": p, "consistent_lines": consistent})
    res = pd.DataFrame(rows).set_index("program")
    pv = res["p_perm"].values
    order = np.argsort(pv)[::-1]
    q = np.empty_like(pv)
    prev = 1.0
    for rank, i in enumerate(order):
        k = len(pv) - rank
        prev = min(prev, pv[i] * len(pv) / k)
        q[i] = prev
    res["q_bh"] = q
    res.attrs["n_strata"] = len(groups)
    return res


def hypoxia_adjusted(scores: pd.DataFrame, meta: pd.DataFrame, target: str = "maturation",
                     covariate: str = "hypoxia", cond_col: str = "condition", ref: str = "ground",
                     test: str = "LEO", strata_cols=("line", "microglia")) -> dict:
    """Does the LEO maturation effect survive removing what hypoxia explains?

    Within each stratum, regress target on covariate using ALL samples (pooled
    within-stratum slope), take residuals, and recompute the LEO effect.
    share_explained = 1 - g_adjusted / g_raw.
    """
    meta = meta.loc[scores.index]
    cond = meta[cond_col].astype(str)
    strata = meta[list(strata_cols)].astype(str).agg("|".join, axis=1)
    y, x = scores[target], scores[covariate]
    yc = y - y.groupby(strata).transform("mean")
    xc = x - x.groupby(strata).transform("mean")
    slope = float((xc * yc).sum() / (xc ** 2).sum()) if (xc ** 2).sum() > 0 else 0.0
    resid = y - slope * x
    g_raw = stratified_effect(y, cond, strata, ref, test)
    g_adj = stratified_effect(resid, cond, strata, ref, test)
    r = float(np.corrcoef(xc, yc)[0, 1]) if xc.std() > 0 and yc.std() > 0 else np.nan
    share = 1 - g_adj / g_raw if g_raw and not np.isnan(g_raw) and g_raw != 0 else np.nan
    return {"target": target, "covariate": covariate, "within_strata_r": r,
            "slope": slope, "g_raw": g_raw, "g_adjusted": g_adj,
            "share_explained": share,
            "verdict": ("not assessable" if not np.isfinite(share)
                        else "hypoxia explains most of the effect" if share > 0.5
                        else "effect largely independent of hypoxia score")}


def classify(res: pd.DataFrame) -> tuple[str, str]:
    """Pre-registered evidence classification (see docs/PREREGISTRATION.md)."""
    # A program counts as changed only if: |g| >= 0.5 AND the bootstrap CI of the
    # difference excludes 0 AND Benjamini-Hochberg q < Q_MAX across all programs.
    def up(p):
        r = res.loc[p] if p in res.index else None
        return (r is not None and r.hedges_g >= G_MEANINGFUL and r.ci_low > 0 and r.q_bh < Q_MAX
                and bool(r.get("consistent_lines", True)))

    def down(p):
        r = res.loc[p] if p in res.index else None
        return (r is not None and r.hedges_g <= -G_MEANINGFUL and r.ci_high < 0 and r.q_bh < Q_MAX
                and bool(r.get("consistent_lines", True)))

    mat_up = up("maturation")
    fid_down = down("da_identity") or down("fidelity")
    dys_up = [p for p in DYSREGULATION if up(p)]
    if not mat_up:
        return "insufficient evidence", "Maturation increase is not robust (g < 0.5 or CI includes 0)."
    if fid_down:
        return "dysregulated state", "Maturation up but dopaminergic identity/fidelity down."
    if dys_up:
        return "mixed state", "Maturation up together with: " + ", ".join(dys_up)
    return "promising maturation", "Maturation up, identity preserved, no dysregulation program up."
