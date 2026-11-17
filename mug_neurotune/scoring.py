"""Single-sample gene-program scores.

Primary method: a rank-based score (in the spirit of singscore, Foroutan et al. 2018).
Within each sample, expressed genes are ranked; the score of a program is the mean
normalised rank of its genes minus 0.5, so it lies in [-0.5, 0.5]. Rank scores do
not depend on library size or platform scaling, which is what lets one score work
on bulk RNA-seq (ISS) and on scRNA pseudobulk (time-courses, fetal references).

Sensitivity method: mean of per-gene z-scores across samples (dataset-relative).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genesets import GENESETS, gene_lists


def _background(logcpm: pd.DataFrame) -> pd.Index:
    program_genes = {g for v in gene_lists().values() for g in v}
    bg = logcpm.index.difference(program_genes)
    return bg if len(bg) >= 100 else logcpm.index


def rank_score(logcpm: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean background-referenced rank of the program's genes, minus 0.5.

    Each gene is ranked against BACKGROUND genes only (genes in no program), never
    against other program genes. Otherwise a rise in one program (e.g. 18 maturation
    genes) pushes every other gene down the within-sample ranking and creates fake
    decreases in unrelated programs (caught by tests/test_pipeline.py::test_no_rank_leakage).
    """
    present = [g for g in genes if g in logcpm.index]
    if not present:
        return pd.Series(np.nan, index=logcpm.columns)
    bg = _background(logcpm)
    out = {}
    for s in logcpm.columns:
        ref = np.sort(logcpm.loc[bg, s].values)
        v = logcpm.loc[present, s].values
        r = (np.searchsorted(ref, v, "left") + np.searchsorted(ref, v, "right")) / 2 / len(ref)
        out[s] = r.mean() - 0.5
    return pd.Series(out)


def z_score(logcpm: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in logcpm.index]
    if not present:
        return pd.Series(np.nan, index=logcpm.columns)
    sub = logcpm.loc[present]
    sd = sub.std(axis=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
    return z.mean(axis=0)


METHODS = {"rank": rank_score, "z": z_score}


def coverage(logcpm: pd.DataFrame, sets: dict[str, list[str]] | None = None) -> pd.DataFrame:
    sets = sets or gene_lists()
    rows = []
    for name, genes in sets.items():
        present = [g for g in genes if g in logcpm.index]
        rows.append({"program": name, "n_genes": len(genes), "n_present": len(present),
                     "missing": ",".join(sorted(set(genes) - set(present)))})
    return pd.DataFrame(rows).set_index("program")


def score_all(logcpm: pd.DataFrame, method: str = "rank",
              sets: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Return samples x programs score table, plus derived 'fidelity' and (when the clock
    genes are present) the validated neuronal 'age_clock'."""
    default_sets = sets is None
    sets = sets or gene_lists()
    fn = METHODS[method]
    out = pd.DataFrame({name: fn(logcpm, genes) for name, genes in sets.items()})
    if "da_identity" in out and "off_target" in out:
        out["fidelity"] = out["da_identity"] - out["off_target"]
    if method == "rank" and default_sets:
        from .maturity import maturity_panel
        panel = maturity_panel(logcpm)
        for col in panel.columns:
            if panel[col].notna().all():
                out[col] = panel[col]
    return out


def leave_one_gene_out(logcpm: pd.DataFrame, program: str, group: pd.Series,
                       a: str, b: str, method: str = "rank") -> pd.DataFrame:
    """Effect (mean b - mean a) when each gene is dropped in turn.

    If one gene flips the sign, the program-level conclusion is fragile.
    """
    genes = [g for g in GENESETS[program]["genes"] if g in logcpm.index]
    fn = METHODS[method]
    rows = []
    full = fn(logcpm, genes)
    base = full[group == b].mean() - full[group == a].mean()
    for g in genes:
        s = fn(logcpm, [x for x in genes if x != g])
        eff = s[group == b].mean() - s[group == a].mean()
        rows.append({"dropped": g, "effect": eff, "sign_flip": np.sign(eff) != np.sign(base)})
    df = pd.DataFrame(rows)
    df.attrs["full_effect"] = base
    return df
