"""Step B - anchor organoids to human fetal ventral-midbrain development.

Input: a reference matrix (genes x reference samples, log-CPM of fetal VM pseudobulk)
with known age in post-conception weeks (PCW).

Method (deliberately simple and explainable):
 1. Age-informative genes = genes whose reference expression correlates with PCW
    (|Spearman rho| largest), top `n_genes`.
 2. For each query sample, Spearman correlation with every reference sample
    over those genes (rank-based -> robust to bulk vs pseudobulk).
 3. Estimated PCW = softmax(correlation * temperature)-weighted mean of reference ages.
 4. 'Reference similarity' = max correlation. Low similarity means the organoid
    does not resemble ANY fetal midbrain stage (a fidelity warning).

Published cross-check: in an integrated midbrain atlas, day-40 organoids aligned to
~11 PCW and day-70/120 organoids to ~15 PCW (Nat Commun 2025, s41467-025-67779-1).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def age_informative_genes(ref: pd.DataFrame, ref_age: pd.Series, n_genes: int = 300,
                          exclude: set[str] | None = None) -> list[str]:
    exclude = exclude or set()
    a = ref_age.loc[ref.columns].values
    rows = []
    for g, v in ref.iterrows():
        if g in exclude or np.std(v.values) == 0:
            continue
        rows.append((g, stats.spearmanr(v.values, a)[0]))
    s = pd.Series(dict(rows)).dropna()
    return list(s.abs().sort_values(ascending=False).index[:n_genes])


def estimate_age(query: pd.DataFrame, ref: pd.DataFrame, ref_age: pd.Series,
                 genes: list[str] | None = None, temperature: float = 25.0) -> pd.DataFrame:
    genes = genes or age_informative_genes(ref, ref_age)
    genes = [g for g in genes if g in query.index and g in ref.index]
    if len(genes) < 20:
        raise ValueError(f"only {len(genes)} shared age-informative genes; check gene IDs")
    q = query.loc[genes].rank(axis=0)
    r = ref.loc[genes].rank(axis=0)
    corr = pd.DataFrame(np.corrcoef(q.T.values, r.T.values)[: q.shape[1], q.shape[1]:],
                        index=query.columns, columns=ref.columns)
    ages = ref_age.loc[ref.columns].values.astype(float)
    out = []
    for s, row in corr.iterrows():
        w = np.exp(temperature * (row.values - row.values.max()))
        w /= w.sum()
        out.append({"sample": s, "est_pcw": float((w * ages).sum()),
                    "best_match_pcw": float(ages[row.values.argmax()]),
                    "ref_similarity": float(row.values.max()), "n_genes": len(genes)})
    return pd.DataFrame(out).set_index("sample")
