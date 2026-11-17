"""Synthetic datasets with KNOWN ground truth, used for tests and the demo.

They are not biological results. They exist to prove the pipeline recovers effects
that were planted on purpose (and does not invent effects that were not).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genesets import gene_lists

N_BACKGROUND = 3000


def _base(rng, n_genes):
    mu = rng.lognormal(mean=3.0, sigma=1.6, size=n_genes)
    return mu


def _genes():
    sets = gene_lists()
    named = sorted({g for v in sets.values() for g in v})
    return named + [f"BG{i:05d}" for i in range(N_BACKGROUND)], sets


def _draw(rng, mu_matrix):
    # negative binomial counts, dispersion 0.1
    disp = 0.1
    p = 1 / (1 + disp * mu_matrix)
    return rng.negative_binomial(1 / disp, p)


def timecourse(days=(15, 30, 60, 90, 120), reps=4, seed=0, maturation_slope=1.2):
    """Organoids of known age. Maturation & DA genes rise with log(day); progenitor falls."""
    rng = np.random.default_rng(seed)
    genes, sets = _genes()
    mu = pd.Series(_base(rng, len(genes)), index=genes)
    cols, ages, mats = [], [], []
    for d in days:
        t = np.log(d / days[0]) / np.log(days[-1] / days[0])  # 0..1
        for r in range(reps):
            m = mu.copy()
            m[sets["maturation"]] *= np.exp(maturation_slope * 2 * t)
            m[sets["da_identity"]] *= np.exp(1.5 * t)
            m[sets["progenitor_cycle"]] *= np.exp(-1.5 * t)
            m *= rng.lognormal(0, 0.15, len(m))
            cols.append(f"d{d}_r{r}")
            ages.append(d)
            mats.append(m.values * rng.uniform(0.5, 2.0))
    counts = pd.DataFrame(_draw(rng, np.array(mats).T), index=genes, columns=cols)
    return counts, pd.Series(ages, index=cols, name="day")


def fetal_reference(pcw=(6, 7, 8, 9, 10, 11, 13, 15, 17, 20), seed=1):
    """Synthetic fetal reference whose 'age genes' drift smoothly with PCW."""
    rng = np.random.default_rng(seed)
    genes, sets = _genes()
    mu = pd.Series(_base(rng, len(genes)), index=genes)
    slopes = pd.Series(rng.normal(0, 1.0, len(genes)), index=genes)
    slopes[rng.random(len(genes)) > 0.15] = 0  # 15% of genes are age-informative
    cols, ages, mats = [], [], []
    for w in pcw:
        t = (w - pcw[0]) / (pcw[-1] - pcw[0])
        for r in range(2):
            m = mu * np.exp(2.0 * slopes * t) * rng.lognormal(0, 0.1, len(mu))
            cols.append(f"PCW{w}_r{r}")
            ages.append(w)
            mats.append(m.values)
    counts = pd.DataFrame(_draw(rng, np.array(mats).T), index=genes, columns=cols)
    return counts, pd.Series(ages, index=cols, name="pcw"), mu, slopes, pcw


def query_at_pcw(mu, slopes, pcw_range, target_pcw, n=3, seed=2):
    rng = np.random.default_rng(seed)
    t = (target_pcw - pcw_range[0]) / (pcw_range[-1] - pcw_range[0])
    mats = [(mu * np.exp(2.0 * slopes * t) * rng.lognormal(0, 0.1, len(mu))).values for _ in range(n)]
    cols = [f"q{target_pcw}_{i}" for i in range(n)]
    return pd.DataFrame(_draw(rng, np.array(mats).T), index=mu.index, columns=cols)


def iss_like(seed=3, maturation_effect=0.8, hypoxia_effect=0.0, hypoxia_drives_maturation=0.0,
             identity_effect=0.0, n_per_cell=3):
    """Mimics CRS-19 design: 2 DA lines x microglia(+/-) x {ground, LEO}.

    hypoxia_drives_maturation > 0 makes part of the apparent maturation come from a
    per-organoid hypoxia level (the D hypothesis), so the adjusted test can be checked.
    """
    rng = np.random.default_rng(seed)
    genes, sets = _genes()
    mu = pd.Series(_base(rng, len(genes)), index=genes)
    line_shift = {L: pd.Series(rng.normal(0, 0.3, len(genes)), index=genes) for L in ("PD", "CTRL")}
    rows, mats = [], []
    i = 0
    for line in ("PD", "CTRL"):
        for mg in ("noMG", "MG"):
            for cond in ("ground", "LEO"):
                for _ in range(n_per_cell):
                    m = mu * np.exp(line_shift[line])
                    leo = cond == "LEO"
                    hyp = rng.normal(hypoxia_effect if leo else 0.0, 0.4)
                    m[sets["hypoxia"]] *= np.exp(hyp)
                    m[sets["maturation"]] *= np.exp((maturation_effect if leo else 0.0)
                                                    + hypoxia_drives_maturation * hyp)
                    m[sets["progenitor_cycle"]] *= np.exp(-(maturation_effect if leo else 0.0)
                                                          - hypoxia_drives_maturation * hyp)
                    m[sets["da_identity"]] *= np.exp(identity_effect if leo else 0.0)
                    if mg == "MG":
                        m[sets["microglia_qc"]] *= 20
                    m = m * rng.lognormal(0, 0.15, len(m))
                    name = f"S{i:02d}"
                    rows.append({"sample": name, "line": line, "microglia": mg, "condition": cond})
                    mats.append(m.values * rng.uniform(0.5, 2.0))
                    i += 1
    counts = pd.DataFrame(_draw(rng, np.array(mats).T), index=genes,
                          columns=[r["sample"] for r in rows])
    meta = pd.DataFrame(rows).set_index("sample")
    return counts, meta
