"""Neuronal age clock.

v0.3 (default): genes whose age trend in fetal ventral-midbrain NEURONS has the same sign in two
independent datasets (Birtele et al. 2022, GSE192405; La Manno et al. 2016, GSE76381). Frozen in git
before a held-out test; validated on three organoid time courses from three labs (rho 0.74-0.99).
v0.2 (kept for history): trained on GSE192405 only; failed the second organoid dataset (rho 0.49).

score = mean background-referenced rank of 50 age-UP genes
      - mean background-referenced rank of 50 age-DOWN genes
Training scripts: examples/train_fetal_clock.py (v0.2), examples/train_fetal_clock_v03.py (v0.3).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .genesets import gene_lists

_DATA = os.path.join(os.path.dirname(__file__), "data")
_PATH = os.path.join(_DATA, "neuronal_age_clock_v03.json")
_VERSIONS = {"0.3": "neuronal_age_clock_v03.json", "0.2": "neuronal_age_clock_v02.json"}
# HGNC renames between older 10x references and current NCBI/Ensembl annotations
ALIASES = {"ATP5MPL": "ATP5MJ", "ATP5MD": "ATP5MK", "H3F3A": "H3-3A", "ATP5ME": "ATP5ME",
           "C1orf61": "CROC4"}


def load_clock(version: str = "0.3") -> dict:
    with open(os.path.join(_DATA, _VERSIONS[version])) as fh:
        return json.load(fh)


def _resolve(genes, index):
    out = []
    for g in genes:
        if g in index:
            out.append(g)
        elif ALIASES.get(g) in index:
            out.append(ALIASES[g])
    return out


def clock_score(logcpm: pd.DataFrame) -> pd.DataFrame:
    ck = load_clock()
    up, down = _resolve(ck["up"], logcpm.index), _resolve(ck["down"], logcpm.index)
    if len(up) < 0.8 * len(ck["up"]) or len(down) < 0.8 * len(ck["down"]):
        out = pd.DataFrame(np.nan, index=logcpm.columns, columns=["age_clock", "clock_up", "clock_down"])
        out.attrs["coverage"] = {"up": f"{len(up)}/{len(ck['up'])}", "down": f"{len(down)}/{len(ck['down'])}"}
        return out
    exclude = set(up) | set(down) | {g for v in gene_lists().values() for g in v}
    bg = logcpm.index.difference(list(exclude))
    rows = {}
    for s in logcpm.columns:
        ref = np.sort(logcpm.loc[bg, s].values)

        def mr(genes):
            v = logcpm.loc[genes, s].values
            return ((np.searchsorted(ref, v, "left") + np.searchsorted(ref, v, "right")) / 2 / len(ref)).mean()
        rows[s] = {"age_clock": mr(up) - mr(down), "clock_up": mr(up) - 0.5, "clock_down": mr(down) - 0.5}
    out = pd.DataFrame(rows).T
    out.attrs["coverage"] = {"up": f"{len(up)}/{len(ck['up'])}", "down": f"{len(down)}/{len(ck['down'])}"}
    return out
