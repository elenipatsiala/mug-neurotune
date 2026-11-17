"""Validated maturity panel (v0.3): three read-outs that each tracked culture age in three independent
midbrain-organoid time courses from three labs (Lund GSE168323, BrainSTEM GSE281535, held-out GSE277032):

  age_clock        neuronal age clock v0.3 (fetal-trained consensus genes)      rho 0.74 / 0.99 / 0.95
  fetal_corr_age   whole-transcriptome correlation to fetal VM tissue (PCW)       rho 0.93 / 0.84 / 0.84
  progenitor_loss  minus the background-ranked expression of SOX2 + NES          rho 0.94 / 0.98 / 0.95

Commonly used 'maturity' markers (RBFOX3/NeuN, SYP, TH, MAP2) failed in 2 of 3 datasets because they
track how many neurons an organoid contains, not how old it is (see docs/BENCHMARK.md).

A maturation change between conditions is called only when the read-outs AGREE
(see compare.classify): at least 2 of 3 significantly shifted in the same direction and none
significantly shifted the other way.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import anchor, clock, scoring

_DATA = os.path.join(os.path.dirname(__file__), "data")
PANEL = ("age_clock", "fetal_corr_age", "progenitor_loss")


def progenitor_loss(logcpm: pd.DataFrame) -> pd.Series:
    if not {"SOX2", "NES"} <= set(logcpm.index):
        return pd.Series(np.nan, index=logcpm.columns)
    return -scoring.rank_score(logcpm, ["SOX2", "NES"])


def load_fetal_reference():
    ref = pd.read_csv(os.path.join(_DATA, "fetal_vm_reference.csv.gz"), index_col=0)
    age = pd.read_csv(os.path.join(_DATA, "fetal_vm_reference_age.csv"), index_col=0).iloc[:, 0]
    from .io import log_cpm
    return log_cpm(ref), age


def fetal_corr_age(logcpm: pd.DataFrame, n_genes: int = 300) -> pd.Series:
    ref, age = load_fetal_reference()
    shared = logcpm.index.intersection(ref.index)
    if len(shared) < 2000:
        return pd.Series(np.nan, index=logcpm.columns)
    expressed = shared[(logcpm.loc[shared] > 1).mean(axis=1) > 0.5]
    genes = anchor.age_informative_genes(ref.loc[expressed], age, n_genes=n_genes)
    return anchor.estimate_age(logcpm, ref, age, genes=genes)["est_pcw"]


def maturity_panel(logcpm: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({"age_clock": clock.clock_score(logcpm)["age_clock"],
                         "fetal_corr_age": fetal_corr_age(logcpm),
                         "progenitor_loss": progenitor_loss(logcpm)})
