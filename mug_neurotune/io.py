"""Loading counts (bulk or pseudobulk), gene-ID mapping and normalisation."""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd


def _read_table(path: str) -> pd.DataFrame:
    sep = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=sep, index_col=0, comment="#")


def load_counts(path: str, annotation: str | None = None,
                symbol_col: str = "Symbol") -> pd.DataFrame:
    """Load a genes x samples raw-count matrix.

    * Rows may be gene symbols, Ensembl IDs or NCBI GeneIDs.
    * If `annotation` is given (e.g. GEO's `Human.GRCh38.p13.annot.tsv.gz`, or a
      BioMart export) the row IDs are mapped to `symbol_col`.
    * featureCounts output (Chr/Start/End/Strand/Length columns) is handled.
    Duplicate symbols are summed.
    """
    df = _read_table(path)
    meta_cols = [c for c in ("Chr", "Start", "End", "Strand", "Length", "gene_name",
                             "Symbol", "symbol") if c in df.columns]
    if annotation is None and any(c in df.columns for c in ("gene_name", "Symbol", "symbol")):
        col = next(c for c in ("gene_name", "Symbol", "symbol") if c in df.columns)
        df.index = df[col].astype(str)
    df = df.drop(columns=meta_cols, errors="ignore")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if annotation is not None:
        ann = _read_table(annotation)
        if symbol_col not in ann.columns:
            raise ValueError(f"annotation has no column {symbol_col!r}; columns: {list(ann.columns)[:10]}")
        ann.index = ann.index.astype(str)
        df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)  # strip Ensembl version
        mapping = ann[symbol_col].astype(str)
        df = df[df.index.isin(mapping.index)]
        df.index = mapping.reindex(df.index).values
    df.index = df.index.astype(str)
    df = df.groupby(level=0).sum()
    return df


def load_metadata(path: str, samples: list[str] | None = None) -> pd.DataFrame:
    meta = _read_table(path)
    meta.index = meta.index.astype(str)
    if samples is not None:
        missing = set(samples) - set(meta.index)
        if missing:
            raise ValueError(f"metadata lacks {len(missing)} samples, e.g. {sorted(missing)[:3]}")
        meta = meta.loc[samples]
    return meta


def log_cpm(counts: pd.DataFrame, prior: float = 1.0) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + prior)


def filter_expressed(logcpm: pd.DataFrame, min_logcpm: float = 1.0,
                     min_frac: float = 0.2) -> pd.DataFrame:
    keep = (logcpm > min_logcpm).mean(axis=1) >= min_frac
    return logcpm[keep]


def pseudobulk_10x(mtx_dir: str, cell_groups: pd.Series) -> pd.DataFrame:
    """Sum single-cell counts into one profile per sample (pseudobulk).

    mtx_dir: 10x folder with matrix.mtx(.gz), features.tsv(.gz)/genes.tsv, barcodes.tsv(.gz)
    cell_groups: Series indexed by barcode -> sample/organoid id.
    Pseudobulk makes scRNA time-courses comparable with bulk RNA-seq (ISS data).
    """
    from scipy.io import mmread

    def _open(name):
        for cand in (name, name + ".gz"):
            p = os.path.join(mtx_dir, cand)
            if os.path.exists(p):
                return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
        raise FileNotFoundError(name)

    with _open("matrix.mtx") as fh:
        m = mmread(fh).tocsr()  # genes x cells
    try:
        fh = _open("features.tsv")
    except FileNotFoundError:
        fh = _open("genes.tsv")
    feats = pd.read_csv(fh, sep="\t", header=None)
    symbols = feats[1] if feats.shape[1] > 1 else feats[0]
    barcodes = pd.read_csv(_open("barcodes.tsv"), header=None)[0].astype(str)
    groups = cell_groups.reindex(barcodes.values)
    out = {}
    for g in groups.dropna().unique():
        idx = np.where(groups.values == g)[0]
        out[g] = np.asarray(m[:, idx].sum(axis=1)).ravel()
    df = pd.DataFrame(out, index=symbols.astype(str).values)
    return df.groupby(level=0).sum()
