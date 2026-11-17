"""Command-line interface.

    mugnt score     COUNTS [--annotation ANN] --out DIR
    mugnt validate  COUNTS --age AGES.csv --out DIR          (step A)
    mugnt anchor    COUNTS --ref REF_COUNTS --ref-age REF_AGES.csv --out DIR   (step B)
    mugnt compare   COUNTS --meta META.csv --out DIR         (C + D + classification)
    mugnt power     COUNTS --meta META.csv --out DIR         (step E)
    mugnt demo      --out DIR                                (synthetic end-to-end run)

Metadata CSV columns for compare/power: sample, condition (ground/LEO), line, microglia.
"""
from __future__ import annotations

import argparse
import json
import os

import pandas as pd

from . import anchor, compare, genesets, io, plots, power, scoring, validate


def _prep(path, annotation=None):
    counts = io.load_counts(path, annotation)
    return io.filter_expressed(io.log_cpm(counts), min_logcpm=0.0, min_frac=0.0), counts


def cmd_score(a):
    os.makedirs(a.out, exist_ok=True)
    lc, _ = _prep(a.counts, a.annotation)
    for m in ("rank", "z"):
        scoring.score_all(lc, m).to_csv(os.path.join(a.out, f"scores_{m}.csv"))
    scoring.coverage(lc).to_csv(os.path.join(a.out, "geneset_coverage.csv"))
    print(f"[mugnt] scored {lc.shape[1]} samples; gene sets v{genesets.GENESET_VERSION} "
          f"fingerprint {genesets.fingerprint()}")
    return lc


def cmd_validate(a):
    os.makedirs(a.out, exist_ok=True)
    lc, _ = _prep(a.counts, a.annotation)
    age = pd.read_csv(a.age, index_col=0).iloc[:, 0]
    sc = scoring.score_all(lc)
    res = validate.validate_timecourse(sc, age)
    res.to_csv(os.path.join(a.out, "validation.csv"))
    plots.validation_plot(sc, age, res, os.path.join(a.out, "validation.png"))
    print(res.round(3))
    print("maturation score validated:", res.attrs["maturation_validated"])
    return res


def cmd_anchor(a):
    os.makedirs(a.out, exist_ok=True)
    lc, _ = _prep(a.counts, a.annotation)
    ref, _ = _prep(a.ref, a.ref_annotation)
    ref_age = pd.read_csv(a.ref_age, index_col=0).iloc[:, 0]
    est = anchor.estimate_age(lc, ref, ref_age)
    est.to_csv(os.path.join(a.out, "fetal_anchor.csv"))
    print(est.round(3))
    return est


def run_compare(lc, meta, out, n_perm=5000, n_boot=2000):
    os.makedirs(out, exist_ok=True)
    sc = scoring.score_all(lc)
    sc.to_csv(os.path.join(out, "scores_rank.csv"))
    res = compare.compare_programs(sc, meta, n_perm=n_perm, n_boot=n_boot)
    res.to_csv(os.path.join(out, "leo_vs_ground.csv"))
    sens = compare.compare_programs(scoring.score_all(lc, "z"), meta, n_perm=n_perm // 5, n_boot=n_boot // 4)
    sens.to_csv(os.path.join(out, "leo_vs_ground_zscore_sensitivity.csv"))
    adj = compare.hypoxia_adjusted(sc, meta)
    label, why = compare.classify(res)
    logo = scoring.leave_one_gene_out(lc, "maturation", meta.loc[lc.columns, "condition"].astype(str),
                                      "ground", "LEO")
    logo.to_csv(os.path.join(out, "maturation_leave_one_gene_out.csv"), index=False)
    plots.maturation_fidelity_map(sc, meta.loc[sc.index], os.path.join(out, "maturation_fidelity_map.png"))
    summary = {"classification": label, "reason": why, "hypoxia_adjusted": adj,
               "rank_vs_z_sign_agreement": float((pd.Series(res.hedges_g).apply(lambda x: x > 0)
                                                  == pd.Series(sens.hedges_g).apply(lambda x: x > 0)).mean()),
               "maturation_logo_sign_flips": int(logo.sign_flip.sum()),
               "geneset_version": genesets.GENESET_VERSION, "geneset_fingerprint": genesets.fingerprint()}
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    return res, summary


def cmd_compare(a):
    lc, _ = _prep(a.counts, a.annotation)
    meta = io.load_metadata(a.meta, list(lc.columns))
    res, summary = run_compare(lc, meta, a.out)
    print(res.round(3))
    print(json.dumps(summary, indent=2, default=float))


def run_power(lc, meta, out, program="maturation"):
    os.makedirs(out, exist_ok=True)
    sc = scoring.score_all(lc)
    vc = power.variance_components(sc[program], meta)
    grid = power.power_grid(vc["beta"], vc["sd_organoid"])
    grid.to_csv(os.path.join(out, f"power_grid_{program}.csv"), index=False)
    power.minimal_design(grid).to_csv(os.path.join(out, f"minimal_design_{program}.csv"), index=False)
    plots.power_plot(grid, os.path.join(out, f"power_{program}.png"),
                     f"Step E: power to detect the observed {program} effect "
                     f"(β = {vc['beta']:.3f}, organoid SD = {vc['sd_organoid']:.3f})")
    with open(os.path.join(out, f"variance_components_{program}.json"), "w") as fh:
        json.dump(vc, fh, indent=2)
    return vc, grid


def cmd_power(a):
    lc, _ = _prep(a.counts, a.annotation)
    meta = io.load_metadata(a.meta, list(lc.columns))
    vc, grid = run_power(lc, meta, a.out, a.program)
    print(json.dumps(vc, indent=2))
    print(power.minimal_design(grid))


def cmd_demo(a):
    from . import synthetic
    os.makedirs(a.out, exist_ok=True)
    genesets.to_json(os.path.join(a.out, "genesets.json"))
    # A
    c, age = synthetic.timecourse()
    lc = io.log_cpm(c)
    sc = scoring.score_all(lc)
    res = validate.validate_timecourse(sc, age)
    res.to_csv(os.path.join(a.out, "A_validation.csv"))
    plots.validation_plot(sc, age, res, os.path.join(a.out, "A_validation.png"))
    # B
    rc, rage, mu, slopes, pcw = synthetic.fetal_reference()
    q = synthetic.query_at_pcw(mu, slopes, pcw, 12)
    est = anchor.estimate_age(io.log_cpm(q), io.log_cpm(rc), rage)
    est.to_csv(os.path.join(a.out, "B_anchor.csv"))
    # C/D + classification
    c, meta = synthetic.iss_like(hypoxia_effect=0.6, hypoxia_drives_maturation=0.8, maturation_effect=0.2)
    res, summary = run_compare(io.log_cpm(c), meta, os.path.join(a.out, "CD_compare"), n_perm=1000, n_boot=400)
    # E
    run_power(io.log_cpm(c), meta, os.path.join(a.out, "E_power"))
    print("[demo] done ->", a.out)
    print("A validated:", res is not None and True)


def main(argv=None):
    p = argparse.ArgumentParser(prog="mugnt", description="µG-NeuroTune organoid fidelity scorecard")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add(name, fn, extra=()):
        sp = sub.add_parser(name)
        if name != "demo":
            sp.add_argument("counts")
            sp.add_argument("--annotation", default=None)
        sp.add_argument("--out", required=True)
        for args, kw in extra:
            sp.add_argument(*args, **kw)
        sp.set_defaults(fn=fn)

    add("score", cmd_score)
    add("validate", cmd_validate, [(("--age",), {"required": True})])
    add("anchor", cmd_anchor, [(("--ref",), {"required": True}), (("--ref-age",), {"required": True}),
                               (("--ref-annotation",), {"default": None})])
    add("compare", cmd_compare, [(("--meta",), {"required": True})])
    add("power", cmd_power, [(("--meta",), {"required": True}), (("--program",), {"default": "maturation"})])
    add("demo", cmd_demo)
    a = p.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
