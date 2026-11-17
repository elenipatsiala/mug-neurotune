"""Ground-truth tests: the pipeline must recover planted effects and not invent absent ones."""
import numpy as np
import pandas as pd

from mug_neurotune import anchor, compare, genesets, io, power, scoring, synthetic, validate


def test_genesets_disjoint_and_wellformed():
    sets = genesets.gene_lists()
    seen = {}
    for name, genes in sets.items():
        assert len(genes) == len(set(genes)), f"duplicate gene inside {name}"
        for g in genes:
            assert g not in seen, f"{g} in both {seen.get(g)} and {name}"
            seen[g] = name
            assert g == g.upper() and " " not in g


def test_rank_score_bounds_and_library_size_invariance():
    c, _ = synthetic.timecourse(reps=2)
    lc = io.log_cpm(c)
    s1 = scoring.score_all(lc)
    s2 = scoring.score_all(io.log_cpm(c * 7))
    assert s1.drop(columns="fidelity").abs().max().max() <= 0.5
    assert np.allclose(s1.values, s2.values, atol=1e-6)


def test_no_rank_leakage():
    """Raising maturation+hypoxia genes must not move programs that were not touched."""
    c, _ = synthetic.iss_like(seed=7, maturation_effect=0.0)
    boosted = c.copy().astype(float)
    sets = genesets.gene_lists()
    boosted.loc[sets["maturation"] + sets["hypoxia"]] *= 4      # same organoids, only two programs changed
    a = scoring.score_all(io.log_cpm(c))
    b = scoring.score_all(io.log_cpm(boosted))
    # Library size changes rescale every untouched gene AND every background gene by the same
    # factor (a monotone transform), so their relative ranks - and hence scores - are identical.
    for prog in ("stress_glycolysis", "stress_er", "p53_dna_damage", "glial_reactive", "off_target",
                 "da_identity", "progenitor_cycle"):
        assert np.allclose(a[prog], b[prog], atol=1e-12), prog
    assert (b["maturation"] > a["maturation"]).all() and (b["hypoxia"] > a["hypoxia"]).all()


def test_A_validation_passes_on_true_maturation():
    c, age = synthetic.timecourse()
    res = validate.validate_timecourse(scoring.score_all(io.log_cpm(c)), age)
    assert res.attrs["maturation_validated"]
    assert res.loc["progenitor_cycle", "rho"] < 0


def test_A_validation_fails_without_maturation():
    c, age = synthetic.timecourse(maturation_slope=0.0)
    res = validate.validate_timecourse(scoring.score_all(io.log_cpm(c)), age)
    assert not res.attrs["maturation_validated"]


def test_B_anchor_recovers_age():
    rc, rage, mu, slopes, pcw = synthetic.fetal_reference()
    ref = io.log_cpm(rc)
    for target in (8, 12, 16):
        q = synthetic.query_at_pcw(mu, slopes, pcw, target, seed=target)
        est = anchor.estimate_age(io.log_cpm(q), ref, rage)
        assert abs(est.est_pcw.mean() - target) < 2.0, (target, est.est_pcw.mean())


def test_compare_detects_planted_maturation_and_classifies():
    c, meta = synthetic.iss_like(maturation_effect=0.8)
    sc = scoring.score_all(io.log_cpm(c))
    res = compare.compare_programs(sc, meta, n_perm=500, n_boot=300)
    assert res.loc["maturation", "hedges_g"] > 0.5 and res.loc["maturation", "ci_low"] > 0
    assert res.loc["progenitor_cycle", "hedges_g"] < 0
    label, _ = compare.classify(res)
    assert label == "promising maturation"


def test_compare_null_has_no_maturation_effect():
    c, meta = synthetic.iss_like(maturation_effect=0.0, seed=11)
    sc = scoring.score_all(io.log_cpm(c))
    res = compare.compare_programs(sc, meta, n_perm=500, n_boot=300)
    label, _ = compare.classify(res)
    assert label == "insufficient evidence"


def test_D_mixed_state_and_hypoxia_adjustment():
    c, meta = synthetic.iss_like(maturation_effect=0.1, hypoxia_effect=0.8, hypoxia_drives_maturation=0.9)
    sc = scoring.score_all(io.log_cpm(c))
    res = compare.compare_programs(sc, meta, n_perm=500, n_boot=300)
    label, _ = compare.classify(res)
    assert label == "mixed state"
    adj = compare.hypoxia_adjusted(sc, meta)
    assert adj["share_explained"] > 0.5


def test_D_adjustment_leaves_independent_effect():
    c, meta = synthetic.iss_like(maturation_effect=0.8, hypoxia_effect=0.0, seed=5)
    adj = compare.hypoxia_adjusted(scoring.score_all(io.log_cpm(c)), meta)
    assert adj["share_explained"] < 0.5


def test_E_power_monotone_in_donors():
    p2 = power.simulate_power(0.05, 0.05, 0.05, n_donors=2, n_org=8, n_sim=2000)
    p8 = power.simulate_power(0.05, 0.05, 0.05, n_donors=8, n_org=8, n_sim=2000)
    assert p8 > p2


def test_E_variance_components_reasonable():
    c, meta = synthetic.iss_like(maturation_effect=0.8)
    sc = scoring.score_all(io.log_cpm(c))
    vc = power.variance_components(sc["maturation"], meta)
    assert vc["beta"] > 0 and vc["sd_organoid"] > 0
