"""Calibration: how often is the evidence classification right on CRS-19-sized data?

Each simulated mission has the CRS-19 DA design: 2 lines x microglia(+/-) x ground/LEO
x 3 organoids = 24 organoids. Four known truths:
  null        no LEO effect                              -> 'insufficient evidence'
  true_mat    real maturation, nothing else              -> 'promising maturation'
  hyp_driven  maturation produced by a hypoxia response  -> 'mixed state'
  id_loss     maturation + loss of DA identity           -> 'dysregulated state'
Two FDR cut-offs are compared. Selection rule fixed in advance: keep the cut-off with
the best WORST-CASE accuracy across the four truths (chosen on simulation only, before
any real data are scored).
"""
import json
import sys

import numpy as np

from mug_neurotune import compare, io, scoring, synthetic

SCEN = {
    "null": dict(maturation_effect=0.0),
    "true_mat": dict(maturation_effect=0.8),
    "hyp_driven": dict(maturation_effect=0.1, hypoxia_effect=0.8, hypoxia_drives_maturation=0.9),
    "id_loss": dict(maturation_effect=0.8, identity_effect=-0.8),
}
EXPECT = {"null": "insufficient evidence", "true_mat": "promising maturation",
          "hyp_driven": "mixed state", "id_loss": "dysregulated state"}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/calibration.json"
CUTS = (0.10, 0.05)

labels = {q: {s: [] for s in SCEN} for q in CUTS}
for name, kw in SCEN.items():
    for seed in range(100, 100 + N):
        c, meta = synthetic.iss_like(seed=seed, **kw)
        res = compare.compare_programs(scoring.score_all(io.log_cpm(c)), meta, n_perm=2000, n_boot=1000)
        for q in CUTS:
            compare.Q_MAX = q
            labels[q][name].append(compare.classify(res)[0])
compare.Q_MAX = 0.10

out = {"n_missions_per_truth": N, "design": "2 lines x MG+/- x ground/LEO x 3 organoids", "cutoffs": {}}
for q in CUTS:
    block = {}
    for name in SCEN:
        L = labels[q][name]
        acc = sum(x == EXPECT[name] for x in L) / N
        lo, hi = acc - 1.96 * np.sqrt(acc * (1 - acc) / N), acc + 1.96 * np.sqrt(acc * (1 - acc) / N)
        block[name] = {"expected": EXPECT[name], "accuracy": round(acc, 3),
                       "ci95": [round(max(lo, 0), 3), round(min(hi, 1), 3)],
                       "labels": {k: L.count(k) for k in sorted(set(L))}}
    block["worst_case_accuracy"] = min(b["accuracy"] for b in block.values())
    out["cutoffs"][f"q<{q}"] = block
best = max(out["cutoffs"], key=lambda k: out["cutoffs"][k]["worst_case_accuracy"])
out["selected_cutoff"] = best
json.dump(out, open(OUT, "w"), indent=2)
print(json.dumps(out, indent=2))
