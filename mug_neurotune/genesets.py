"""Frozen gene-program definitions for µG-NeuroTune (version 0.1.0).

RULES (pre-registration):
  * These sets are frozen BEFORE the ISS dopaminergic data are scored.
  * Any change after scoring must bump GENESET_VERSION and be logged in
    docs/CHANGELOG_GENESETS.md together with the reason.
  * Sets are kept non-overlapping with each other on purpose, so that no gene
    can push two programs at once (checked by tests/test_genesets.py).
  * Microglia markers are QC only (they confirm +/- microglia status) and are
    never counted as "dysregulation".

Every symbol is an HGNC approved symbol. Before freezing v1.0, re-check all
symbols against HGNC / MSigDB (see docs/PREREGISTRATION.md, step 1).
"""
from __future__ import annotations

import hashlib
import json

GENESET_VERSION = "0.1.0"

GENESETS: dict[str, dict] = {
    # ------------------------------------------------------------------ axis 1
    "maturation": {
        "role": "axis",
        "description": "[v0.2: NOT VALIDATED as a maturity clock - failed Step A on Fiorenzano "
                       "2021 (rho = -0.16 whole organoid, 0.41 neurons only); kept as a descriptive "
                       "synaptic-marker score. Use 'age_clock' (mug_neurotune/clock.py).] "
                       "Pan-neuronal maturation: synaptic machinery, cytoskeleton of "
                       "post-mitotic neurons, mature ion channels. Deliberately excludes "
                       "dopaminergic genes so maturation and DA identity stay independent.",
        "source": "Synaptic/neuronal maturation genes used to stage brain organoids "
                  "(e.g. Gordon et al. 2021 Nat Neurosci; Bhaduri et al. 2020 Nature).",
        "genes": ["RBFOX3", "SYP", "SYN1", "SNAP25", "SYT1", "SV2A", "DLG4", "STMN2",
                  "MAPT", "NEFL", "NEFM", "GRIA2", "SCN2A", "KCNQ2", "CACNA1A",
                  "ATP1A3", "VAMP2", "NRXN1"],
    },
    # ------------------------------------------------------------------ axis 2
    "da_identity": {
        "role": "axis",
        "description": "Midbrain dopaminergic identity: floor-plate lineage TFs and "
                       "dopamine-synthesis/handling genes.",
        "source": "Canonical mDA markers (Fiorenzano et al. 2021 Nat Commun; "
                  "Birtele et al. 2022 Development; Marotta et al. 2024 reported TH, DDC).",
        "genes": ["TH", "DDC", "SLC6A3", "SLC18A2", "NR4A2", "LMX1A", "LMX1B",
                  "FOXA2", "EN1", "PITX3", "KCNJ6", "ALDH1A1"],
    },
    "off_target": {
        "role": "fidelity_penalty",
        "description": "Non-midbrain regional identity (forebrain, hindbrain). A rise "
                       "means loss of lineage fidelity. Fidelity = da_identity - off_target.",
        "source": "Regional patterning markers commonly used for organoid "
                  "mis-specification checks (FOXG1/EMX/DLX forebrain; HOX/PHOX2B hindbrain).",
        "genes": ["FOXG1", "EMX1", "EMX2", "DLX2", "DLX5", "HOXA2", "HOXB2", "PHOX2B"],
    },
    # ------------------------------------------------------------- state layers
    "progenitor_cycle": {
        "role": "state",
        "description": "Proliferation / progenitor state (expected to fall with maturation).",
        "source": "Cell-cycle and neural-progenitor markers; Marotta et al. 2024 report "
                  "lower VIM, SOX2, PAX6, CCND1 in LEO.",
        "genes": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA", "MCM2", "AURKB",
                  "SOX2", "PAX6", "VIM", "NES", "HES1", "HES5", "CCND1"],
    },
    "stress_glycolysis": {
        "role": "dysregulation",
        "description": "Glycolytic stress program characteristic of organoid stressed cells.",
        "source": "Bhaduri et al. 2020 Nature (PGK1 as stress marker); "
                  "Vertesy et al. 2022 EMBO J (Gruffi, glycolysis GO term).",
        "genes": ["PGK1", "ENO1", "LDHA", "ALDOA", "TPI1", "PFKP", "PKM", "HK2", "GPI"],
    },
    "stress_er": {
        "role": "dysregulation",
        "description": "ER stress / unfolded protein response.",
        "source": "Bhaduri et al. 2020 (ARCN1, GORASP2); Vertesy et al. 2022 (ER stress GO).",
        "genes": ["ARCN1", "GORASP2", "DDIT3", "ATF4", "HSPA5", "XBP1", "HERPUD1",
                  "SEL1L", "DNAJB9", "PDIA4", "ATF6"],
    },
    "hypoxia": {
        "role": "dysregulation",
        "description": "HIF-target hypoxia response. Glycolytic enzymes are EXCLUDED on "
                       "purpose (they live in stress_glycolysis) so hypoxia is not counted twice.",
        "source": "Canonical HIF targets (cf. MSigDB HALLMARK_HYPOXIA; Buffa et al. 2010 "
                  "Br J Cancer hypoxia metagene).",
        "genes": ["VEGFA", "SLC2A1", "BNIP3", "BNIP3L", "PDK1", "CA9", "ADM", "NDRG1",
                  "ANKRD37", "EGLN3", "P4HA1", "P4HA2", "ERO1A", "DDIT4", "BHLHE40",
                  "ANGPTL4", "PFKFB3", "PFKFB4", "KDM3A", "MXI1"],
    },
    "p53_dna_damage": {
        "role": "dysregulation",
        "description": "p53 / DNA-damage response (radiation is a LEO co-exposure).",
        "source": "Canonical p53 targets; Marotta et al. 2024 report lower CDKN1A, GADD45A in LEO.",
        "genes": ["CDKN1A", "GADD45A", "MDM2", "BAX", "BBC3", "TP53I3", "SESN1",
                  "SESN2", "FAS", "RRM2B"],
    },
    "glial_reactive": {
        "role": "dysregulation",
        "description": "Astroglial / reactive glia shift.",
        "source": "Astrocyte and reactive-astrocyte markers.",
        "genes": ["GFAP", "AQP4", "S100B", "SLC1A3", "SOX9", "SERPINA3", "CHI3L1",
                  "LCN2", "C3"],
    },
    # ------------------------------------------------------------------- QC only
    "microglia_qc": {
        "role": "qc",
        "description": "Confirms microglia presence; NOT a dysregulation program.",
        "source": "Microglia identity markers.",
        "genes": ["AIF1", "P2RY12", "CX3CR1", "TMEM119", "C1QA", "PTPRC", "CSF1R"],
    },
}

AXES = ("maturation", "da_identity")
DYSREGULATION = tuple(k for k, v in GENESETS.items() if v["role"] == "dysregulation")


def gene_lists() -> dict[str, list[str]]:
    return {k: list(v["genes"]) for k, v in GENESETS.items()}


def fingerprint() -> str:
    """SHA-256 of the gene sets; printed in every report so results are traceable."""
    blob = json.dumps(gene_lists(), sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def to_json(path: str) -> None:
    with open(path, "w") as fh:
        json.dump({"version": GENESET_VERSION, "fingerprint": fingerprint(),
                   "genesets": GENESETS}, fh, indent=2)
