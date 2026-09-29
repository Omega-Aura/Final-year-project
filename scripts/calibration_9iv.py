#!/usr/bin/env python
"""B4 -- the 9IV matched-pair calibration.

The manuscript claims a compound distinguishes TTBK1 from TTBK2 by ~1.7 kcal/mol. Before
that number can mean anything, the protocol has to be shown capable of measuring a
paralog difference of that size without inventing one. 9IV is the control that makes this
testable: it is crystallised in BOTH proteins (7Q8V = TTBK1, 7Q8Y = TTBK2) and is close to
equipotent against them, so its true dDG is near zero. Whatever margin we measure for 9IV
is therefore the protocol's own systematic bias between these two receptors.

Both receptors were prepared by one run of prep_receptor.sh and are docked here with ONE
ligand file, so the measured margin cannot come from ligand preparation.

    margin = consensus(7Q8V) - consensus(7Q8Y)        [kcal/mol]

    margin ~ 0  -> protocol unbiased; selectivity margins stand as measured
    margin != 0 -> systematic bias; subtract it from every reported selectivity margin
"""
import glob
import json
import os
import re

import numpy as np
import pandas as pd

RT = 0.5925  # kcal/mol at 298 K

PAIR = {"7Q8V": "TTBK1", "7Q8Y": "TTBK2"}

# Primary source, both values from ONE paper and ONE assay format:
#   Nozal et al. 2022, J Med Chem 65(2):1585-1607, PMID 34978799,
#   DOI 10.1021/acs.jmedchem.1c01942 -- compound 42 / VNG2.73, PDB ligand 9IV.
#   Inhibition of recombinant human TTBK1 / TTBK2, RICDLHDDEEDEAMSITA substrate.
#   TTBK1 IC50 330 nM, TTBK2 IC50 490 nM. Curated as ChEMBL5200069, whose canonical
#   SMILES matches 01_smiles/references.csv exactly (C18H13ClN4O).
#
# CORRECTED 2026-09-27. This previously read (430.0, 330.0, 530.0) for TTBK1: the 430 was
# the midpoint of BindingDB's 330-530 nM range, which AGGREGATES TWO ASSAYS FROM DIFFERENT
# PAPERS. Taking a cross-paper midpoint for one arm and a single-paper value for the other
# is exactly the asymmetry this project keeps catching elsewhere -- the matched pair only
# means anything if both numbers come from the same measurement. The 330 nM value is the
# one Nozal reports alongside the 490 nM TTBK2 figure, so it is the arm-symmetric choice.
# No range is carried now, because a single assay per arm does not define one.
IC50_NM = {"7Q8V": (330.0, 330.0, 330.0), "7Q8Y": (490.0, 490.0, 490.0)}

SEEDS = [11, 22, 33]


def per_seed_best(receptor, seed):
    """Best (most negative) Vina score for 9IV in one receptor at one seed."""
    # Globbed rather than hardcoded to 04_docking/native_redock/ so that regrouping the
    # run directories cannot silently turn this calibration into a no-op. More than one
    # match means two groups hold the same run, which must not be averaged over silently.
    hits = glob.glob(f"04_docking/*/{receptor}_native_9IV_seed{seed}/native_9IV.log")
    if len(hits) > 1:
        raise SystemExit(f"[abort] {receptor} seed{seed}: {len(hits)} copies of "
                         f"native_9IV.log across groups: {hits}")
    if not hits:
        return None
    log = hits[0]
    scores = [float(x) for x in
              re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", open(log).read(), re.M)]
    return min(scores) if scores else None


rows = []
for rec in PAIR:
    for s in SEEDS:
        v = per_seed_best(rec, s)
        if v is None:
            raise SystemExit(f"missing docking result: {rec} seed{s}")
        rows.append({"receptor": rec, "protein": PAIR[rec], "seed": s, "best": v})

df = pd.DataFrame(rows)
summary = (df.groupby(["receptor", "protein"])["best"]
             .agg(consensus="mean", sd="std", n="count", best_overall="min")
             .reset_index())

c = {r.receptor: r.consensus for r in summary.itertuples()}
sd = {r.receptor: r.sd for r in summary.itertuples()}

margin = c["7Q8V"] - c["7Q8Y"]
# independent seeds, so propagate as quadrature of the two SEMs
sem = {k: sd[k] / np.sqrt(len(SEEDS)) for k in sd}
margin_sem = float(np.hypot(sem["7Q8V"], sem["7Q8Y"]))

# experimental ddG from the measured IC50s (TTBK1 minus TTBK2)
mid = RT * np.log(IC50_NM["7Q8V"][0] / IC50_NM["7Q8Y"][0])

# The two IC50s are single values from one assay and the paper publishes no error on either, so
# the experimental ddG has no measured uncertainty. Reporting the bias as a bare "1.018" would
# therefore put three significant figures on the number that gates every selectivity claim in
# this project, with no error bar -- which is precisely the failure this project documents
# elsewhere (see the SEM finding in 07_mmgbsa/README.md). Instead, propagate an EXPLICITLY
# ASSUMED IC50 precision so the reader can see what the threshold actually rests on.
#
# IC50_REL_SD is an assumption, not a measurement. 0.30 is a conventional within-assay figure
# for enzymatic IC50 replicates; both values here come from the same assay in the same paper, so
# within-assay precision is the relevant scale rather than the larger inter-laboratory spread.
# Change it and the band moves -- that is the point of having it named.
IC50_REL_SD = 0.30

# ddG = RT ln(A/B), so fractional errors on A and B add in quadrature through the log.
exp_sd = RT * np.hypot(IC50_REL_SD, IC50_REL_SD)
exp_lo, exp_hi = mid - exp_sd, mid + exp_sd

bias = margin - mid  # what the protocol adds on top of the real difference
# Dominated by the experimental term: the docking SEM is seed noise and ~17x smaller.
bias_sd = float(np.hypot(margin_sem, exp_sd))

print(df.to_string(index=False))
print()
print(summary.to_string(index=False))
print()
print(f"measured margin (7Q8V TTBK1 - 7Q8Y TTBK2) = {margin:+.3f} +/- {margin_sem:.3f} kcal/mol"
      f"   (seed noise only)")
print(f"experimental ddG from IC50                = {mid:+.3f} +/- {exp_sd:.3f} kcal/mol"
      f"   (assumed {IC50_REL_SD:.0%} IC50 precision, NOT measured)")
print(f"systematic protocol bias                  = {bias:+.3f} +/- {bias_sd:.3f} kcal/mol"
      f"   -> quote as ~{bias:.1f}")
print(f"  the bias uncertainty is {exp_sd / margin_sem:.0f}x the docking SEM, so it is set almost"
      f" entirely by the\n  experimental reference, not by the docking. Margins below"
      f" ~{bias + bias_sd:.1f} cannot be called selectivity.")
print()
if abs(margin) <= 2 * margin_sem or abs(bias) < 0.5:
    print("VERDICT: margin is not distinguishable from the near-zero experimental value.")
    print("         Protocol shows no meaningful TTBK1/TTBK2 bias at this magnitude;")
    print("         selectivity margins in 3.7 stand as measured.")
else:
    print("VERDICT: margin is systematically non-zero. Subtract the bias above from every")
    print(f"         reported TTBK1/TTBK2 selectivity margin (correction {-bias:+.3f} kcal/mol).")

os.makedirs("05_validation", exist_ok=True)
df.to_csv("05_validation/calibration_9IV.csv", index=False)
json.dump({"margin_kcal_mol": round(float(margin), 3),
           "margin_sem": round(margin_sem, 3),
           "margin_sem_basis": "spread across 3 velocity seeds; seed noise only",
           "experimental_ddG_kcal_mol": round(float(mid), 3),
           "experimental_ddG_sd": round(float(exp_sd), 3),
           "experimental_ddG_range": [round(float(exp_lo), 3), round(float(exp_hi), 3)],
           "experimental_ddG_basis": (
               f"TTBK1 {IC50_NM['7Q8V'][0]:.0f} nM vs TTBK2 {IC50_NM['7Q8Y'][0]:.0f} nM, both "
               "single values from one assay in Nozal 2022 (PMID 34978799). The paper publishes "
               f"no error on either, so the +/- is an ASSUMED {IC50_REL_SD:.0%} within-assay IC50 "
               "precision propagated through the log, not a measured uncertainty."),
           "systematic_bias_kcal_mol": round(float(bias), 3),
           "systematic_bias_sd": round(bias_sd, 3),
           "systematic_bias_quote_as": f"~{bias:.1f} kcal/mol",
           "interpretability_floor_kcal_mol": round(float(bias + bias_sd), 1),
           "per_receptor": {r.receptor: {"protein": r.protein,
                                         "consensus": round(r.consensus, 3),
                                         "sd": round(r.sd, 3)}
                            for r in summary.itertuples()}},
          open("05_validation/calibration_9IV_margin.json", "w"), indent=2)
print("\nwrote 05_validation/calibration_9IV.csv and calibration_9IV_margin.json")
