#!/usr/bin/env python
"""Wet-vs-dry water-shell test: is the 5 A native-ligand water shell helping or hurting?

Arm 1 (strong): the 56 candidates on 7JXX wet vs dry. Vina scores normally correlate
NEGATIVELY with ligand size -- more atoms, more contacts, better score. If the wet receptor
inverts that while the dry one does not, the waters are sterically excluding ligands rather
than shaping chemistry.

Arm 2 (weak, reported with its caveat): references against measured IC50. Only 7 compounds
carry use_in_correlation=yes and they split across targets, so n per receptor is 1-3. Stated,
not spun.
"""
import glob, os, re, sys
import numpy as np
import pandas as pd
from scipy import stats


def consensus(receptor, ligandset):
    """Mean over seeds of each ligand's best pose, using only logs newer than the receptor."""
    rec = f"03_receptors/{receptor}/receptor.pdbqt"
    if not os.path.exists(rec):
        return pd.DataFrame(columns=["ligand", "score", "n"])
    rmt = os.path.getmtime(rec)
    rows = []
    for log in glob.glob(f"04_docking/{receptor}_{ligandset}_seed*/*.log"):
        if os.path.getmtime(log) <= rmt:
            continue
        s = [float(x) for x in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", open(log).read(), re.M)]
        if s:
            rows.append({"ligand": os.path.basename(log)[:-4], "best": min(s)})
    if not rows:
        return pd.DataFrame(columns=["ligand", "score", "n"])
    return (pd.DataFrame(rows).groupby("ligand")["best"]
            .agg(score="mean", n="count").reset_index())


def r_line(label, x, y):
    if len(x) < 3:
        return f"  {label:28s} n={len(x)}  (too few points to correlate)"
    r, p = stats.pearsonr(x, y)
    return f"  {label:28s} n={len(x)}  r={r:+.3f}  p={p:.2e}"


print("=" * 72)
print("ARM 1  -- size bias on the 56 candidates, 7JXX wet vs dry")
print("=" * 72)

props = pd.read_csv("01_smiles/candidates_56.csv")[["id", "MW", "TPSA", "Rings"]] \
          .rename(columns={"id": "ligand"})
wet = consensus("7JXX", "candidates_56").rename(columns={"score": "wet"})
dry = consensus("7JXXdry", "candidates_56").rename(columns={"score": "dry"})

if dry.empty:
    print("  no dry candidate runs yet -- run scripts/run_water_test.sh first")
else:
    m = wet[["ligand", "wet"]].merge(dry[["ligand", "dry"]], on="ligand").merge(props, on="ligand")
    print(f"\n{len(m)} candidates docked against both receptors\n")
    print("  score vs ligand size (negative = Vina's normal size bias, bigger scores better)")
    for col in ("MW", "TPSA", "Rings"):
        print(r_line(f"WET  score vs {col}", m[col], m.wet))
        print(r_line(f"DRY  score vs {col}", m[col], m.dry))
        print()
    m["delta"] = m.wet - m.dry
    print("  effect of adding the waters (delta = wet - dry; + = waters made it worse)")
    print(f"    mean {m.delta.mean():+.3f}   median {m.delta.median():+.3f}   "
          f"sd {m.delta.std():.3f}   range {m.delta.min():+.3f} .. {m.delta.max():+.3f}")
    print(r_line("delta vs MW", m.MW, m.delta))
    print(r_line("delta vs TPSA", m.TPSA, m.delta))
    rho = m.wet.corr(m.dry, method="spearman")
    top = len(set(m.nsmallest(15, "wet").ligand) & set(m.nsmallest(15, "dry").ligand))
    print(f"\n  Spearman(wet, dry) = {rho:.3f}    top-15 overlap = {top}/15")
    m.sort_values("delta", ascending=False).to_csv("08_analysis/water_test_candidates.csv",
                                                   index=False)
    print("  wrote 08_analysis/water_test_candidates.csv")

print()
print("=" * 72)
print("ARM 2  -- correlation with measured IC50 (WEAK: n=1-3 per target)")
print("=" * 72)

refs = pd.read_csv("01_smiles/references.csv")
refs = refs[refs.use_in_correlation.astype(str).str.lower() == "yes"].copy()
to_nM = {"nm": 1.0, "um": 1000.0, "µm": 1000.0}
refs["nM"] = refs.apply(
    lambda r: r.measured_value * to_nM.get(str(r.measured_unit).strip().lower(), np.nan), axis=1)
refs = refs.dropna(subset=["nM"])
refs["pIC50"] = 9 - np.log10(refs.nM)

TARGET_RECEPTOR = {"MAO-B": "2V5Z", "MAO-A": "2Z5X", "TTBK1": "7JXX"}
rows = []
for target, rec in TARGET_RECEPTOR.items():
    sub = refs[refs.target == target]
    if sub.empty:
        continue
    w = consensus(rec, "references").rename(columns={"score": "wet"})
    d = consensus(rec + "dry", "references").rename(columns={"score": "dry"})
    if w.empty or d.empty:
        print(f"  {target} ({rec}): missing runs (wet={len(w)}, dry={len(d)}) -- skipped")
        continue
    j = sub[["id", "pIC50"]].rename(columns={"id": "ligand"}) \
           .merge(w[["ligand", "wet"]], on="ligand").merge(d[["ligand", "dry"]], on="ligand")
    j["target"] = target
    rows.append(j)
    print(f"\n  {target} ({rec}), n={len(j)}")
    for _, x in j.iterrows():
        print(f"    {x.ligand:16s} pIC50 {x.pIC50:5.2f}   wet {x.wet:+.3f}   dry {x.dry:+.3f}")

if rows:
    allj = pd.concat(rows)
    print("\n  pooled across targets (pooling different proteins is itself questionable):")
    # more negative score should mean more potent, i.e. NEGATIVE correlation with pIC50
    print(r_line("WET  score vs pIC50", allj.pIC50, allj.wet))
    print(r_line("DRY  score vs pIC50", allj.pIC50, allj.dry))
    allj.to_csv("08_analysis/water_test_references.csv", index=False)
    print("\n  wrote 08_analysis/water_test_references.csv")
    print("\n  CAVEAT: with n of this size neither correlation is decisive on its own.")
    print("  Arm 1 carries the weight; this arm is a consistency check.")
else:
    print("\n  no reference runs available on both wet and dry receptors yet")
