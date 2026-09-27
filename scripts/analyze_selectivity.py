#!/usr/bin/env python
"""Selectivity margins for manuscript Table 3.7, on the minimal-water receptors.

Margin convention follows the existing manuscript: margin = on_target - anti_target, so a
POSITIVE margin means the anti-target binds more strongly, i.e. UNFAVOURABLE selectivity.

Receptors (minimal-water protocol adopted 2026-09-25):
    TTBK1   7JXXdry   0 waters      MAO-B   2V5Zdry   0 waters
    TTBK2   7Q8Ybrg   5 bridging    MAO-A   2Z5Xdry   0 waters
    TTBK1   7Q8V      1 water (required; second TTBK1 point, B4 cross-check pair with 7Q8Y)
"""
import glob, os, re, sys
import numpy as np
import pandas as pd

SHORTLIST_N = 15

# A Vina run that cannot place the ligand still reports a "best" mode -- often near zero or
# frankly positive (cand_013 on 2Z5Xdry: -0.449, then +0.71, +1.36, +12.99, +45.17, and only
# 5 modes instead of 9). That is a steric failure, not a weak affinity, and averaging it into
# a margin manufactures a selectivity number out of a non-measurement. Such ligands are
# reported separately as "no viable pose" -- which is itself a real result, just a censored one.
FAIL_THRESHOLD = -5.0


def consensus(receptor, ligandset="candidates_56"):
    rec = f"03_receptors/{receptor}/receptor.pdbqt"
    if not os.path.exists(rec):
        return None
    rmt = os.path.getmtime(rec)
    rows, stale = [], 0
    for log in glob.glob(f"04_docking/{receptor}_{ligandset}_seed*/*.log"):
        if os.path.getmtime(log) <= rmt:
            stale += 1
            continue
        s = [float(x) for x in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", open(log).read(), re.M)]
        if s:
            rows.append({"ligand": os.path.basename(log)[:-4], "best": min(s)})
    if not rows:
        return None
    df = (pd.DataFrame(rows).groupby("ligand")["best"]
          .agg(score="mean", sd="std", n="count").reset_index())
    if stale:
        print(f"  [note] {receptor}: ignored {stale} stale log(s)", file=sys.stderr)
    return df


def margin_block(label, on_name, on_rec, anti_name, anti_rec, shortlist):
    on, anti = consensus(on_rec), consensus(anti_rec)
    if on is None or anti is None:
        have = f"{on_rec}={'ok' if on is not None else 'MISSING'}, " \
               f"{anti_rec}={'ok' if anti is not None else 'MISSING'}"
        print(f"\n{label}: not available ({have})")
        return None
    m = on.merge(anti, on="ligand", suffixes=("_on", "_anti"))
    m = m[m.ligand.str.startswith("cand_")].copy()
    m["margin"] = m.score_on - m.score_anti
    incomplete = m[(m.n_on < 3) | (m.n_anti < 3)]
    m["favorable"] = m.margin < 0

    fail_on = m.score_on > FAIL_THRESHOLD
    fail_anti = m.score_anti > FAIL_THRESHOLD
    failed = m[fail_on | fail_anti]
    ok = m[~(fail_on | fail_anti)]

    print(f"\n{label}   [{on_name} {on_rec} vs {anti_name} {anti_rec}]   n={len(m)}"
          + (f"   ({len(incomplete)} ligand(s) with <3 seeds)" if len(incomplete) else ""))
    if len(failed):
        print(f"  EXCLUDED {len(failed)} ligand(s) with no viable pose "
              f"(best > {FAIL_THRESHOLD} kcal/mol) -- steric failure, not weak binding:")
        for _, r in failed.iterrows():
            where = anti_name if r.score_anti > FAIL_THRESHOLD else on_name
            print(f"      {r.ligand}: {on_name} {r.score_on:+.2f}  "
                  f"{anti_name} {r.score_anti:+.2f}   (fails in {where})")
        print(f"      -> these are sterically excluded from {where}, which is a qualitative "
              f"selectivity result but not a delta-G")
    sub = ok[ok.ligand.isin(shortlist)] if shortlist is not None else ok
    print(f"  quantifiable n={len(ok)}   mean margin {ok.margin.mean():+.2f}   "
          f"favourable {ok.favorable.sum()}/{len(ok)}")
    if shortlist is not None and len(sub):
        print(f"  shortlist   n={len(sub)}   mean margin {sub.margin.mean():+.2f}   "
              f"favourable {sub.favorable.sum()}/{len(sub)}")
    m["no_pose"] = fail_on | fail_anti
    return m[["ligand", "score_on", "score_anti", "margin", "favorable"]] \
        .rename(columns={"score_on": f"{on_name}", "score_anti": f"{anti_name}",
                         "margin": f"margin_{label.replace(' ', '_')}"})


print("=" * 78)
print("SELECTIVITY MARGINS -- minimal-water receptors")
print("positive margin = anti-target binds more strongly = UNFAVOURABLE")
print("=" * 78)

ttbk1 = consensus("7JXXdry")
maob = consensus("2V5Zdry")
if ttbk1 is None or maob is None:
    sys.exit("on-target runs missing (7JXXdry / 2V5Zdry) -- nothing to compare against")

# Shortlist: top N by the dual-target sum, the manuscript's selection basis.
dual = ttbk1[["ligand", "score"]].merge(maob[["ligand", "score"]], on="ligand",
                                        suffixes=("_ttbk1", "_maob"))
dual = dual[dual.ligand.str.startswith("cand_")].copy()
dual["dual"] = dual.score_ttbk1 + dual.score_maob
shortlist = set(dual.nsmallest(SHORTLIST_N, "dual").ligand)
print(f"\nshortlist = top {SHORTLIST_N} by TTBK1+MAO-B combined score")
print("  " + ", ".join(sorted(shortlist)))

blocks = []
b = margin_block("TTBK1 vs TTBK2", "TTBK1", "7JXXdry", "TTBK2", "7Q8Ybrg", shortlist)
blocks.append(b)
b = margin_block("MAO-B vs MAO-A", "MAO-B", "2V5Zdry", "MAO-A", "2Z5Xdry", shortlist)
blocks.append(b)
b = margin_block("TTBK1(7Q8V) vs TTBK2", "TTBK1_7Q8V", "7Q8V", "TTBK2", "7Q8Ybrg", shortlist)
blocks.append(b)

# B4 cross-check: the +0.784 offset was measured on 9IV alone, native to both TTBK structures.
# If the candidate margins cluster near it, the offset is a property of the receptor pair; if
# they scatter, it is a property of 9IV and should not be applied as a blanket correction.
xc = blocks[2]
if xc is not None:
    col = [c for c in xc.columns if c.startswith("margin_")][0]
    v = xc[col]
    print(f"\nB4 cross-check (7Q8V - 7Q8Ybrg across 56 candidates vs +0.784 measured on 9IV):")
    print(f"  mean {v.mean():+.3f}   sd {v.std():.3f}   median {v.median():+.3f}   "
          f"range {v.min():+.3f} .. {v.max():+.3f}")
    print(f"  within +/-0.5 of +0.784: {(v.sub(0.784).abs() <= 0.5).sum()}/{len(v)}")

out = None
for b in blocks:
    if b is None:
        continue
    keep = ["ligand"] + [c for c in b.columns if c != "ligand"]
    out = b[keep] if out is None else out.merge(
        b[["ligand"] + [c for c in b.columns if c.startswith("margin_")]], on="ligand", how="outer")
if out is not None:
    out["shortlist"] = out.ligand.isin(shortlist)
    out.to_csv("08_analysis/selectivity_margins.csv", index=False)
    print("\nwrote 08_analysis/selectivity_margins.csv")
