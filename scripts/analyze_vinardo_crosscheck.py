#!/usr/bin/env python
"""Vina vs Vinardo on the validated receptors: does the docking rank order survive a change of
scoring function?

The prior-phase answer was no -- Pearson r = 0.56 at TTBK1 and -0.28 at MAO-B over the top 15 --
but that was measured on 4NFM/2V5Z, the earlier receptor set, whose TTBK1 structure was apo and
carried no passing redocking validation. This script re-measures it on the receptors actually
reported (7JXXdry, 2V5Zdry), over all 56 candidates rather than 15.

Two things are deliberately NOT done here:

  * Vina and Vinardo scores are never pooled, averaged or plotted on a shared axis. They are
    different functions on different scales (vina roughly -7 to -12 on this set, vinardo -5 to
    -9.5). Only the CORRELATION and the RANK AGREEMENT between them are meaningful.
  * No "consensus of the two functions" is computed. Averaging two scoring functions is how the
    prior phase nearly picked a different lead: a weak vina score got masked rather than
    corroborated. Agreement is tested, not averaged.

Reads vina from 04_docking/<group>/<run>/ and vinardo from 04_docking/crosscheck/vinardo/<run>/.
Writes 08_analysis/vinardo_crosscheck.csv.

Usage: python scripts/analyze_vinardo_crosscheck.py
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

import project_paths

SET = "candidates_56"
RECEPTORS = [("7JXXdry", "TTBK1"), ("2V5Zdry", "MAO-B")]
OUT_CSV = "08_analysis/vinardo_crosscheck.csv"
TOP_N = 15          # the slice the prior-phase cross-check covered, kept for comparability


def consensus(pattern, receptor):
    """Mean of the per-seed best pose, per ligand, dropping logs older than the receptor.

    Same freshness rule as collect_results.py: a log older than the receptor.pdbqt it names was
    produced against a different preparation and is not comparable with a new one.
    """
    rec = project_paths.receptor_file(receptor)
    if rec is None or not os.path.exists(rec):
        return None, "no receptor.pdbqt for %s" % receptor
    rmt = os.path.getmtime(rec)
    rows, stale = [], 0
    for log in glob.glob(pattern):
        if os.path.getmtime(log) <= rmt:
            stale += 1
            continue
        lig = os.path.basename(log)[:-4]
        if not lig.startswith("cand_"):
            continue
        with open(log) as fh:
            scores = [float(x) for x in
                      re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", fh.read(), re.M)]
        if scores:
            rows.append({"ligand": lig, "best": min(scores)})
    if not rows:
        return None, "no fresh logs matched %s (%d stale)" % (pattern, stale)
    df = (pd.DataFrame(rows).groupby("ligand")["best"]
          .agg(score="mean", sd="std", n="count").reset_index())
    return df, ("%d stale log(s) ignored" % stale if stale else "")


print("=" * 78)
print("VINA vs VINARDO on the validated receptors -- %s" % SET)
print("=" * 78)

frames, summary = [], []
for rec, label in RECEPTORS:
    vina, msg_a = consensus("04_docking/*/%s_%s_seed*/*.log" % (rec, SET), rec)
    vdo, msg_b = consensus("04_docking/crosscheck/vinardo/%s_%s_seed*/*.log" % (rec, SET), rec)
    for m in (msg_a, msg_b):
        if m:
            print("  [note] %s: %s" % (rec, m), file=sys.stderr)
    if vina is None or vdo is None:
        print("\n%s (%s): NOT AVAILABLE -- vina=%s vinardo=%s"
              % (label, rec, "ok" if vina is not None else msg_a,
                 "ok" if vdo is not None else msg_b))
        continue

    m = vina.merge(vdo, on="ligand", suffixes=("_vina", "_vinardo"))
    short = m[(m.n_vina < 3) | (m.n_vinardo < 3)]
    pear = stats.pearsonr(m.score_vina, m.score_vinardo)
    spear = stats.spearmanr(m.score_vina, m.score_vinardo)

    # Rank agreement on the slice a reader would actually act on.
    top_vina = set(m.nsmallest(TOP_N, "score_vina").ligand)
    top_vdo = set(m.nsmallest(TOP_N, "score_vinardo").ligand)
    overlap = len(top_vina & top_vdo)
    best_vina = m.loc[m.score_vina.idxmin(), "ligand"]
    best_vdo = m.loc[m.score_vinardo.idxmin(), "ligand"]
    r_vina = m.score_vina.rank()
    r_vdo = m.score_vinardo.rank()
    c3_vina = int(r_vina[m.ligand == "cand_003"].iloc[0]) if (m.ligand == "cand_003").any() else None
    c3_vdo = int(r_vdo[m.ligand == "cand_003"].iloc[0]) if (m.ligand == "cand_003").any() else None

    print("\n%s (%s)   n=%d%s" % (label, rec, len(m),
                                  "   (%d ligand(s) with <3 seeds)" % len(short) if len(short) else ""))
    print("  vina    range %7.2f .. %7.2f   mean inter-seed SD %.3f"
          % (m.score_vina.min(), m.score_vina.max(), m.sd_vina.mean()))
    print("  vinardo range %7.2f .. %7.2f   mean inter-seed SD %.3f"
          % (m.score_vinardo.min(), m.score_vinardo.max(), m.sd_vinardo.mean()))
    print("  Pearson  r = %+.3f  (p = %.3g)" % (pear[0], pear[1]))
    print("  Spearman r = %+.3f  (p = %.3g)" % (spear[0], spear[1]))
    print("  top-%d overlap: %d/%d" % (TOP_N, overlap, TOP_N))
    print("  best by vina: %s      best by vinardo: %s%s"
          % (best_vina, best_vdo, "   (SAME)" if best_vina == best_vdo else "   (DIFFERENT)"))
    if c3_vina is not None:
        print("  the lead cand_003: vina rank %d/%d, vinardo rank %d/%d"
              % (c3_vina, len(m), c3_vdo, len(m)))

    m.insert(0, "receptor", rec)
    m.insert(1, "target", label)
    m["rank_vina"] = r_vina.astype(int)
    m["rank_vinardo"] = r_vdo.astype(int)
    frames.append(m)
    summary.append({"receptor": rec, "target": label, "n": len(m),
                    "pearson_r": round(float(pear[0]), 3),
                    "pearson_p": float("%.3g" % pear[1]),
                    "spearman_r": round(float(spear[0]), 3),
                    "spearman_p": float("%.3g" % spear[1]),
                    "top15_overlap": overlap,
                    "best_vina": best_vina, "best_vinardo": best_vdo,
                    "cand_003_rank_vina": c3_vina, "cand_003_rank_vinardo": c3_vdo})

if not frames:
    sys.exit("no receptor had both vina and vinardo results -- nothing to compare")

out = pd.concat(frames, ignore_index=True)
out.to_csv(OUT_CSV, index=False)
print("\nwrote %s (%d rows)" % (OUT_CSV, len(out)))

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
s = pd.DataFrame(summary)
print(s.to_string(index=False))
print("""
How to read this. A high correlation means the two functions agree on the ORDER, not that they
agree on a value -- the scales differ and are never compared directly. A low or negative
correlation is a real negative result: it says the fine-grained rank within this high-scoring
slice is a property of the scoring function, not of the chemistry, and that selecting a single
best compound on one function's ranking is unsafe. Either way, absolute scores from the two
functions must not be pooled, averaged or plotted on a shared axis.""")
