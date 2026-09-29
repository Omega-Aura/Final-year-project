# -*- coding: utf-8 -*-
"""Check every Vinardo number in the manuscript against 08_analysis/vinardo_crosscheck.csv.

Written as a separate pass from audit_manuscript.py, which covers MD/MM-GBSA. Recomputes the
statistics from the per-ligand table rather than trusting the summary the analysis script printed --
the point is to catch a number that was transcribed wrong or has drifted, so re-deriving it is the
whole job.

A failing check here is a HYPOTHESIS about a defect, not a defect. This session has already produced
three verifier bugs that each looked like a real finding, so anything that fails gets checked by hand
before any document is touched.
"""
import io
import sys

import pandas as pd
from scipy import stats

MS = io.open("09_manuscript/manuscript_dualtarget_flavonol.md", encoding="utf-8").read()
LB = io.open("LOGBOOK.md", encoding="utf-8").read()
df = pd.read_csv("08_analysis/vinardo_crosscheck.csv")

fails, checks = [], 0


def check(label, ok, detail=""):
    global checks
    checks += 1
    if not ok:
        fails.append("%s  %s" % (label, detail))


check("csv has 112 rows (56 x 2 receptors)", len(df) == 112, "got %d" % len(df))
check("both receptors present", set(df.receptor) == {"7JXXdry", "2V5Zdry"},
      "got %s" % sorted(set(df.receptor)))

EXPECT = {
    "7JXXdry": {"label": "TTBK1", "pear_all": 0.489, "spear_all": 0.463,
                "pear_top": 0.354, "overlap": 8, "bv": "cand_013", "bd": "cand_002",
                "c3v": 8, "c3d": 12, "spread_all": 1.87, "spread_top": 0.60,
                "vdo_min": -7.72, "vdo_max": -4.99, "sd_vina": 0.065, "sd_vdo": 0.012},
    "2V5Zdry": {"label": "MAO-B", "pear_all": 0.690, "spear_all": 0.516,
                "pear_top": 0.142, "overlap": 6, "bv": "cand_043", "bd": "cand_001",
                "c3v": 2, "c3d": 10, "spread_all": 2.92, "spread_top": 0.77,
                "vdo_min": -9.42, "vdo_max": -4.55, "sd_vina": 0.012, "sd_vdo": 0.009},
}

for rec, e in EXPECT.items():
    m = df[df.receptor == rec]
    lab = e["label"]
    check("%s n=56" % lab, len(m) == 56, "got %d" % len(m))

    pear = stats.pearsonr(m.score_vina, m.score_vinardo)[0]
    spear = stats.spearmanr(m.score_vina, m.score_vinardo)[0]
    check("%s Pearson all56 = %.3f" % (lab, e["pear_all"]),
          abs(pear - e["pear_all"]) < 0.0005, "recomputed %.4f" % pear)
    check("%s Spearman all56 = %.3f" % (lab, e["spear_all"]),
          abs(spear - e["spear_all"]) < 0.0005, "recomputed %.4f" % spear)

    t = m.nsmallest(15, "score_vina")
    pear_t = stats.pearsonr(t.score_vina, t.score_vinardo)[0]
    check("%s Pearson top15 = %.3f" % (lab, e["pear_top"]),
          abs(pear_t - e["pear_top"]) < 0.0005, "recomputed %.4f" % pear_t)
    # the claim is non-significance, so assert the p value, not just the r
    p_t = stats.pearsonr(t.score_vina, t.score_vinardo)[1]
    check("%s top15 p > 0.05 (claimed n.s.)" % lab, p_t > 0.05, "p = %.3f" % p_t)

    ov = len(set(t.ligand) & set(m.nsmallest(15, "score_vinardo").ligand))
    check("%s top-15 overlap = %d" % (lab, e["overlap"]), ov == e["overlap"], "recomputed %d" % ov)

    bv = m.loc[m.score_vina.idxmin(), "ligand"]
    bd = m.loc[m.score_vinardo.idxmin(), "ligand"]
    check("%s best by vina = %s" % (lab, e["bv"]), bv == e["bv"], "got %s" % bv)
    check("%s best by vinardo = %s" % (lab, e["bd"]), bd == e["bd"], "got %s" % bd)
    check("%s the two disagree" % lab, bv != bd, "both %s -- manuscript claims DIFFERENT" % bv)

    c3 = m[m.ligand == "cand_003"]
    check("%s cand_003 vina rank = %d" % (lab, e["c3v"]),
          int(c3.rank_vina.iloc[0]) == e["c3v"], "got %d" % int(c3.rank_vina.iloc[0]))
    check("%s cand_003 vinardo rank = %d" % (lab, e["c3d"]),
          int(c3.rank_vinardo.iloc[0]) == e["c3d"], "got %d" % int(c3.rank_vinardo.iloc[0]))

    sa = m.score_vina.max() - m.score_vina.min()
    st = t.score_vina.max() - t.score_vina.min()
    check("%s vina spread all56 = %.2f" % (lab, e["spread_all"]),
          abs(sa - e["spread_all"]) < 0.005, "recomputed %.3f" % sa)
    check("%s vina spread top15 = %.2f" % (lab, e["spread_top"]),
          abs(st - e["spread_top"]) < 0.005, "recomputed %.3f" % st)
    check("%s range restriction is real" % lab, st < sa, "top15 spread not narrower")

    check("%s vinardo min = %.2f" % (lab, e["vdo_min"]),
          abs(m.score_vinardo.min() - e["vdo_min"]) < 0.005,
          "recomputed %.3f" % m.score_vinardo.min())
    check("%s vinardo max = %.2f" % (lab, e["vdo_max"]),
          abs(m.score_vinardo.max() - e["vdo_max"]) < 0.005,
          "recomputed %.3f" % m.score_vinardo.max())
    check("%s all vinardo scores favourable" % lab, (m.score_vinardo < 0).all(),
          "%d non-negative" % (m.score_vinardo >= 0).sum())

    check("%s mean inter-seed SD vina = %.3f" % (lab, e["sd_vina"]),
          abs(m.sd_vina.mean() - e["sd_vina"]) < 0.0005, "recomputed %.4f" % m.sd_vina.mean())
    check("%s mean inter-seed SD vinardo = %.3f" % (lab, e["sd_vdo"]),
          abs(m.sd_vinardo.mean() - e["sd_vdo"]) < 0.0005, "recomputed %.4f" % m.sd_vinardo.mean())
    check("%s 3 seeds per ligand" % lab, (m.n_vina == 3).all() and (m.n_vinardo == 3).all(),
          "some ligand has fewer than 3 seeds")

    # the disagreement must exceed the noise for the manuscript's claim to hold
    check("%s seed noise << disagreement" % lab, m.sd_vina.mean() < st / 3.0,
          "SD %.3f vs top15 spread %.3f" % (m.sd_vina.mean(), st))

# --- cand_003's worst rank across both functions AND both targets ---
c3 = df[df.ligand == "cand_003"]
worst = int(max(c3.rank_vina.max(), c3.rank_vinardo.max()))
check("cand_003 worst rank = 12 (manuscript says top 12 of 56)", worst == 12, "got %d" % worst)
check("top ~21%% is right for rank 12/56", abs(worst / 56.0 * 100 - 21.4) < 0.6,
      "12/56 = %.1f%%" % (worst / 56.0 * 100))

# --- the numbers must actually be in the documents, as written ---
for s in ["+0.489", "+0.690", "+0.463", "+0.516", "+0.354", "+0.142",
          "8 of 15", "6 of 15", "`cand_013`", "`cand_043`", "1.87", "0.60", "2.92", "0.77"]:
    check("manuscript contains %s" % s, s in MS, "not found")
check("manuscript states non-significance", "n.s." in MS)
check("manuscript still warns rank is not reproducible",
      "fine-grained" in MS and "scoring function" in MS)
check("manuscript does NOT claim the gap is merely closed",
      "confirms rather than resolves" in MS or "confirms the caveat" in MS)
check("no pooled/averaged consensus claimed", "consensus of two scoring functions" not in MS
      or "masks a weak score" in MS)
check("prior-phase figures marked superseded", "superseded" in MS)
check("logbook entry N present", "N (the Vinardo cross-check" in LB)

print("checks run: %d" % checks)
if fails:
    print("\nFAILED (%d):" % len(fails))
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("ALL PASS -- every Vinardo number recomputes from 08_analysis/vinardo_crosscheck.csv")
