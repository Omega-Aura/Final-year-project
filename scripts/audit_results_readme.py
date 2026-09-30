# -*- coding: utf-8 -*-
"""Check 10_results/README.md against the primary files and against the manuscript.

A failing check is a hypothesis about a defect, not a defect -- this session has already produced
three verifier bugs that each looked like a real finding.
"""
import csv
import io
import sys

import pandas as pd
from scipy import stats

R = io.open("10_results/README.md", encoding="utf-8").read()
MS = io.open("09_manuscript/manuscript_dualtarget_flavonol.md", encoding="utf-8").read()
md = list(csv.DictReader(io.open("07_mmgbsa/md_mmgbsa_summary.csv", encoding="utf-8")))
vd = pd.read_csv("08_analysis/vinardo_crosscheck.csv")

fails, checks = [], 0


def check(label, ok, detail=""):
    global checks
    checks += 1
    if not ok:
        fails.append("%s  %s" % (label, detail))


# ---- pose verdict counts must match the source table ----
stable = sum(1 for r in md if r["pose_verdict"] == "stable")
diss = sum(1 for r in md if r["pose_verdict"] == "dissociates")
check("17 systems", len(md) == 17, "got %d" % len(md))
check("9 stable matches 'Nine of seventeen'", stable == 9 and u"Nine of seventeen" in R,
      "stable=%d" % stable)
check("2 dissociate matches 'two leave the site'", diss == 2 and u"two leave the site" in R,
      "dissociates=%d" % diss)

# ---- every arm now has a replicate pair: limitation 4's new claim ----
names = set(r["system_dir"].split("/")[-1] for r in md)
for pair in ("system_MAOA_p3_r2", "system_MAOB_p2_r2", "system_TTBK1_r2", "system_TTBK2_p2_r2"):
    check("replicate exists: %s" % pair, pair in names, "not in summary table")
check("limitation 4 no longer claims TTBK lacks a replicate",
      "no on-pose replicate pair beyond TTBK1 pose 1" not in R)
check("limitation 4 names all four arms",
      "MAO-A" in R and "TTBK2 pose 2" in R and "Two replicates per arm at most" in R)

# ---- Vinardo numbers recomputed, not copied ----
EXP = {"7JXXdry": (0.489, 0.354, 8, "cand_013", "cand_002", 8, 12, 1.87, 0.60),
       "2V5Zdry": (0.690, 0.142, 6, "cand_043", "cand_001", 2, 10, 2.92, 0.77)}
for rec, (pa, pt, ov, bv, bd, c3v, c3d, sa, st) in EXP.items():
    m = vd[vd.receptor == rec]
    r_all = stats.pearsonr(m.score_vina, m.score_vinardo)[0]
    t = m.nsmallest(15, "score_vina")
    r_top, p_top = stats.pearsonr(t.score_vina, t.score_vinardo)
    check("%s all56 r = %.3f" % (rec, pa), abs(r_all - pa) < 0.0005, "recomputed %.4f" % r_all)
    check("%s top15 r = %.3f" % (rec, pt), abs(r_top - pt) < 0.0005, "recomputed %.4f" % r_top)
    check("%s top15 not significant" % rec, p_top > 0.05, "p=%.3f" % p_top)
    o = len(set(t.ligand) & set(m.nsmallest(15, "score_vinardo").ligand))
    check("%s overlap %d" % (rec, ov), o == ov, "recomputed %d" % o)
    check("%s best vina %s" % (rec, bv), m.loc[m.score_vina.idxmin(), "ligand"] == bv)
    check("%s best vinardo %s" % (rec, bd), m.loc[m.score_vinardo.idxmin(), "ligand"] == bd)
    c3 = m[m.ligand == "cand_003"]
    check("%s cand_003 %d vs %d" % (rec, c3v, c3d),
          int(c3.rank_vina.iloc[0]) == c3v and int(c3.rank_vinardo.iloc[0]) == c3d,
          "got %d vs %d" % (int(c3.rank_vina.iloc[0]), int(c3.rank_vinardo.iloc[0])))
    s_all = m.score_vina.max() - m.score_vina.min()
    s_top = t.score_vina.max() - t.score_vina.min()
    check("%s spread %.2f -> %.2f" % (rec, sa, st),
          abs(s_all - sa) < 0.005 and abs(s_top - st) < 0.005,
          "recomputed %.3f -> %.3f" % (s_all, s_top))

worst = int(max(vd[vd.ligand == "cand_003"].rank_vina.max(),
                vd[vd.ligand == "cand_003"].rank_vinardo.max()))
check("cand_003 worst rank 12 (README says top 12 of 56)", worst == 12, "got %d" % worst)

# ---- Result 5 present and correctly framed ----
check("Result 5 section exists", "## Result 5" in R)
check("Result 5 states non-reproducibility", "not reproducible across scoring functions" in R)
check("Result 5 warns against the all-56 reading", "wrong one to read" in R)
check("Result 5 rules out pooling", "never pooled or averaged" in R)
check("Result 5 rules out a third function", "third" in R and "adjudicated" in R)
check("claim 11 present", "| 11 |" in R)
check("claim 2 cross-references claim 11", "see claim 11" in R)
check("limitation 7 present", "7. **One docking engine, two scoring functions" in R)
check("lead-compound ranks carry the caveat", "Read those ranks with Result 5" in R)

# ---- consistency with the manuscript ----
check("five false positives, matching manuscript 4.3",
      "**five** documented instances" in R and "Five separate false positives" in MS)
check("README and manuscript agree on 9 of 17",
      "Nine of seventeen" in R and "Nine of seventeen" in MS)
for n in ["+0.489", "+0.690", "+0.354", "+0.142", "8 of 15", "6 of 15"]:
    check("manuscript and README share %s" % n, n in R and n in MS, "missing in one")

# ---- the removed Done bullets must not have taken data with them ----
for n in ["0.26", "0.24", "3.97", "0.53", "1.45"]:
    check("%s still present after bullet removal" % n, n in R, "lost from the file")
check("Done bullets gone", "**Done 2026-09-2" not in R)

# ---- the forward-looking section holds only forward-looking items ----
sec = R.split("## What would change the answer")[1].split("## Relevant files")[0]
check("no completed work left in 'What would change the answer'",
      "Done 2026" not in sec and "DONE" not in sec)

# ---- Result 6 (strict-preparation redock): recomputed from the Vina logs and crystal coordinates ----
import glob
import re

import numpy as np

CD = "04_docking/cand003_redock"


def top_score(log):
    for l in io.open(log, encoding="utf-8", errors="replace"):
        m = re.match(r"\s+1\s+(-?\d+\.\d+)\s", l)
        if m:
            return float(m.group(1))


best6 = {}
for r in ("2V5Z", "2Z5X", "7JXX", "7Q8Y"):
    sc = [top_score(f) for f in sorted(glob.glob("%s/%s_seed*.log" % (CD, r)))]
    check("Result 6: %s has 3 seed logs with a mode-1 score" % r, len(sc) == 3 and None not in sc, str(sc))
    best6[r] = min(sc)
    check("Result 6: %s seeds within 0.06 kcal/mol" % r, max(sc) - min(sc) <= 0.061, str(sc))
for r, want in (("2V5Z", -11.15), ("2Z5X", -11.44), ("7JXX", -8.49), ("7Q8Y", -9.87)):
    check("Result 6: %s best %.2f" % (r, want), round(best6[r], 2) == want, "got %.3f" % best6[r])

# earlier same-ligand dry (cofactor-kept) scores that the table quotes
for d, want in (("2V5Zdry", -11.42), ("2Z5Xdry", -8.24), ("7JXXdry", -8.41), ("7Q8Ydry", -9.89)):
    sc = [top_score("04_docking/candidates/%s_candidates_56_seed%d/cand_003.log" % (d, s)) for s in (11, 22, 33)]
    check("Result 6: earlier dry %s = %.2f" % (d, want), round(min(sc) - 1e-9, 2) == want or round(min(sc), 2) == want,
          str(sc))

# MAO-A: per-mode overlap with the crystal FAD, recomputed from the docked poses
fad = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in io.open("03_receptors/mao/2Z5X/raw.pdb", encoding="utf-8")
                if l.startswith("HETATM") and l[17:20] == "FAD" and l[21] == "A" and l[16] in " A"])
modes = []
for f in sorted(glob.glob(CD + "/2Z5X_seed*_out.pdbqt")):
    seed, xs, aff, mode = int(re.search(r"seed(\d+)", f).group(1)), [], None, None

    def flush():
        if xs:
            d = np.linalg.norm(np.array(xs)[:, None] - fad[None], axis=2).min(1)
            modes.append((seed, mode, aff, int((d < 2.5).sum()), float(d.min())))
    for l in io.open(f, encoding="utf-8", errors="replace"):
        if l.startswith("MODEL"):
            flush(); mode, xs = int(l.split()[1]), []
        elif l.startswith("REMARK VINA RESULT"):
            aff = float(l.split()[3])
        elif l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"):
            xs.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    flush()
top5 = [m for m in modes if m[1] <= 5]
check("Result 6: MAO-A modes 1-5 of every seed overlap FAD", len(top5) == 15 and all(m[3] > 0 for m in top5),
      "n=%d" % len(top5))
check("Result 6: overlap counts span 4-12 atoms", min(m[3] for m in top5) == 4 and max(m[3] for m in top5) == 12,
      str(sorted(set(m[3] for m in top5))))
check("Result 6: closest FAD approach 0.3-0.56 A", 0.25 < min(m[4] for m in top5) < 0.35 and max(m[4] for m in top5) < 0.57,
      "%.2f..%.2f" % (min(m[4] for m in top5), max(m[4] for m in top5)))
free = min((m for m in modes if m[3] == 0), key=lambda m: m[2])
check("Result 6: best FAD-clash-free MAO-A mode is -8.77", round(free[2], 2) == -8.77, str(free))
seeds_free = {m[0] for m in modes if m[3] == 0 and round(m[2], 2) == -8.77}
check("Result 6: -8.77 clash-free in all 3 seeds", seeds_free == {11, 22, 33}, str(seeds_free))
check("Result 6: MAO-B leads clean MAO-A by 2.38", round(best6["2V5Z"] - free[2], 2) == -2.38,
      "%.3f" % (best6["2V5Z"] - free[2]))
check("Result 6: raw ranking would reverse claim 5", best6["2Z5X"] < best6["2V5Z"])
check("Result 6: TTBK2 leads TTBK1 by 1.38", round(best6["7Q8Y"] - best6["7JXX"], 2) == -1.38,
      "%.3f" % (best6["7Q8Y"] - best6["7JXX"]))

# README wording for Result 6
check("Result 6 section exists", "## Result 6" in R)
check("claim 12 present", "| 12 |" in R)
check("limitation 8 present", "8. **Result 6 is one ligand" in R)
check("Result 6 says the MAO-A top score is an artefact", "headline score is an artefact" in R)
check("Result 6 states the −8.77 / 2.38 / 3.18 figures",
      all(s in R for s in (u"−8.77", "2.38", "3.18")))
check("'Nothing new was computed' is qualified", "except Result 6" in R)

print("checks run: %d" % checks)
if fails:
    print("\nFAILED (%d):" % len(fails))
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("ALL PASS -- 10_results/README.md reconciles with the primary files and the manuscript")
