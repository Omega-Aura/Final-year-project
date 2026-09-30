#!/usr/bin/env python
"""Collect the strict-receptor redock of the reference inhibitors (04_docking/cand003_redock/refs/docking/), cand_003 included, into one table, treated identically.

Strict preparation deletes FAD, so on MAO-A/MAO-B a pose can score well by sitting inside the cofactor's
crystal position. The same test that was applied to cand_003 (10_results/cand_003_redock_clean/
mao_a_fad_clash_scan.py: any heavy atom < 2.5 A from a crystal FAD atom, chain A of raw.pdb) is applied to
every mode of every seed here. Output: 08_analysis/refs_strict_redock.csv with, per receptor and ligand,
  best_raw          top score of the single run (what Vina reports)
  best_fad_free     best score over all 9 modes with no FAD overlap (MAO receptors only)
  top_mode_fad_atoms  heavy atoms of the seed-best top mode inside the FAD position (MAO receptors only)
"""
import csv, glob, io, os, re, statistics as st
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(ROOT, "04_docking", "cand003_redock")
FAD_RECEPTORS = {"2V5Z": "mao/2V5Z", "2Z5X": "mao/2Z5X"}


def fad_atoms(rec):
    p = os.path.join(ROOT, "03_receptors", FAD_RECEPTORS[rec], "raw.pdb")
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in io.open(p, encoding="utf-8")
                     if l.startswith("HETATM") and l[17:20] == "FAD" and l[21] == "A" and l[16] in " A"])


def modes(path):
    """[(mode, affinity, heavy-atom xyz array)]"""
    out, xs, aff, mode = [], [], None, None
    for l in io.open(path, encoding="utf-8", errors="replace"):
        if l.startswith("MODEL"):
            if xs: out.append((mode, aff, np.array(xs)))
            mode, xs = int(l.split()[1]), []
        elif l.startswith("REMARK VINA RESULT"):
            aff = float(l.split()[3])
        elif l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"):
            xs.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if xs: out.append((mode, aff, np.array(xs)))
    return out


files = {}
for p in glob.glob(os.path.join(D, "refs", "docking", "*_seed*_out.pdbqt")):
    r, lig, s = re.match(r"(\w{4})_(.+)_seed(\d+)_out\.pdbqt", os.path.basename(p)).groups()
    files.setdefault((r, lig), []).append((int(s), p))

# cand_003 comes from its own redock (04_docking/cand003_redock/, same receptor recipe); use the seed-11 run only
for p in glob.glob(os.path.join(D, "*_seed11_out.pdbqt")):
    files.setdefault((os.path.basename(p)[:4], "cand_003"), []).append((11, p))

fads = {r: fad_atoms(r) for r in FAD_RECEPTORS}
rows = []
for (r, lig), lst in sorted(files.items()):
    tops, free, top_fad, top_seed = [], [], [], []
    for s, p in sorted(lst):
        ms = modes(p)
        tops.append(ms[0][1]); top_seed.append(s)
        if r in fads:
            for mode, aff, xyz in ms:
                n = int((np.linalg.norm(xyz[:, None] - fads[r][None], axis=2).min(1) < 2.5).sum())
                if mode == 1: top_fad.append(n)
                if n == 0: free.append((aff, s, mode))
    row = dict(receptor=r, ligand=lig, n_seeds=len(tops), best_raw=min(tops),
               mean_top=round(st.mean(tops), 3), sd_top=round(st.pstdev(tops), 3))
    if r in fads:
        b = min(free) if free else None
        row.update(best_fad_free=b[0] if b else "", fad_free_seed=b[1] if b else "",
                   fad_free_mode=b[2] if b else "", top_mode_fad_atoms=max(top_fad))
    else:
        row.update(best_fad_free=min(tops), fad_free_seed=top_seed[tops.index(min(tops))], fad_free_mode=1,
                   top_mode_fad_atoms="")
    rows.append(row)
    print(r, lig, row["n_seeds"], row["best_raw"], row["best_fad_free"], row["top_mode_fad_atoms"])

with open(os.path.join(ROOT, "08_analysis", "refs_strict_redock.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
