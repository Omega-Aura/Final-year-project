#!/usr/bin/env python
"""Collect the cand_003 clean-receptor redock: best pose per receptor, Discovery Studio files, contacts.

Runs in the WSL `ligprep` env (RDKit + Meeko are blocked natively on this machine).
Run from the project root:  bash scripts/wsl_run.sh python scripts/cand003_redock_export.py

Inputs : 04_docking/cand003_redock/{R}_seed{S}.log, {R}_seed{S}_out.pdbqt, {R}_fixed.pdb
Outputs: 10_results/cand_003_redock_clean/  (complex PDBs, ligand PDBs, 9-mode PDBs, CSVs)
"""
import csv, glob, math, os, re, subprocess, sys
import numpy as np
from rdkit import Chem

SRC = "04_docking/cand003_redock"
OUT = "10_results/cand_003_redock_clean"
NAMES = {"2V5Z": "MAO-B", "2Z5X": "MAO-A", "7JXX": "TTBK1", "7Q8Y": "TTBK2"}
NATIVE = {"2V5Z": "SAG", "2Z5X": "HRM", "7JXX": "VP7", "7Q8Y": "9IV"}
GROUP = {"2V5Z": "mao", "2Z5X": "mao", "7JXX": "ttbk", "7Q8Y": "ttbk"}
os.makedirs(OUT, exist_ok=True)


def modes(log):
    rows = []
    for l in open(log):
        m = re.match(r"\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$", l)
        if m:
            rows.append((int(m[1]), float(m[2])))
    return rows


def ligand_mols(pdbqt, tag):
    """PDBQT -> RDKit mols with correct bond orders, via Meeko's own exporter (obabel would emit
    wildcard atoms for Meeko glue atoms)."""
    sdf = f"{OUT}/_{tag}.sdf"
    subprocess.run(["mk_export.py", pdbqt, "-s", sdf], check=True, capture_output=True)
    mols = [m for m in Chem.SDMolSupplier(sdf, removeHs=False) if m is not None]
    os.remove(sdf)
    return mols


def pdb_lines(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def write_complex(rec_pdb, mol, out, resname="C03"):
    """Receptor (ATOM records) + TER + ligand HETATM + CONECT. Serials continue past the receptor's."""
    rec = pdb_lines(rec_pdb)
    last = int(rec[-1][6:11])
    blk = Chem.MolToPDBBlock(mol).splitlines()   # default flavor: CONECT for every bond (flavor 2 drops single bonds)
    het, con = [], []
    for l in blk:
        if l.startswith(("ATOM", "HETATM")):
            n = int(l[6:11]) + last
            l = "HETATM" + f"{n:5d}" + l[11:17] + f"{resname:>3s}" + " L" + f"{1:4d}" + l[26:]
            het.append(l)
        elif l.startswith("CONECT"):
            nums = [int(l[i:i + 5]) + last for i in range(6, len(l.rstrip()), 5)]
            con.append("CONECT" + "".join(f"{n:5d}" for n in nums))
    with open(out, "w") as fh:
        fh.writelines(rec)
        fh.write("TER\n")
        fh.write("\n".join(het) + "\n" + "\n".join(con) + "\nEND\n")
    return het


def contacts(rec_lines, mol):
    """Residues within 4.0 A of the ligand, and heavy-atom N/O ... N/O pairs <= 3.3 A as H-bond candidates."""
    conf = mol.GetConformer()
    lig = [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    res, hb = {}, []
    for l in rec_lines:
        el = l[76:78].strip() or l[12:16].strip()[0]
        if el == "H":
            continue
        p = xyz(l)
        key = f"{l[17:20]}{int(l[22:26])}"
        for i, (s, q) in enumerate(lig):
            d = float(np.linalg.norm(p - q))
            if d <= 4.0:
                res[key] = min(res.get(key, 9), d)
            if d <= 3.3 and s in "NO" and el in ("N", "O"):
                hb.append((f"{key}:{l[12:16].strip()}", f"lig {s}{i + 1}", round(d, 2)))
    return res, hb


rows, summary = [], []
for R in NAMES:
    seeds = {}
    for log in sorted(glob.glob(f"{SRC}/{R}_seed*.log")):
        s = int(re.search(r"seed(\d+)", log)[1])
        seeds[s] = modes(log)
    best_seed = min(seeds, key=lambda s: seeds[s][0][1])
    scores = [seeds[s][0][1] for s in sorted(seeds)]
    mols = ligand_mols(f"{SRC}/{R}_seed{best_seed}_out.pdbqt", f"{R}")
    rec_pdb = f"{SRC}/{R}_fixed.pdb"
    write_complex(rec_pdb, mols[0], f"{OUT}/{R}_{NAMES[R]}_cand_003_best_pose_complex.pdb")
    # ligand alone, every mode as a MODEL (Discovery Studio can step through the 9 poses)
    with open(f"{OUT}/{R}_{NAMES[R]}_cand_003_all9_poses.pdb", "w") as fh:
        for i, m in enumerate(mols, 1):
            fh.write(f"MODEL     {i:4d}\n")
            fh.write("".join(l + "\n" for l in Chem.MolToPDBBlock(m).splitlines() if l.startswith(("HETATM", "ATOM", "CONECT"))))
            fh.write("ENDMDL\n")
    Chem.MolToPDBFile(mols[0], f"{OUT}/{R}_{NAMES[R]}_cand_003_best_pose_ligand.pdb")
    Chem.MolToMolFile(mols[0], f"{OUT}/{R}_{NAMES[R]}_cand_003_best_pose_ligand.mol")

    # where is it relative to the crystallographic ligand?
    nat = [xyz(l) for l in open(f"03_receptors/{GROUP[R]}/{R}/native_{NATIVE[R]}.pdb") if l.startswith("HETATM")]
    lig_c = mols[0].GetConformer().GetPositions().mean(0)
    off = float(np.linalg.norm(np.mean(nat, 0) - lig_c))
    res, hb = contacts(pdb_lines(rec_pdb), mols[0])
    summary.append(dict(receptor=R, target=NAMES[R], best_seed=best_seed, seeds=seeds, off_native=off,
                        contacts=res, hbonds=hb, mode1=seeds[best_seed][0][1]))
    for s in sorted(seeds):
        for mode, aff in seeds[s]:
            rows.append(dict(receptor=R, target=NAMES[R], seed=s, mode=mode, affinity_kcal_mol=aff))

with open(f"{OUT}/all_modes_scores.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

import json
json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1, default=str)
print(f"{'receptor':8} {'target':6} {'best':>7} {'seeds (11/22/33)':>26} {'centroid off native':>20}  H-bond candidates")
for s in summary:
    sc = "/".join(f"{s['seeds'][k][0][1]:.2f}" for k in sorted(s["seeds"]))
    print(f"{s['receptor']:8} {s['target']:6} {s['mode1']:7.2f} {sc:>26} {s['off_native']:17.2f} A  "
          + "; ".join(f"{a}-{b} {d}" for a, b, d in s["hbonds"]))
    print("         contact residues:", " ".join(sorted(s["contacts"], key=lambda k: int(re.sub(r'\D', '', k)))))

# ---- MAO-A supplement -------------------------------------------------------------------------------
# FAD was deleted from the receptor as instructed, so nothing stopped the top MAO-A modes from sitting
# inside the cofactor's crystal position (mao_a_fad_clash_scan.py: modes 1-5 of every seed overlap it).
# The best mode with zero ligand atoms within 2.5 A of the crystal FAD is seed 11, mode 6 (-8.77).
mols = ligand_mols(f"{SRC}/2Z5X_seed11_out.pdbqt", "2Z5X_m6")
write_complex(f"{SRC}/2Z5X_fixed.pdb", mols[5], f"{OUT}/2Z5X_MAO-A_cand_003_FADclashfree_pose_complex.pdb")
Chem.MolToPDBFile(mols[5], f"{OUT}/2Z5X_MAO-A_cand_003_FADclashfree_pose_ligand.pdb")
print("MAO-A FAD-clash-free pose written (seed 11 mode 6)")
