#!/usr/bin/env python
"""Export the reference redock (single Vina run per receptor/ligand, exhaustiveness 32) as PDB files for
Discovery Studio. Run in WSL (env ligprep) from the project root:
    python scripts/refs_export_pdb.py

Output, in 04_docking/cand003_redock/refs/pdb_for_DS/:
    <REC>_<lig>_complex_mode1.pdb        prepared receptor (protein + H) + best pose as HETATM, with CONECT bonds
    <REC>_<lig>_ligand_all_modes.pdb     the ligand only, every Vina mode as a MODEL
    <REC>_receptor.pdb                   the prepared receptor alone
    (MAO only, when mode 1 overlaps the deleted FAD)
    <REC>_<lig>_complex_FADfree_mode<N>.pdb   best mode that does not sit in the cofactor's crystal position
    scores.csv                           every mode of every run
"""
import csv, glob, io, os, re, shutil, subprocess
import numpy as np
from rdkit import Chem

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
B = f"{ROOT}/04_docking/cand003_redock/refs"
OUT = f"{B}/pdb_for_DS"
os.makedirs(OUT, exist_ok=True)
FAM = {"2V5Z": "mao", "2Z5X": "mao", "7JXX": "ttbk", "7Q8Y": "ttbk"}


def fad(rec):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in io.open(f"{ROOT}/03_receptors/mao/{rec}/raw.pdb", encoding="utf-8")
                     if l.startswith("HETATM") and l[17:20] == "FAD" and l[21] == "A" and l[16] in " A"])


def scores(log):
    return [(int(m[1]), float(m[2])) for l in open(log)
            if (m := re.match(r"\s+(\d+)\s+(-?\d+\.\d+)\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s*$", l))]


def complex_pdb(rec_pdb, mol, out, resname):
    rec = [l for l in open(rec_pdb) if l.startswith(("ATOM", "HETATM"))]
    last = int(rec[-1][6:11])
    het, con = [], []
    for l in Chem.MolToPDBBlock(mol).splitlines():
        if l.startswith(("ATOM", "HETATM")):
            het.append("HETATM" + f"{int(l[6:11]) + last:5d}" + l[11:17] + f"{resname:>3s}" + " L" + f"{1:4d}" + l[26:])
        elif l.startswith("CONECT"):
            nums = [int(l[i:i + 5]) + last for i in range(6, len(l.rstrip()), 5)]
            con.append("CONECT" + "".join(f"{n:5d}" for n in nums))
    with open(out, "w") as fh:
        fh.writelines(rec); fh.write("TER\n")
        fh.write("\n".join(het) + "\n" + "\n".join(con) + "\nEND\n")


rows = []
for pq in sorted(glob.glob(f"{B}/docking/*_seed11_out.pdbqt")):
    rec, lig = re.match(r"(\w{4})_(.+)_seed11_out\.pdbqt", os.path.basename(pq)).groups()
    rec_pdb = f"{B}/receptors/{rec}/{rec}_fixed.pdb"
    shutil.copy(rec_pdb, f"{OUT}/{rec}_receptor.pdb")
    sdf = f"{OUT}/_tmp.sdf"
    subprocess.run(["mk_export.py", pq, "-s", sdf], check=True, capture_output=True)
    mols = [m for m in Chem.SDMolSupplier(sdf, removeHs=False) if m is not None]
    os.remove(sdf)
    sc = scores(f"{B}/docking/{rec}_{lig}_seed11.log")
    resname = lig[:3].upper()
    complex_pdb(rec_pdb, mols[0], f"{OUT}/{rec}_{lig}_complex_mode1.pdb", resname)
    with open(f"{OUT}/{rec}_{lig}_ligand_all_modes.pdb", "w") as fh:
        for i, m in enumerate(mols, 1):
            fh.write(f"MODEL     {i:4d}\n")
            fh.write("".join(l + "\n" for l in Chem.MolToPDBBlock(m).splitlines()
                             if l.startswith(("HETATM", "ATOM", "CONECT"))))
            fh.write("ENDMDL\n")
    note = ""
    if rec in ("2V5Z", "2Z5X"):
        f = fad(rec)
        ov = [int((np.linalg.norm(m.GetConformer().GetPositions()[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]][:, None]
                                  - f[None], axis=2).min(1) < 2.5).sum()) for m in mols]
        free = [i for i, n in enumerate(ov) if n == 0]
        note = f"mode1 overlaps FAD with {ov[0]} atoms" if ov[0] else "mode1 clear of FAD"
        if ov[0] and free:
            k = free[0]
            complex_pdb(rec_pdb, mols[k], f"{OUT}/{rec}_{lig}_complex_FADfree_mode{k + 1}.pdb", resname)
            note += f"; FAD-free pose written: mode {k + 1} ({sc[k][1]:.2f})"
        elif ov[0]:
            note += "; NO FAD-free mode among the 9"
    for mode, aff in sc:
        rows.append(dict(receptor=rec, ligand=lig, mode=mode, affinity_kcal_mol=aff))
    print(f"{rec} {lig:11s} mode1 {sc[0][1]:7.2f}  {note}")

with open(f"{OUT}/scores.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
