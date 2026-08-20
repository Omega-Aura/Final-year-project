#!/usr/bin/env python
"""Build 3D ligands from SMILES (or from an extracted native ligand).
   --csv  : CSV with columns id,smiles
   --native : an .sdf extracted from a crystal structure
Outputs <out>/sdf/<id>.sdf and <out>/pdbqt/<id>.pdbqt
"""
import argparse, os, shutil, subprocess, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

p = argparse.ArgumentParser()
p.add_argument("--csv")
p.add_argument("--native")
p.add_argument("-o", required=True)
p.add_argument("--ph", type=float, default=7.4)
p.add_argument("--nconf", type=int, default=20)
a = p.parse_args()

os.makedirs(f"{a.o}/sdf", exist_ok=True)
os.makedirs(f"{a.o}/pdbqt", exist_ok=True)

if a.csv:
    df = pd.read_csv(a.csv)
    items = list(zip(df["id"].astype(str), df["smiles"]))
else:
    m = Chem.MolFromMolFile(a.native)
    items = [(os.path.basename(a.native).replace(".sdf", ""), Chem.MolToSmiles(m))]

# Identify mk_prepare_ligand executable
mk_cmd = (shutil.which("mk_prepare_ligand") or 
          shutil.which("mk_prepare_ligand.exe") or 
          shutil.which("mk_prepare_ligand.py") or 
          "mk_prepare_ligand")

obabel_cmd = shutil.which("obabel") or shutil.which("obabel.exe") or "obabel"

for lid, smi in items:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print(f"[FAIL parse] {lid} {smi}", file=sys.stderr)
        continue
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    ids = AllChem.EmbedMultipleConfs(m, numConfs=a.nconf, params=ps)
    if not len(ids):
        print(f"[FAIL embed] {lid}", file=sys.stderr)
        continue
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = min(range(len(res)), key=lambda i: res[i][1])       # lowest MMFF94 energy
    sdf = f"{a.o}/sdf/{lid}.sdf"
    Chem.SDWriter(sdf).write(m, confId=ids[best])
    # protonate at target pH, then convert
    subprocess.run([obabel_cmd, sdf, "-O", sdf, "-p", str(a.ph)],
                   check=True, capture_output=True)
    subprocess.run([mk_cmd, "-i", sdf,
                    "-o", f"{a.o}/pdbqt/{lid}.pdbqt"], check=True)
    print(f"[ok] {lid}")
