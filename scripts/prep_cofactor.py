#!/usr/bin/env python
"""Prepare a cofactor as rigid receptor atoms and merge it into receptor.pdbqt.

Why this is needed
------------------
Meeko's polymer path templates a receptor residue by residue. FAD is not in its residue
library, so it tries to build a template on the fly -- and that fails whenever the crystal
models the 8alpha-S-cysteinyl linkage short enough to be perceived as a covalent bond:

    2V5Z (MAO-B)  Cys397 SG--FAD C8M = 2.31 A  -> not perceived as a bond -> FAD survives
    2Z5X (MAO-A)  Cys406 SG--FAD C8M = 1.65 A  -> perceived as a bond     -> FAD DROPPED

Both are the same chemistry; only the refined bond length differs. Because
`mk_prepare_receptor` is run with `-a/--allow_bad_res`, the 2Z5X failure was silent, and
MAO-A was docked for weeks against an active site with no flavin in it.

Relying on that 2.31 A accident would also violate rule 1.5: MAO-B and MAO-A are compared
to each other, so they must be prepared the same way. This script therefore takes the
cofactor out of Meeko's polymer path entirely and prepares it through the SAME ligand
toolchain every other molecule in this project goes through (RCSB ideal chemistry ->
bond-order transfer onto crystal coordinates -> obabel protonation at pH 7.4 ->
mk_prepare_ligand), then strips the torsion tree and appends the atoms to receptor.pdbqt as
rigid receptor. Applied identically to both MAO receptors.

Usage: prep_cofactor.py <receptor_dir> <RESNAME> [--ph 7.4]
"""
import argparse
import os
import shutil
import subprocess
import sys
import urllib.request

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ap = argparse.ArgumentParser()
ap.add_argument("recdir")
ap.add_argument("resname")
ap.add_argument("--ph", type=float, default=7.4)
a = ap.parse_args()

res = a.resname.upper()
cache = "03_receptors/_cofactors"
os.makedirs(cache, exist_ok=True)

# 1. crystal coordinates of the cofactor, from the already chain/altloc-filtered file
src = f"{a.recdir}/clean_noH.pdb"
lines = [l for l in open(src)
         if l.startswith("HETATM") and l[17:20].strip() == res]
if not lines:
    sys.exit(f"[FAIL] no {res} atoms in {src}")
pdb_frag = f"{a.recdir}/cofactor_{res}.pdb"
open(pdb_frag, "w").writelines(lines + ["END\n"])

# 2. ideal chemistry from the RCSB Chemical Component Dictionary (cached)
ideal = f"{cache}/{res}_ideal.sdf"
if not os.path.exists(ideal):
    url = f"https://files.rcsb.org/ligands/download/{res}_ideal.sdf"
    urllib.request.urlretrieve(url, ideal)
    print(f"[fetch] {url}")
tmpl = Chem.MolFromMolFile(ideal)
if tmpl is None:
    sys.exit(f"[FAIL] could not parse {ideal}")
tmpl = Chem.RemoveHs(tmpl)

pose = Chem.MolFromPDBFile(pdb_frag, removeHs=False, sanitize=False)
if pose is None:
    sys.exit(f"[FAIL] could not read {pdb_frag}")
if pose.GetNumAtoms() != tmpl.GetNumAtoms():
    sys.exit(f"[FAIL] {res}: crystal has {pose.GetNumAtoms()} heavy atoms, "
             f"RCSB ideal has {tmpl.GetNumAtoms()}")

fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, pose)
Chem.SanitizeMol(fixed)
sdf = f"{a.recdir}/cofactor_{res}.sdf"
w = Chem.SDWriter(sdf)
w.write(fixed)
w.close()
print(f"[chem] {res}: {fixed.GetNumAtoms()} heavy atoms, "
      f"formal charge {Chem.GetFormalCharge(fixed)} (pre-protonation)")

# 3. protonate at pH and type through the project's normal ligand path
obabel = shutil.which("obabel") or shutil.which("obabel.exe") or "obabel"
mk = (shutil.which("mk_prepare_ligand") or shutil.which("mk_prepare_ligand.exe")
      or "mk_prepare_ligand")
subprocess.run([obabel, sdf, "-O", sdf, "-p", str(a.ph)], check=True, capture_output=True)
lig_pdbqt = f"{a.recdir}/cofactor_{res}.pdbqt"
subprocess.run([mk, "-i", sdf, "-o", lig_pdbqt], check=True, capture_output=True)

# 4. merge as RIGID receptor atoms: keep coordinates/types, drop the torsion tree
keep = [l for l in open(lig_pdbqt) if l.startswith(("ATOM", "HETATM"))]
if not keep:
    sys.exit(f"[FAIL] no atoms in {lig_pdbqt}")
# stamp the real residue name/number back on so the guard in prep_receptor.sh can see it
stamped = []
for l in keep:
    stamped.append(l[:17] + f"{res:>3}" + l[20:])

rec = f"{a.recdir}/receptor.pdbqt"
body = [l for l in open(rec) if not l.startswith(("ROOT", "ENDROOT", "BRANCH",
                                                  "ENDBRANCH", "TORSDOF"))]
already = sum(1 for l in body if l[17:20].strip() == res)
if already:
    sys.exit(f"[skip] {res} already present in {rec} ({already} atoms)")

with open(rec, "w") as fh:
    fh.writelines(body)
    fh.writelines(stamped)
print(f"[merge] {res}: {len(stamped)} atoms appended to {rec} as rigid receptor")
