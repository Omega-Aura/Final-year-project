#!/usr/bin/env python
"""Rebuild each receptor's native-ligand reference SDF with correct bond orders.

Why this exists
---------------
`prep_receptor.sh` extracts the native ligand straight out of the crystal PDB, which
carries coordinates and element symbols but NO bond orders and NO hydrogens. Open Babel
then has to guess the chemistry from interatomic distances alone, and for fused
heteroaromatics it guesses wrong. Measured on this project's own structures:

    7Q8V 9IV  ->  Clc1ccccc1Oc1ccc(NC2=NC=NC3=NC=C[C@@H]32)cc1     WRONG
    correct   ->  Clc1ccccc1Oc1ccc(Nc2ncnc3[nH]ccc23)cc1

i.e. the pyrrolo[2,3-d]pyrimidine came back non-aromatic, the pyrrole N-H was dropped,
and a stereocentre was invented at a carbon that is planar in reality. A reference pose
with that connectivity cannot be graph-matched against a correctly-perceived docked pose,
so symmetry-corrected RMSD either fails outright or silently falls back to a bad mapping.

The fix is to stop guessing: take the bond orders from the compound's verified SMILES
(01_smiles/references.csv, each row already cross-checked against the RCSB Chemical
Component Dictionary) and transfer them onto the crystal coordinates via RDKit's
AssignBondOrdersFromTemplate. Coordinates are untouched -- this only corrects chemistry.

Output: 03_receptors/<family>/<PDB>/native_<LIG>_ref.sdf  (the RMSD reference to actually use)
"""
import glob
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

REFS = "01_smiles/references.csv"
rows = pd.read_csv(REFS)
rows = rows[rows["native_of_pdb"].notna()]

fail = 0
for _, r in rows.iterrows():
    pdb = str(r["native_of_pdb"]).strip()
    hits = (glob.glob(f"03_receptors/{pdb}/native_*.pdb") +
            glob.glob(f"03_receptors/*/{pdb}/native_*.pdb"))
    hits = [h for h in hits if not h.endswith("_ref.pdb")]
    if len(hits) != 1:
        print(f"[skip] {pdb}: expected 1 native_*.pdb, found {len(hits)}")
        fail += 1
        continue
    src = hits[0]
    lig = os.path.basename(src)[len("native_"):-len(".pdb")]

    # proximityBonding gives us connectivity from geometry; bond ORDERS come from template
    pose = Chem.MolFromPDBFile(src, removeHs=False, sanitize=False)
    if pose is None:
        print(f"[FAIL read] {pdb} {lig}")
        fail += 1
        continue

    tmpl = Chem.MolFromSmiles(str(r["smiles"]))
    if tmpl is None:
        print(f"[FAIL template] {pdb} {lig} <- {r['id']}")
        fail += 1
        continue

    n_pose = pose.GetNumAtoms()
    n_tmpl = tmpl.GetNumAtoms()
    if n_pose != n_tmpl:
        print(f"[FAIL atomcount] {pdb} {lig}: pose has {n_pose} heavy atoms, "
              f"template ({r['id']}) has {n_tmpl}. "
              f"Alternate conformations or a partial-occupancy copy left in the extract?")
        fail += 1
        continue

    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, pose)
    except Exception as e:  # noqa: BLE001 - want the reason printed, not a traceback
        print(f"[FAIL assign] {pdb} {lig}: {e}")
        fail += 1
        continue

    Chem.SanitizeMol(fixed)
    fixed.SetProp("_Name", f"{lig}_from_{pdb}")
    # Write beside the structure's own PDB, wherever the family grouping put it.
    out = os.path.join(os.path.dirname(hits[0]), f"native_{lig}_ref.sdf")
    w = Chem.SDWriter(out)
    w.write(fixed)
    w.close()

    got = Chem.MolToSmiles(Chem.RemoveHs(fixed))
    want = Chem.MolToSmiles(tmpl)
    flag = "ok " if got == want else "MISMATCH"
    print(f"[{flag}] {pdb} {lig:4s} n={n_pose:3d}  {got}")

sys.exit(1 if fail else 0)
