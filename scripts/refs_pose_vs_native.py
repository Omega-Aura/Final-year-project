#!/usr/bin/env python
"""Sanity check of the reference redock: how far do the docked poses of a receptor's own crystal ligand lie from the
crystal pose? Symmetry-aware heavy-atom RMSD with NO re-alignment (poses stay in the receptor frame).
Run in WSL (env ligprep) from the project root:
    python scripts/refs_pose_vs_native.py
"""
import glob, os, re
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = f"{ROOT}/04_docking/cand003_redock/refs/docking"
# receptor, family, docked ligand id, native ligand code
CASES = [("7Q8Y", "ttbk", "9IV", "9IV"), ("7JXX", "ttbk", "VP7", "VP7"), ("2Z5X", "mao", "harmine", "HRM")]
for rec, fam, lid, code in CASES:
    ref = Chem.MolFromMolFile(f"{ROOT}/03_receptors/{fam}/{rec}/native_{code}_ref.sdf", removeHs=True)
    per_seed, allm = [], []
    for f in sorted(glob.glob(f"{D}/{rec}_{lid}_seed*_out.pdbqt")):
        seed = int(re.search(r"seed(\d+)", f).group(1))
        pm = PDBQTMolecule.from_file(f, skip_typing=True)
        m = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
        for conf in [c.GetId() for c in m.GetConformers()]:
            r = rdMolAlign.CalcRMS(m, ref, prbId=conf)
            allm.append((r, seed, conf + 1, pm._pose_data["free_energies"][conf]))
            if conf == 0:
                per_seed.append(round(r, 2))
    b = min(allm)
    print(f"{rec} {lid}: mode-1 RMSD to the crystal pose per seed {per_seed} A; closest of all modes "
          f"{b[0]:.2f} A (seed {b[1]}, mode {b[2]}, score {b[3]:.2f})")
