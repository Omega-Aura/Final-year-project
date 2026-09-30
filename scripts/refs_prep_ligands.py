#!/usr/bin/env python
"""Prepare each reference inhibitor (and cand_003) SEPARATELY for the reference redock.

One directory per ligand, nothing shared or reused between ligands and nothing taken from 02_ligands/:
    refs/ligands/<id>/<id>_3d.sdf     RDKit ETKDG conformers -> lowest MMFF94 energy, explicit H
    refs/ligands/<id>/<id>_ph74.sdf   Open Babel protonation at pH 7.4
    refs/ligands/<id>/<id>.pdbqt      Meeko: AutoDock types, polar-H merge, Gasteiger charges, torsion tree
    refs/ligands/<id>/report.json     validation record
Every ligand is validated before it may be docked (exit status non-zero on any failure):
  input SMILES parses; one fragment; stereocentres of the 3D result match the input SMILES; heavy-atom count
  survives; PDBQT is one ligand with a ROOT and a TORSDOF; PDBQT partial charges sum to the formal charge of the
  pH 7.4 protonated molecule; every heteroatom-bound hydrogen is kept as HD; no atom type is missing.
Run in WSL (env ligprep) from the project root:
    python scripts/refs_prep_ligands.py
"""
import collections, json, os, subprocess, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "04_docking", "cand003_redock", "refs", "ligands")
REFS = pd.read_csv(os.path.join(ROOT, "01_smiles", "references.csv")).set_index("id")
CAND = pd.read_csv(os.path.join(ROOT, "01_smiles", "candidates_56.csv")).set_index("id")
# docking id -> (SMILES, source)
LIGS = {i: (REFS.loc[i, "smiles"], "01_smiles/references.csv") for i in
        ["safinamide", "selegiline", "rasagiline", "lazabemide", "isatin", "kaempferol", "quercetin",
         "clorgyline", "harmine", "VP7", "DTQ"]}
LIGS["9IV"] = (REFS.loc["9IV_ttbk1", "smiles"], "01_smiles/references.csv (9IV_ttbk1)")
LIGS["cand_003"] = (CAND.loc["cand_003", "smiles"], "01_smiles/candidates_56.csv")

ok_all = True


def cips(m):
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    # carbon centres only: a protonated tertiary amine N is flagged as a stereocentre by RDKit but inverts
    # freely at room temperature, so it is not a real stereochemical feature
    return sorted((i, c) for i, c in Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
                  if m.GetAtomWithIdx(i).GetSymbol() == "C")


for lid, (smi, src) in LIGS.items():
    d = f"{OUT}/{lid}"
    os.makedirs(d, exist_ok=True)
    fails = []
    m0 = Chem.MolFromSmiles(smi)
    if m0 is None:
        print(f"[FAIL] {lid}: SMILES does not parse"); ok_all = False; continue
    if len(Chem.GetMolFrags(m0)) != 1:
        fails.append("input has more than one fragment")
    m = Chem.AddHs(m0)
    ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=20, params=ps))
    if not cids:
        print(f"[FAIL] {lid}: embedding failed"); ok_all = False; continue
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = min(range(len(res)), key=lambda i: res[i][1])
    w = Chem.SDWriter(f"{d}/{lid}_3d.sdf"); w.write(m, confId=cids[best]); w.close()
    subprocess.run(["obabel", f"{d}/{lid}_3d.sdf", "-O", f"{d}/{lid}_ph74.sdf", "-p", "7.4"],
                   check=True, capture_output=True)
    subprocess.run(["mk_prepare_ligand.py", "-i", f"{d}/{lid}_ph74.sdf", "-o", f"{d}/{lid}.pdbqt"],
                   check=True, capture_output=True)

    p = Chem.MolFromMolFile(f"{d}/{lid}_ph74.sdf", removeHs=False)
    if p is None:
        fails.append("protonated SDF does not parse"); p = m
    formal = Chem.GetFormalCharge(p)
    # stereochemistry: 3D perceived centres of the final molecule must equal the input's
    Chem.AssignStereochemistryFrom3D(p)
    if [c for _, c in cips(m0)] != [c for _, c in cips(Chem.RemoveHs(p))]:
        fails.append(f"stereo changed: input {cips(m0)} vs prepared {cips(Chem.RemoveHs(p))}")
    heavy_in = m0.GetNumHeavyAtoms()
    if Chem.RemoveHs(p).GetNumHeavyAtoms() != heavy_in:
        fails.append("heavy-atom count changed")

    txt = open(f"{d}/{lid}.pdbqt").read().splitlines()
    atoms = [l for l in txt if l.startswith(("ATOM", "HETATM"))]
    types = collections.Counter(l[77:79].strip() for l in atoms)
    q = sum(float(l[70:76]) for l in atoms)
    # Meeko breaks ring closures of macrocycles with a pair of pseudo-atoms (G, CG) that carry no chemistry
    n_pseudo = types["G"] + types["CG"]
    n_heavy = sum(1 for l in atoms if l[77:79].strip() not in ("H", "HD", "G", "CG"))
    if n_heavy != heavy_in: fails.append(f"PDBQT heavy atoms {n_heavy} != {heavy_in}")
    if sum(1 for l in txt if l.startswith("ROOT")) != 1 or not any(l.startswith("TORSDOF") for l in txt):
        fails.append("PDBQT torsion tree missing")
    if any(l.startswith("MODEL") for l in txt):
        fails.append("PDBQT holds more than one model")
    if abs(q - formal) > 0.1:
        fails.append(f"PDBQT charge sum {q:.3f} != formal charge {formal}")
    polar_h = sum(1 for a in p.GetAtoms() if a.GetAtomicNum() == 1 and
                  a.GetNeighbors()[0].GetAtomicNum() in (7, 8, 16))
    if types["HD"] != polar_h:
        fails.append(f"polar H in PDBQT {types['HD']} != polar H in molecule {polar_h}")
    prot = [f"{a.GetSymbol()}{a.GetIdx()}" for a in p.GetAtoms() if a.GetFormalCharge() != 0]
    rep = dict(id=lid, source=src, smiles_in=smi, formal_charge_pH74=formal, charged_atoms=prot,
               heavy_atoms=heavy_in, pdbqt_atoms=len(atoms), pdbqt_types=dict(types), pseudo_atoms=n_pseudo,
               pdbqt_charge_sum=round(q, 3), mmff_energy=round(res[best][1], 2), failures=fails)
    json.dump(rep, open(f"{d}/report.json", "w"), indent=1)
    print(f"[{'ok' if not fails else 'FAIL'}] {lid:11s} net charge {formal:+d}  heavy {heavy_in}  "
          f"pdbqt {len(atoms)} atoms  q_sum {q:+.2f}  HD {types['HD']}  {fails if fails else ''}")
    ok_all &= not fails

print("LIGANDS VALIDATED" if ok_all else "LIGAND VALIDATION FAILED")
sys.exit(0 if ok_all else 1)
