#!/usr/bin/env python
"""Strict receptor clean-up for the cand_003 redock: chain only, no water, no heteroatoms.

Steps, in the order the user specified them:
  1. remove water            2. remove every heteroatom (ligand, cofactor, ions, buffer, glycerol)
  3. find and rebuild missing atoms in existing residues; report missing residues
  4. add hydrogens at pH 7.4 (polar H come from this step; MGLTools later merges the non-polar ones)

Deliberately NOT the project's prep_receptor.sh: that one keeps the 5 A pocket water shell and the FAD
cofactor, which is what makes the ATP-site receptors validate. This script is the opposite on purpose
and its outputs are kept in 04_docking/cand003_redock/, never mixed into the main docking tree.

Runs in WSL env `recprep` (pdbfixer + openmm).
Usage: python scripts/prep_receptor_strict.py PDBID CHAIN RAW.pdb OUTDIR
"""
import json, sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

pdb, chain, raw, outdir = sys.argv[1:5]
fx = PDBFixer(filename=raw)

# keep only the docking chain (the other chains are symmetry copies in these crystals)
fx.removeChains(chainIds=[c.id for c in fx.topology.chains() if c.id != chain])

fx.findMissingResidues()
missing_res = {f"{k[0]}:{k[1]}": v for k, v in fx.missingResidues.items()}
# Missing residues = unmodelled loops/termini. Building them would be invented coordinates, far
# from the site in most cases; only missing ATOMS of residues that exist are rebuilt.
fx.missingResidues = {}

fx.findNonstandardResidues()
nonstd = [(r.name, r.id) for r, _ in fx.nonstandardResidues]
fx.replaceNonstandardResidues()

het = sorted({r.name for r in fx.topology.residues()
              if r.name not in {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
                                "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}})
fx.removeHeterogens(keepWater=False)   # waters AND every non-protein residue
fx.findMissingAtoms()
n_missing = sum(len(v) for v in fx.missingAtoms.values())
n_term = sum(len(v) for v in fx.missingTerminals.values())
detail = {f"{r.name}{r.id}": [a.name for a in v] for r, v in fx.missingAtoms.items()}
fx.addMissingAtoms()
fx.addMissingHydrogens(7.4)

out = f"{outdir}/{pdb}_fixed.pdb"
with open(out, "w") as fh:
    PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)

# net formal charge implied by the protonation states chosen (for the charge sanity check later)
pos = sum(1 for r in fx.topology.residues() if r.name in ("ARG","LYS") or
          (r.name in ("HIP",)))
res = [r.name for r in fx.topology.residues()]
nH = {r.name: {a.name for a in r.atoms()} for r in fx.topology.residues()}
q = 0
for r in fx.topology.residues():
    names = {a.name for a in r.atoms()}
    if r.name == "ARG": q += 1
    elif r.name == "LYS" and "HZ3" in names: q += 1
    elif r.name == "ASP" and "HD2" not in names: q -= 1
    elif r.name == "GLU" and "HE2" not in names: q -= 1
    elif r.name == "HIS" and {"HD1","HE2"} <= names: q += 1
n_atoms = sum(1 for _ in fx.topology.atoms())
rep = dict(pdb=pdb, chain=chain, hetero_and_water_removed=het, nonstandard_replaced=nonstd,
           missing_residues_not_built=missing_res, missing_heavy_atoms_rebuilt=n_missing,
           missing_terminal_atoms_rebuilt=n_term, missing_atoms_detail=detail,
           residues=len(res), atoms_with_H=n_atoms, expected_net_charge=q)
json.dump(rep, open(f"{outdir}/{pdb}_prep_report.json", "w"), indent=1)
print(json.dumps({k: v for k, v in rep.items() if k != "missing_atoms_detail"}, indent=1))
