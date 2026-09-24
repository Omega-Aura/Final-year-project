#!/usr/bin/env python3
"""
REINVENT4 ExternalProcess-compatible docking-in-the-loop scoring script.

Reads SMILES from stdin (one per line). For each molecule: embed a 3D
conformer (RDKit ETKDGv3 + MMFF), protonate at pH 7.4 (obabel), convert to
PDBQT (meeko mk_prepare_ligand), dock with AutoDock Vina at the given
receptor/box, and report the best pose affinity.

Mirrors this week's already-validated local docking pattern
(scripts/prep_ligands.py, scripts/dock.sh) so RL-time scores are on the same
footing as offline consensus docking. Vina reports affinity as negative
kcal/mol (more negative = better); dock_score is reported as the positive
magnitude (-best_affinity) so higher = better, matching this component's
sigmoid(low=4, high=12) reward transform in the staged-learning TOML.

Failed molecules (bad SMILES, embed failure, docking crash) score 0.0
(worst case, never crashes the RL batch).

Prints {"version":1,"payload":{"dock_score":[...]}} to stdout, matching
scripts/bbb_score.py's I/O contract.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

# Meeko's CLI (mk_prepare_ligand) lives in this interpreter's own conda env
# (Scripts/, per this week's LOGBOOK finding) but is invoked here as a bare
# subprocess with no shell activation, so PATH won't include it unless we
# add it ourselves -- same fix scripts/dock.sh applies via $CONDA_PREFIX.
_ENV_ROOT = os.path.dirname(os.path.abspath(sys.executable))
os.environ["PATH"] = os.pathsep.join([
    os.path.join(_ENV_ROOT, "Library", "bin"),
    os.path.join(_ENV_ROOT, "Scripts"),
    _ENV_ROOT,
    os.environ.get("PATH", ""),
])

from rdkit import Chem
from rdkit.Chem import AllChem

VINA_CMD = shutil.which("vina") or shutil.which("vina.exe") or "vina"
OBABEL_CMD = shutil.which("obabel") or shutil.which("obabel.exe") or "obabel"
MK_CMD = (shutil.which("mk_prepare_ligand") or shutil.which("mk_prepare_ligand.exe")
          or "mk_prepare_ligand")

AFFINITY_RE = re.compile(r"^\s+\d+\s+(-?\d+\.\d+)", re.M)


def dock_one(args):
    smi, receptor, center, box_size, exhaustiveness, num_modes = args
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return 0.0
    try:
        mol = Chem.AddHs(mol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xC0FFEE
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=3, params=ps)
        if not len(ids):
            return 0.0
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
        best = min(range(len(res)), key=lambda i: res[i][1])

        with tempfile.TemporaryDirectory() as td:
            sdf = os.path.join(td, "lig.sdf")
            pdbqt = os.path.join(td, "lig.pdbqt")
            out_pdbqt = os.path.join(td, "out.pdbqt")

            Chem.SDWriter(sdf).write(mol, confId=ids[best])
            subprocess.run([OBABEL_CMD, sdf, "-O", sdf, "-p", "7.4"],
                            check=True, capture_output=True, timeout=60)
            subprocess.run([MK_CMD, "-i", sdf, "-o", pdbqt],
                            check=True, capture_output=True, timeout=60)

            proc = subprocess.run(
                [VINA_CMD, "--receptor", receptor, "--ligand", pdbqt,
                 "--center_x", str(center[0]), "--center_y", str(center[1]),
                 "--center_z", str(center[2]),
                 "--size_x", str(box_size[0]), "--size_y", str(box_size[1]),
                 "--size_z", str(box_size[2]),
                 "--exhaustiveness", str(exhaustiveness),
                 "--num_modes", str(num_modes),
                 "--out", out_pdbqt],
                capture_output=True, text=True, timeout=300)

            scores = [float(x) for x in AFFINITY_RE.findall(proc.stdout)]
            if not scores:
                return 0.0
            return -min(scores)
    except Exception:
        return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--center", nargs=3, type=float, required=True)
    ap.add_argument("--box_size", nargs=3, type=float, required=True)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--num_modes", type=int, default=5)
    ap.add_argument("--n_workers", type=int, default=8)
    a = ap.parse_args()

    smilies = [s.strip() for s in sys.stdin.readlines() if s.strip()]
    jobs = [(smi, a.receptor, a.center, a.box_size, a.exhaustiveness, a.num_modes)
            for smi in smilies]

    if not jobs:
        scores = []
    elif a.n_workers > 1:
        with Pool(a.n_workers) as pool:
            scores = pool.map(dock_one, jobs)
    else:
        scores = [dock_one(j) for j in jobs]

    print(json.dumps({"version": 1, "payload": {"dock_score": scores}}))


if __name__ == "__main__":
    main()
