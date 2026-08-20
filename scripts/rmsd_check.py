#!/usr/bin/env python
"""Symmetry-corrected RMSD of docked poses vs a crystal reference.
   Naive atom-order RMSD inflates correct poses on symmetric rings — always use this."""
import argparse, glob
from spyrmsd import io, rmsd

p = argparse.ArgumentParser()
p.add_argument("--ref", required=True)
p.add_argument("--poses", required=True)
a = p.parse_args()

ref = io.loadmol(a.ref)
ref.strip()
best = (None, 1e9)
for f in sorted(glob.glob(a.poses)):
    for i, m in enumerate(io.loadallmols(f)):
        m.strip()
        r = rmsd.symmrmsd(ref.coordinates, m.coordinates,
                          ref.atomicnums, m.atomicnums,
                          ref.adjacency_matrix, m.adjacency_matrix)
        print(f"{f}  pose {i+1}  RMSD {r:.2f}")
        if r < best[1]:
            best = (f"{f}:{i+1}", r)
print(f"\nBEST {best[0]}  {best[1]:.2f} A  ->  {'PASS' if best[1] < 2.0 else 'FAIL'}")
