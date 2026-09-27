#!/usr/bin/env python3
"""
REINVENT4 ExternalProcess-compatible BBB-permeability scoring script.
Reads SMILES from stdin (one per line). For each molecule, computes RDKit
TPSA/WLogP and a SIGNED distance to the digitized Daina & Zoete (2016)
BOILED-Egg BBB ellipse boundary (positive = inside/BBB-permeant region,
negative = outside), using the same ellipse coordinates as the Phase 4
filtering pipeline (filtering/boiled_egg_coords.py) so the RL reward and
the final filter gate are numerically consistent.

A signed distance (rather than a boolean pass/fail) gives RL a gradient to
climb even for molecules currently outside the ellipse, instead of a flat
zero reward for all failing molecules regardless of how close they are.

Prints {"version":1,"payload":{"bbb_signed_dist":[...]}} to stdout.
Invalid SMILES score -50.0 (strongly negative, never crashes the batch).
"""
import sys, os, json

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "filtering"))
from boiled_egg_coords import BBB_COORDS
from shapely.geometry import Point, Polygon
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

BBB_POLY = Polygon(BBB_COORDS)

def score_one(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return -50.0
        tpsa = Descriptors.TPSA(mol)
        wlogp = Crippen.MolLogP(mol)
        pt = Point(tpsa, wlogp)
        dist = pt.distance(BBB_POLY.boundary)
        return float(dist if BBB_POLY.contains(pt) else -dist)
    except Exception:
        return -50.0

def main():
    smilies = [s.strip() for s in sys.stdin.readlines() if s.strip()]
    scores = [score_one(s) for s in smilies]
    print(json.dumps({"version": 1, "payload": {"bbb_signed_dist": scores}}))

if __name__ == "__main__":
    main()
