#!/usr/bin/env python
"""Computed properties + gate outcomes for the reference inhibitors and cand_003, using the SAME
score_one() as the filtering cascade (Lipinski, BOILED-Egg BBB/GI, PAINS+BRENK) so the comparison table
is apples to apples. Output: 08_analysis/reference_vs_cand003_properties.csv"""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filter_cascade import score_one

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
refs = pd.read_csv(os.path.join(ROOT, "01_smiles", "references.csv"))
keep = ["safinamide", "selegiline", "rasagiline", "lazabemide", "kaempferol", "quercetin", "isatin",
        "clorgyline", "harmine", "VP7", "DTQ", "9IV_ttbk1"]
rows = []
for k in keep:
    r = refs[refs.id == k].iloc[0]
    d = score_one(r.smiles); d.update(id=k, smiles=r.smiles)
    rows.append(d)
cand = pd.read_csv(os.path.join(ROOT, "01_smiles", "candidates_56.csv"))
c = cand[cand.iloc[:, 0] == "cand_003"] if "id" not in cand.columns else cand[cand.id == "cand_003"]
smi = c["smiles"].iloc[0]
d = score_one(smi); d.update(id="cand_003", smiles=smi); rows.append(d)
out = pd.DataFrame(rows)[["id", "MW", "WLogP", "TPSA", "HBD", "HBA", "Lipinski_violations",
                          "BBB_pass", "GIA_pass", "n_structural_alerts", "smiles"]]
out.to_csv(os.path.join(ROOT, "08_analysis", "reference_vs_cand003_properties.csv"), index=False)
print(out.drop(columns="smiles").round(2).to_string(index=False))
