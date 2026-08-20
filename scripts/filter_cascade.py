#!/usr/bin/env python
"""
Reconstruction of the Phase 4 filtering cascade (§A4). The original
filtering code was never recovered -- only its output CSVs survived. This
was rebuilt from published/standard definitions and verified to reproduce
all four pass/fail columns with zero mismatches against the recovered
00_library/reinvent4_output/campaign2_v2/rl_v2_filtered_full.csv (405
molecules), including the exact funnel counts (369 Lipinski / 81 BBB / 334
GI / 253 alert-free / 56 passing all). See LOGBOOK.md for the verification.

Filters:
  - Lipinski Ro5: pass if <=1 of {MW>500, WLogP>5, HBD>5, HBA>10}
  - GI absorption / BBB penetration: point-in-polygon against the digitized
    BOILED-Egg ellipses (Daina & Zoete, ChemMedChem 2016), see
    filtering/boiled_egg_coords.py. TPSA computed with includeSandP=True.
  - Structural alerts: RDKit FilterCatalog, PAINS + BRENK combined.

Usage: filter_cascade.py --csv 01_smiles/candidates_56.csv -o OUT.csv
       (--csv file needs an "id" and "smiles" column)
"""
import argparse
import os
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from shapely.geometry import Point, Polygon

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "filtering"))
from boiled_egg_coords import BBB_COORDS, GIA_COORDS

GIA_POLY = Polygon(GIA_COORDS)
BBB_POLY = Polygon(BBB_COORDS)

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
CATALOG = FilterCatalog(params)


def score_one(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    tpsa = Descriptors.TPSA(m, includeSandP=True)
    wlogp = Crippen.MolLogP(m)
    mw = Descriptors.MolWt(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)

    violations = sum([mw > 500, wlogp > 5, hbd > 5, hba > 10])
    lipinski_pass = violations <= 1

    gia_pass = GIA_POLY.contains(Point(tpsa, wlogp))
    bbb_pass = BBB_POLY.contains(Point(tpsa, wlogp))

    matches = CATALOG.GetMatches(m)
    alert_names = ";".join(sorted({match.GetDescription() for match in matches}))
    n_alerts = len(matches)

    return {
        "MW": mw, "WLogP": wlogp, "TPSA": tpsa, "HBD": hbd, "HBA": hba,
        "Lipinski_violations": violations, "Lipinski_pass": lipinski_pass,
        "GIA_pass": gia_pass, "BBB_pass": bbb_pass,
        "n_structural_alerts": n_alerts, "structural_alerts": alert_names,
        "alert_free": n_alerts == 0,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True)
    p.add_argument("-o", required=True)
    a = p.parse_args()

    df = pd.read_csv(a.csv)
    rows = []
    for _, r in df.iterrows():
        result = score_one(r["smiles"])
        if result is None:
            print(f"[FAIL parse] {r.get('id', r['smiles'])}", file=sys.stderr)
            continue
        result["id"] = r.get("id", r["smiles"])
        result["smiles"] = r["smiles"]
        rows.append(result)

    out = pd.DataFrame(rows)
    out.to_csv(a.o, index=False)

    n = len(out)
    print(f"total: {n}")
    print(f"Lipinski_pass: {out['Lipinski_pass'].sum()}")
    print(f"BBB_pass: {out['BBB_pass'].sum()}")
    print(f"GIA_pass: {out['GIA_pass'].sum()}")
    print(f"alert_free: {out['alert_free'].sum()}")
    all_pass = out['Lipinski_pass'] & out['BBB_pass'] & out['GIA_pass'] & out['alert_free']
    print(f"passing all: {all_pass.sum()}")


if __name__ == "__main__":
    main()
