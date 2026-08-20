import warnings, sys, time
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import pandas as pd
import MDAnalysis as mda
import prolif
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbqt_utils import parse_pdbqt_models, build_posed_mol

receptor_pdb_files = {
    "2V5Z": "prep/2V5Z_h.pdb", "4NFM": "prep/4NFM_h.pdb", "2Z5X": "prep/2Z5X_h.pdb",
    "6U0K": "prep/6U0K_h.pdb", "7Q8Y": "prep/7Q8Y_h.pdb",
}
ligmeta = pd.read_csv("ligands/phase1_ligand_protonation.csv")

print("building protein molecules...", flush=True)
protein_mols = {}
for rec_name, pdb_path in receptor_pdb_files.items():
    u = mda.Universe(pdb_path)
    sel = u.select_atoms("protein")
    protein_mols[rec_name] = prolif.Molecule.from_mda(sel, inferrer=None)
    print(f"  {rec_name} OK", flush=True)

print("building embedded 3D ligand templates...", flush=True)
ligand_templates = {}
for _, row in ligmeta.iterrows():
    name = row["name"]
    mol = Chem.MolFromSmiles(row["chosen_protonated_smiles"])
    mol = Chem.AddHs(mol)
    ok = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    if ok != 0:
        AllChem.EmbedMolecule(mol, randomSeed=1, useRandomCoords=True, maxAttempts=200)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    ligand_templates[name] = mol
    print(f"  {name} OK ({mol.GetNumAtoms()} atoms)", flush=True)

baseline_df = pd.read_csv("docking/phase2_baseline_docking_final.csv")
print(f"total receptor-ligand pairs: {len(baseline_df)}", flush=True)

fp = prolif.Fingerprint()
all_rows, failures = [], []

for idx, row in baseline_df.iterrows():
    rec_name, lig_name, pdbqt_path = row["receptor"], row["ligand"], row["pose_file"]
    t0 = time.time()
    try:
        models = parse_pdbqt_models(pdbqt_path)
        posed_mol, matched = build_posed_mol(ligand_templates[lig_name], models[0])  # top pose only
        posed_mol_noH_ok = Chem.RemoveHs(posed_mol, sanitize=False)  # keep polar Hs matter for HBond; use full mol
        lig_mol = prolif.Molecule.from_rdkit(posed_mol)
        fp.run_from_iterable([lig_mol], protein_mols[rec_name])
        ifp_df = fp.to_dataframe()
        interactions = []
        if len(ifp_df) > 0:
            row_data = ifp_df.iloc[0]
            for col in ifp_df.columns:
                if row_data[col]:
                    res = col[1] if isinstance(col, tuple) and len(col) > 1 else str(col)
                    itype = col[2] if isinstance(col, tuple) and len(col) > 2 else ""
                    interactions.append(f"{res}:{itype}")
        all_rows.append({"receptor": rec_name, "ligand": lig_name,
                          "n_interactions": len(interactions), "interactions": ";".join(interactions)})
        print(f"  [{idx+1}/{len(baseline_df)}] {rec_name}/{lig_name}: {len(interactions)} interactions ({time.time()-t0:.1f}s)", flush=True)
    except Exception as e:
        failures.append({"receptor": rec_name, "ligand": lig_name, "error": str(e)[:300]})
        print(f"  [{idx+1}/{len(baseline_df)}] {rec_name}/{lig_name}: FAILED - {str(e)[:150]}", flush=True)

pd.DataFrame(all_rows).to_csv("docking/phase2_interaction_fingerprints.csv", index=False)
pd.DataFrame(failures).to_csv("docking/phase2_ifp_failures.csv", index=False)
print(f"\nDone. {len(all_rows)} succeeded, {len(failures)} failed.", flush=True)
