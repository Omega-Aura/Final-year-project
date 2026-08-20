"""Safe PDBQT pose -> RDKit Mol conversion via coordinate nearest-neighbor matching
(avoids RDKit's AssignBondOrdersFromTemplate substructure search, which can hang
combinatorially on aromatic/symmetric templates when used via prolif.pdbqt_supplier)."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

TYPE_TO_ELEM = {"A":"C","C":"C","N":"N","NA":"N","OA":"O","O":"O","F":"F","HD":"H","Cl":"Cl",
                "CL":"Cl","S":"S","SA":"S","P":"P","Br":"Br","I":"I","B":"B"}

def parse_pdbqt_models(path):
    """Return list of models, each a list of (elem, x,y,z, is_polarH) in file order."""
    models = []
    cur = None
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith(("ATOM","HETATM")):
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            t = line[76:79].strip()
            elem = TYPE_TO_ELEM.get(t, t)
            cur.append((elem, x, y, z))
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
    if cur:
        models.append(cur)
    return models

def template_heavy_polarH_atoms(template_mol_h):
    """From an RDKit mol with explicit Hs (3D embedded), return list of
    (atom_idx, elem) for heavy atoms + polar (O/N-attached) hydrogens, and their coords."""
    conf = template_mol_h.GetConformer()
    atoms = []
    for a in template_mol_h.GetAtoms():
        sym = a.GetSymbol()
        if sym != "H":
            atoms.append(a.GetIdx())
        else:
            nbr = a.GetNeighbors()[0]
            if nbr.GetSymbol() in ("N","O","S"):
                atoms.append(a.GetIdx())
    coords = np.array([list(conf.GetAtomPosition(i)) for i in atoms])
    elems = [template_mol_h.GetAtomWithIdx(i).GetSymbol() for i in atoms]
    return atoms, elems, coords

def match_pose_to_template(pose_atoms, template_atom_idx, template_elems, template_coords):
    """Greedy nearest-neighbor matching by element between pose atoms (from PDBQT,
    order = file order) and template atom subset (heavy+polarH). Returns list of
    template_atom_idx values in POSE FILE ORDER (i.e. perm[i] = template atom idx
    corresponding to pose_atoms[i])."""
    pose_coords = np.array([[a[1],a[2],a[3]] for a in pose_atoms])
    pose_elems = [a[0] for a in pose_atoms]
    n = len(pose_atoms)
    assert n == len(template_atom_idx), f"atom count mismatch: pose={n} template={len(template_atom_idx)}"
    used = set()
    result = [None]*n
    for i in range(n):
        best_j, best_d = None, 1e9
        for j in range(len(template_atom_idx)):
            if j in used or template_elems[j] != pose_elems[i]:
                continue
            d = np.linalg.norm(pose_coords[i] - template_coords[j])
            if d < best_d:
                best_d = d; best_j = j
        if best_j is None:
            raise ValueError(f"no match for pose atom {i} elem {pose_elems[i]}")
        result[i] = template_atom_idx[best_j]
        used.add(best_j)
    return result

def build_posed_mol(template_mol_h, pose_atoms):
    """Build an RDKit mol with the template's bond orders/formal charges but the
    docked pose's 3D coordinates (heavy atoms + polar H only; template's nonpolar
    H's get positions from the ORIGINAL template conformer, translated by the
    matched heavy-atom centroid shift -- approximate but fine for ProLIF's distance/
    angle-based interaction detection which is dominated by heavy-atom geometry)."""
    t_atom_idx, t_elems, t_coords = template_heavy_polarH_atoms(template_mol_h)
    perm = match_pose_to_template(pose_atoms, t_atom_idx, t_elems, t_coords)  # perm[i] -> template atom idx

    mol = Chem.Mol(template_mol_h)
    new_conf = Chem.Conformer(mol.GetNumAtoms())
    # start by copying template conformer (keeps correct geometry for atoms not in the pose,
    # e.g. nonpolar H's whose exact docked position isn't in the PDBQT)
    old_conf = template_mol_h.GetConformer()
    for i in range(mol.GetNumAtoms()):
        new_conf.SetAtomPosition(i, old_conf.GetAtomPosition(i))
    # overwrite with actual docked positions for heavy+polarH atoms
    pose_coords = [(a[1],a[2],a[3]) for a in pose_atoms]
    matched_template_atoms = set()
    for i, t_idx in enumerate(perm):
        new_conf.SetAtomPosition(t_idx, pose_coords[i])
        matched_template_atoms.add(t_idx)
    mol.RemoveAllConformers()
    mol.AddConformer(new_conf, assignId=True)
    return mol, matched_template_atoms
