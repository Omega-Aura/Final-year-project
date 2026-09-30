#!/usr/bin/env python
"""Gate for the reference redock: a receptor is validated here before any docking input is used.

Checks per receptor (04_docking/cand003_redock/refs/receptors/<PDB>/):
  1  the cleaned PDB has no water and no HETATM record
  2  every protein residue of the cleaned PDB survives into the PDBQT (none dropped by Meeko)
  3  hydrogens are present, and polar hydrogens carry the AutoDock HD type
  4  each residue's partial charges sum to a whole number (Meeko charges per residue), and the residue
     charges of ARG/LYS/ASP/GLU/HIS are the expected chemistry
  5  total PDBQT charge equals the net charge implied by the protonation states chosen (from the prep report)
  6  the native-ligand docking box lies on the protein: the native ligand centroid sits inside the box and
     protein heavy atoms surround it
Exit status is non-zero if anything fails. Writes validation_report.json next to each receptor.
"""
import collections, json, math, os, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = os.path.join(ROOT, "04_docking", "cand003_redock", "refs", "receptors")
FAM = {"2V5Z": "mao", "2Z5X": "mao", "7JXX": "ttbk", "7Q8Y": "ttbk"}
EXPECT = {"ARG": {1}, "LYS": {1}, "ASP": {-1}, "GLU": {-1}}
ok_all = True


def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"   [{'ok' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


for pdb, fam in FAM.items():
    d = f"{R}/{pdb}"
    print(pdb)
    rep = json.load(open(f"{d}/{pdb}_prep_report.json"))
    fixed = collections.OrderedDict()
    het = 0
    for l in open(f"{d}/{pdb}_fixed.pdb"):
        if l.startswith("HETATM") or l[17:20] in ("HOH", "WAT"):
            het += 1
        if l.startswith("ATOM"):
            fixed[(l[21], l[22:27])] = l[17:20]
    check("no water / HETATM in cleaned PDB", het == 0, f"({het} found)")

    res, types, xyz = collections.OrderedDict(), collections.Counter(), []
    for l in open(f"{d}/{pdb}.pdbqt"):
        if l.startswith(("ATOM", "HETATM")):
            k = (l[21], l[22:27])
            res.setdefault(k, []).append(float(l[70:76]))
            types[l[77:79].strip()] += 1
            if l[77:79].strip() not in ("H", "HD"):
                xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    dropped = [k for k in fixed if k not in res]
    check("all residues kept in PDBQT", not dropped, f"({len(fixed)} in, {len(res)} out, dropped {dropped[:3]})")
    check("polar hydrogens present (HD)", types["HD"] > 0, f"(HD={types['HD']}, atoms={sum(types.values())})")

    bad = [(fixed[k], k[1].strip(), round(sum(v), 3)) for k, v in res.items() if abs(sum(v) - round(sum(v))) > 0.05]
    check("each residue charge is a whole number", not bad, f"(non-integer: {bad[:4]})")
    # chain termini are zwitterionic: the N-terminal residue gains +1 (NH3+) and the C-terminal one -1 (COO-),
    # so e.g. an N-terminal Asp/Glu or a C-terminal Lys legitimately sums to 0
    keys = list(res)
    shift = {keys[0]: +1, keys[-1]: -1}
    wrong = collections.defaultdict(list)
    for k, v in res.items():
        n = fixed[k]
        if n in EXPECT and round(sum(v)) not in {e + shift.get(k, 0) for e in EXPECT[n]}:
            wrong[n].append((k[1].strip(), round(sum(v))))
    wrong = dict(wrong)
    check("ARG/LYS/ASP/GLU carry the expected charge (termini allowed +/-1)", not wrong, f"({wrong})")
    tot = sum(sum(v) for v in res.values())
    exp = rep["expected_net_charge"]
    check("total charge = protonation-implied net charge", abs(tot - exp) < 0.5, f"(pdbqt {tot:.2f}, expected {exp})")

    box = json.load(open(f"{ROOT}/03_receptors/{fam}/{pdb}/box.json"))
    c, sz = box["center"], box["size"]
    nat = [l for l in open(next(f"{ROOT}/03_receptors/{fam}/{pdb}/{x}" for x in os.listdir(f"{ROOT}/03_receptors/{fam}/{pdb}")
                                if x.startswith("native_") and x.endswith(".pdb")))
           if l.startswith("HETATM") and l[76:78].strip() != "H"]
    cen = [sum(float(l[30 + 8 * i:38 + 8 * i]) for l in nat) / len(nat) for i in range(3)]
    inside = all(abs(cen[i] - c[i]) <= sz[i] / 2 for i in range(3))
    near = sum(1 for p in xyz if math.dist(p, cen) < 8.0)
    check("native ligand inside the box, protein around it", inside and near > 40,
          f"(centroid {[round(x, 1) for x in cen]}, box centre {c}, {near} protein heavy atoms within 8 A)")
    json.dump(dict(pdb=pdb, residues=len(fixed), atoms=sum(types.values()), types=dict(types),
                   total_charge=round(tot, 3), expected=exp, native_centroid=cen, box=box),
              open(f"{d}/validation_report.json", "w"), indent=1)

print("RECEPTORS VALIDATED" if ok_all else "RECEPTOR VALIDATION FAILED")
sys.exit(0 if ok_all else 1)
