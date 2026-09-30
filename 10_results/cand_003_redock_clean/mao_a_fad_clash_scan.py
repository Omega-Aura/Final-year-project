"""MAO-A only: score and FAD overlap for every mode of every seed (FAD was deleted from the receptor, so nothing in docking stopped a pose from entering its space)."""
import numpy as np, re, glob
fad = np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in open("03_receptors/mao/2Z5X/raw.pdb") if l.startswith("HETATM") and l[17:20]=="FAD" and l[21]=="A" and l[16] in " A"])
rows=[]
for f in sorted(glob.glob("04_docking/cand003_redock/2Z5X_seed*_out.pdbqt")):
    seed = re.search(r"seed(\d+)",f)[1]
    mode=None; xs=[]; aff=None
    def flush():
        if xs:
            d=np.linalg.norm(np.array(xs)[:,None]-fad[None],axis=2).min(1)
            rows.append((int(seed),mode,aff,int((d<2.5).sum()),round(float(d.min()),2)))
    for l in open(f):
        if l.startswith("MODEL"): flush(); mode=int(l.split()[1]); xs=[]
        elif l.startswith("REMARK VINA RESULT"): aff=float(l.split()[3])
        elif l.startswith(("ATOM","HETATM")) and l[77:79].strip() not in ("H","HD"): xs.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    flush()
print("seed mode  score  ligand atoms <2.5A of FAD  min dist")
for r in rows: print("%4d %4d %7.2f %10d %14.2f"%r)
ok=[r for r in rows if r[3]==0]
best=min(ok,key=lambda r:r[2]) if ok else None
print("\nbest FAD-clash-free (0 atoms <2.5 A):", best)
