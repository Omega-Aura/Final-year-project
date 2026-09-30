"""Does the best pose sit where a REMOVED heteroatom (FAD, native ligand, waters, ions) is in the crystal?
Distances are from each ligand heavy atom to removed-atom heavy atoms of chain A in the same coordinate frame."""
import json, numpy as np, sys
NAT = {"2V5Z":"mao","2Z5X":"mao","7JXX":"ttbk","7Q8Y":"ttbk"}
NM = {"2V5Z":"MAO-B","2Z5X":"MAO-A","7JXX":"TTBK1","7Q8Y":"TTBK2"}
out = {}
for R in NAT:
    lig = np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in open(f"10_results/cand_003_redock_clean/{R}_{NM[R]}_cand_003_best_pose_ligand.pdb") if l.startswith("HETATM") and l[76:78].strip()!="H"])
    groups = {}
    for l in open(f"03_receptors/{NAT[R]}/{R}/raw.pdb"):
        if l.startswith("HETATM") and l[21]=="A" and l[16] in " A":
            g = l[17:20]; g = "water" if g=="HOH" else g
            groups.setdefault(g, []).append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    res = {}
    for g,xs in groups.items():
        d = np.linalg.norm(lig[:,None,:]-np.array(xs)[None,:,:],axis=2)
        res[g] = dict(n_atoms=len(xs), min_dist=round(float(d.min()),2), lig_atoms_within_2p5=int((d.min(1)<2.5).sum()), lig_atoms_within_3p5=int((d.min(1)<3.5).sum()))
    out[R]=res
    print(R, NM[R]); [print("   %-6s atoms %3d  min dist %5.2f  ligand heavy atoms <2.5A: %d  <3.5A: %d"%(g,v["n_atoms"],v["min_dist"],v["lig_atoms_within_2p5"],v["lig_atoms_within_3p5"])) for g,v in res.items()]
json.dump(out, open("10_results/cand_003_redock_clean/overlap_with_removed_atoms.json","w"), indent=1)
