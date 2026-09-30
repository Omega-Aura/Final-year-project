import collections, sys
for r in ["2V5Z","2Z5X","7JXX","7Q8Y"]:
    fixed = collections.OrderedDict()
    for l in open(f"{r}_fixed.pdb"):
        if l.startswith("ATOM"): fixed[(l[21],l[22:27])]=l[17:20]
    res = collections.OrderedDict(); types=collections.Counter()
    for l in open(f"{r}_meeko.pdbqt"):
        if l.startswith(("ATOM","HETATM")):
            res.setdefault((l[21],l[22:27]),[]).append(float(l[70:76])); types[l[77:79].strip()]+=1
    dropped=[k for k in fixed if k not in res]
    nonint=[(fixed[k],k[1].strip(),round(sum(v),3)) for k,v in res.items() if abs(sum(v)-round(sum(v)))>0.05]
    tot=sum(sum(v) for v in res.values())
    byres=collections.defaultdict(list)
    for k,v in res.items(): byres[fixed[k]].append(round(sum(v)))
    print(r,"residues in/out",len(fixed),len(res),"dropped",dropped[:5],"| total %.3f"%tot,"| non-integer residues:",nonint[:5])
    print("   charged residues:",{n:collections.Counter(v) for n,v in byres.items() if n in("ARG","LYS","ASP","GLU","HIS","ASN","GLN")})
    print("   types",dict(types))
