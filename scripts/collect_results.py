#!/usr/bin/env python
"""Parse Vina logs into a consensus table: best, per-seed mean, inter-seed SD."""
import argparse, glob, os, re
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("root")
p.add_argument("-o", required=True)
a = p.parse_args()

rows = []
for log in glob.glob(f"{a.root}/*/*.log"):
    run = os.path.basename(os.path.dirname(log))
    m = re.match(r"(.+?)_(.+)_seed(\d+)$", run)
    if not m:
        continue
    rec, lset, seed = m.groups()
    scores = [float(x) for x in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)",
                                           open(log).read(), re.M)]
    if scores:
        rows.append({"receptor": rec, "ligandset": lset, "seed": int(seed),
                     "ligand": os.path.basename(log)[:-4], "best": min(scores)})

df = pd.DataFrame(rows)
if not df.empty:
    out = (df.groupby(["receptor", "ligand"])["best"]
             .agg(best_overall="min", consensus="mean", seed_sd="std", n_seeds="count")
             .reset_index())
    out.to_csv(a.o, index=False)
    print(out.to_string(index=False))
else:
    print("No docking results found.")
