#!/usr/bin/env python
"""Parse Vina logs into a consensus table: best, per-seed mean, inter-seed SD.

Only results docked against the CURRENT receptor are collected. A log older than the
receptor.pdbqt it names was produced by a different receptor -- after the Week 2 rebuild
(pocket water shell, 2Z5X FAD restored, 7Q8Y reference de-duplicated) those numbers are not
comparable with new ones, and the water shell alone moves consensus scores by -1.0 to +0.7
kcal/mol. Mixing them silently produces a table that looks fine and is wrong, so stale logs
are dropped and reported rather than averaged in. --allow-stale restores the old behaviour
for the rare case where you really do want everything on disk.
"""
import argparse, glob, os, re, sys
import pandas as pd

import project_paths          # receptor_file(): 03_receptors/<family>/<PDB>/, resolved by glob

p = argparse.ArgumentParser()
p.add_argument("root")
p.add_argument("-o", required=True)
p.add_argument("--receptors", default="03_receptors",
               help="where receptor.pdbqt files live, for the freshness check")
p.add_argument("--allow-stale", action="store_true",
               help="collect every log regardless of receptor mtime (NOT for production tables)")
p.add_argument("--expect-seeds", type=int, default=3)
a = p.parse_args()

rows = []
stale = {}      # run -> count of logs older than their receptor
no_receptor = set()

# 04_docking/<group>/<run>/*.log since the runs were grouped by purpose; the older flat
# 04_docking/<run>/*.log layout still matches, so this reads either. The run name -- and
# with it receptor, ligandset and seed -- always comes from the log's immediate parent.
for log in sorted(glob.glob(f"{a.root}/*/*.log") + glob.glob(f"{a.root}/*/*/*.log")):
    run = os.path.basename(os.path.dirname(log))
    m = re.match(r"(.+?)_(.+)_seed(\d+)$", run)
    if not m:
        continue
    rec, lset, seed = m.groups()

    if not a.allow_stale:
        rpath = project_paths.receptor_file(rec, root=a.receptors)
        if rpath is None or not os.path.exists(rpath):
            no_receptor.add(rec)
            continue
        if os.path.getmtime(log) <= os.path.getmtime(rpath):
            stale[run] = stale.get(run, 0) + 1
            continue

    scores = [float(x) for x in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)",
                                           open(log).read(), re.M)]
    if scores:
        rows.append({"receptor": rec, "ligandset": lset, "seed": int(seed),
                     "ligand": os.path.basename(log)[:-4], "best": min(scores)})

if stale:
    print(f"[stale] dropped {sum(stale.values())} log(s) older than their receptor:",
          file=sys.stderr)
    for run in sorted(stale):
        print(f"          {run}: {stale[run]}", file=sys.stderr)
    print("        re-dock these or pass --allow-stale if you really want them.",
          file=sys.stderr)
if no_receptor:
    print(f"[warn] no receptor.pdbqt for: {', '.join(sorted(no_receptor))} "
          f"-- their logs were skipped", file=sys.stderr)

df = pd.DataFrame(rows)
if df.empty:
    print("No current docking results found.")
    sys.exit(0)

# A (receptor, ligand) group drawing on more than one ligandset is pooling two different
# preparations of nominally the same compound (e.g. a native_* run and the same ligand in
# references/). Averaging across them is not a consensus, so say so.
mixed = (df.groupby(["receptor", "ligand"])["ligandset"].nunique()
           .loc[lambda s: s > 1])
for (rec, lig) in mixed.index:
    sets = sorted(df[(df.receptor == rec) & (df.ligand == lig)]["ligandset"].unique())
    print(f"[warn] {rec}/{lig} pooled across ligandsets {sets}", file=sys.stderr)

out = (df.groupby(["receptor", "ligand"])["best"]
         .agg(best_overall="min", consensus="mean", seed_sd="std", n_seeds="count")
         .reset_index())

short = out[out.n_seeds != a.expect_seeds]
if not short.empty:
    print(f"[warn] {len(short)} ligand(s) not at {a.expect_seeds} seeds "
          f"(consensus/SD are over fewer runs):", file=sys.stderr)
    for _, r in short.iterrows():
        print(f"          {r.receptor}/{r.ligand}: {int(r.n_seeds)}", file=sys.stderr)

out.to_csv(a.o, index=False)
print(out.to_string(index=False))
