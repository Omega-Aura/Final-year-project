#!/usr/bin/env python
"""Rebuild 07_mmgbsa/md_mmgbsa_summary.csv from the primary per-system files.

The summary table is the thing that gets read and quoted, so it must never be hand-edited:
a typo in a transcribed number is invisible and survives review. Everything here is parsed
out of the files cpptraj and MMPBSA.py actually wrote -- `lig_rmsd.dat`, `prot_rmsd.dat`,
`fad_rmsd.dat`, the core-fit variants, and `mmgbsa_results*.dat`.

Stability is judged on the MEAN OF THE FINAL 100 FRAMES, not the whole-run mean: a ligand that
leaves late still shows a low whole-run average. See 06_md/README.md.

The MAO systems additionally get core-fit numbers, because both MAO constructs end in a
solvent-exposed C-terminal tail (MAO-A runs to residue 513, MAO-B to 499) that flails without a
membrane and inflates whole-protein backbone RMSD to as much as 3.44 A while the core is at
1.3 A. Fitting on `:1-496` and re-measuring the ligand separates real ligand motion from the
tail dragging the global fit around.

Run through `conda run`, NOT a bare `python`: the interpreter first on PATH is MGLTools' Python
2.7 (it is there for the AutoDock prep scripts) and this file is Python 3.

Usage: conda run -n docking_project python scripts/collect_md_summary.py          # writes the CSV
       conda run -n docking_project python scripts/collect_md_summary.py --check  # 1 if stale
"""
import argparse
import csv
import io
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "07_mmgbsa", "md_mmgbsa_summary.csv")

# (system_dir, target, structure, start). Order is the reporting order, not disk order.
SYSTEMS = [
    ("06_md/system",            "TTBK1", "7JXX", "pose 1, run 1"),
    ("06_md/system_TTBK1_r2",   "TTBK1", "7JXX", "pose 1, run 2"),
    ("06_md/system_TTBK1_p2",   "TTBK1", "7JXX", "pose 2"),
    ("06_md/system_TTBK1_p3",   "TTBK1", "7JXX", "pose 3"),
    ("06_md/system_TTBK2",      "TTBK2", "7Q8Y", "pose 1, run 1"),
    ("06_md/system_TTBK2m",     "TTBK2", "7Q8Y", "pose 1, run 2"),
    ("06_md/system_TTBK2_p2",   "TTBK2", "7Q8Y", "pose 2"),
    ("06_md/system_TTBK2_p3",   "TTBK2", "7Q8Y", "pose 3"),
    ("06_md/system_MAOA",       "MAO-A", "2Z5X", "pose 1"),
    ("06_md/system_MAOA_p2",    "MAO-A", "2Z5X", "pose 2"),
    ("06_md/system_MAOA_p3",    "MAO-A", "2Z5X", "pose 3"),
    ("06_md/system_MAOB",       "MAO-B", "2V5Z", "pose 1"),
    ("06_md/system_MAOB_p2",    "MAO-B", "2V5Z", "pose 2"),
    ("06_md/system_MAOB_p3",    "MAO-B", "2V5Z", "pose 3"),
]

# Verdict thresholds on the last-100-frame ligand RMSD, in Angstrom. The <=2.3 / >=3.2 bands are
# where the original ten runs fell; the MAO pose scan then produced two runs in the gap between
# them (2.43 and 2.57 A), which is why the middle band exists and is named separately instead of
# being rounded into "stable" or "drifts".
STABLE_MAX = 2.30
LOOSE_MAX = 3.20
DISSOCIATES_MIN = 5.00


def verdict(last100):
    if last100 is None:
        return ""
    if last100 <= STABLE_MAX:
        return "stable"
    if last100 <= LOOSE_MAX:
        return "holds (loose)"
    if last100 < DISSOCIATES_MIN:
        return "drifts"
    return "dissociates"


def read_series(path):
    """Column 2 of a cpptraj .dat file, skipping its '#Frame ...' header."""
    if not os.path.exists(path):
        return None
    vals = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            vals.append(float(line.split()[1]))
    return vals or None


def last100(vals):
    return None if vals is None else sum(vals[-100:]) / len(vals[-100:])


def read_dg(sysdir):
    """(dG, SD, SEM, source_filename) from the LAST 'DELTA TOTAL' line of the MM-GBSA output.

    'mmgbsa_results_full.dat' wins over 'mmgbsa_results.dat' where both exist: in system/ the
    plain file is the same run truncated at 1 ns, kept only to show the estimate converged.
    The last DELTA TOTAL is the one under 'Differences (Complex - Receptor - Ligand)'; earlier
    ones belong to the per-moiety decompositions.
    """
    for name in ("mmgbsa_results_full.dat", "mmgbsa_results.dat"):
        path = os.path.join(sysdir, name)
        if not os.path.exists(path):
            continue
        hit = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("DELTA TOTAL"):
                    hit = line.split()
        if hit:
            return float(hit[2]), float(hit[3]), float(hit[4]), name
    return None, None, None, "not computed (off-pose trajectory)"


def fmt(v, nd=2):
    return "" if v is None else f"{v:.{nd}f}"


def build():
    rows = []
    for rel, target, structure, start in SYSTEMS:
        sysdir = os.path.join(ROOT, rel.replace("06_md/", "06_md" + os.sep))
        lig = read_series(os.path.join(sysdir, "lig_rmsd.dat"))
        if lig is None:
            sys.exit(f"missing primary data: {rel}/lig_rmsd.dat")
        dg, sd, sem, source = read_dg(sysdir)
        rows.append({
            "system_dir": rel,
            "target": target,
            "structure": structure,
            "start": start,
            "frames": len(lig),
            "lig_rmsd_mean": fmt(sum(lig) / len(lig)),
            "lig_rmsd_max": fmt(max(lig)),
            "lig_rmsd_last100": fmt(last100(lig)),
            "lig_corefit_last100": fmt(last100(read_series(os.path.join(sysdir, "lig_corefit_rmsd.dat")))),
            "pose_verdict": verdict(last100(lig)),
            "prot_rmsd_last100": fmt(last100(read_series(os.path.join(sysdir, "prot_rmsd.dat")))),
            "prot_core_rmsd_last100": fmt(last100(read_series(os.path.join(sysdir, "prot_core_rmsd.dat")))),
            "fad_rmsd_last100": fmt(last100(read_series(os.path.join(sysdir, "fad_rmsd.dat")))),
            "dG_kcal_mol": fmt(dg),
            "dG_SD": fmt(sd),
            "dG_SEM": fmt(sem),
            "gbsa_source": source,
        })
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue(), rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the committed CSV differs from the primary data")
    a = ap.parse_args()

    text, rows = build()
    if a.check:
        on_disk = open(OUT).read() if os.path.exists(OUT) else ""
        if on_disk != text:
            print(f"STALE: {os.path.relpath(OUT, ROOT)} does not match the primary .dat files")
            return 1
        print(f"up to date: {os.path.relpath(OUT, ROOT)} ({len(rows)} systems)")
        return 0

    with open(OUT, "w", newline="") as fh:
        fh.write(text)
    computed = sum(1 for r in rows if r["dG_kcal_mol"])
    print(f"wrote {os.path.relpath(OUT, ROOT)}: {len(rows)} systems, {computed} with MM-GBSA")
    return 0


if __name__ == "__main__":
    sys.exit(main())
