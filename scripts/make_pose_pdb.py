#!/usr/bin/env python
"""Write an acpype-template ligand PDB carrying a chosen docked pose's coordinates.

The docking SDF produced by the pipeline has the same atom count and the same atom ORDER as the
acpype template PDB (verified below, and asserted at runtime), so transferring a pose is a
positional coordinate swap. Nothing is re-embedded, so the MD starting geometry is exactly the
geometry the docking produced.

SDF V2000 is parsed directly rather than via RDKit: RDKit stopped loading mid-session under a
Windows Application Control policy (see LOGBOOK 2026-09-25), and the format is fixed-width --
a counts line giving atom and bond counts, then one line per atom with x/y/z and element.

Usage:
  python scripts/make_pose_pdb.py <docked.sdf> <pose_index_1based> <template.pdb> <out.pdb>

Verification: pose 1 must round-trip to the template within --tol (default 0.01 A), which is
what proves the atom order assumption. Pass --pose 1 to run that check explicitly.
"""
import argparse
import math
import sys


def read_sdf_poses(path):
    """Return [[(x,y,z,element), ...], ...], one list per record, in file order."""
    with open(path, 'r', errors='replace') as fh:
        lines = [l.rstrip('\n').rstrip('\r') for l in fh]
    poses, i = [], 0
    while i < len(lines):
        # a record is: title, program, comment, counts, then atoms
        if i + 3 >= len(lines):
            break
        counts = lines[i + 3]
        if 'V2000' not in counts:
            i += 1
            continue
        try:
            natoms = int(counts[0:3])
        except ValueError:
            i += 1
            continue
        atoms = []
        for j in range(i + 4, i + 4 + natoms):
            l = lines[j]
            atoms.append((float(l[0:10]), float(l[10:20]), float(l[20:30]), l[31:34].strip()))
        poses.append(atoms)
        # advance to the record terminator
        while i < len(lines) and lines[i].strip() != '$$$$':
            i += 1
        i += 1
    return poses


def read_template(path):
    """Return (all_lines, [indices of ATOM/HETATM lines])."""
    with open(path, 'r', errors='replace') as fh:
        lines = fh.readlines()
    idx = [n for n, l in enumerate(lines) if l.startswith(('ATOM', 'HETATM'))]
    return lines, idx


def rmsd(a, b):
    n = min(len(a), len(b))
    return math.sqrt(sum((a[i][0] - b[i][0]) ** 2 + (a[i][1] - b[i][1]) ** 2
                         + (a[i][2] - b[i][2]) ** 2 for i in range(n)) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sdf')
    ap.add_argument('pose', type=int, help='1-based pose index')
    ap.add_argument('template')
    ap.add_argument('out')
    ap.add_argument('--tol', type=float, default=0.01,
                    help='max deviation allowed when verifying pose 1 against the template')
    a = ap.parse_args()

    poses = read_sdf_poses(a.sdf)
    if not poses:
        sys.exit('no V2000 records parsed from %s' % a.sdf)
    if not 1 <= a.pose <= len(poses):
        sys.exit('pose %d out of range: %s holds %d poses' % (a.pose, a.sdf, len(poses)))

    lines, idx = read_template(a.template)
    tmpl = []
    for n in idx:
        l = lines[n]
        tmpl.append((float(l[30:38]), float(l[38:46]), float(l[46:54]), l[76:79].strip()))

    sel = poses[a.pose - 1]
    if len(sel) != len(tmpl):
        sys.exit('atom count mismatch: sdf pose has %d, template has %d -- the atom-order '
                 'assumption does not hold, do not proceed' % (len(sel), len(tmpl)))

    # elements must agree position-by-position, else the order assumption is wrong
    bad = [(i, sel[i][3], tmpl[i][3]) for i in range(len(sel))
           if sel[i][3] and tmpl[i][3] and sel[i][3] != tmpl[i][3]]
    if bad:
        for i, e1, e2 in bad[:5]:
            print('  atom %d: sdf %s vs template %s' % (i + 1, e1, e2), file=sys.stderr)
        sys.exit('element sequence mismatch (%d atoms) -- atom orders differ, do not proceed'
                 % len(bad))

    # pose 1 is the pose the template was built from: it must round-trip
    d1 = rmsd(poses[0], tmpl)
    print('pose 1 vs template : %.4f A  (%s)'
          % (d1, 'OK' if d1 <= a.tol else 'FAIL -- atom order assumption broken'))
    if d1 > a.tol:
        sys.exit(1)

    print('pose %d vs pose 1  : %.3f A' % (a.pose, rmsd(sel, poses[0])))

    out = list(lines)
    for k, n in enumerate(idx):
        l = lines[n]
        out[n] = '%s%8.3f%8.3f%8.3f%s' % (l[:30], sel[k][0], sel[k][1], sel[k][2], l[54:])
    with open(a.out, 'w') as fh:
        fh.writelines(out)
    print('wrote %s  (%d atoms, pose %d of %d)' % (a.out, len(idx), a.pose, len(poses)))


if __name__ == '__main__':
    main()
