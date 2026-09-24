#!/usr/bin/env bash
# Usage: prep_receptor.sh PDBID CHAIN LIGRESNAME "COFACTORS" [WATER_SHELL_A]
# Example: prep_receptor.sh 2V5Z A SAG "FAD"
#
# WATER_SHELL_A (default 5.0) retains crystallographic waters whose oxygen lies within
# that distance of the native ligand, and discards all bulk solvent. This is not a
# cosmetic choice -- it is what makes the ATP-site receptors validate. Measured on this
# project's own structures, redocking the native ligand (best of 9 poses, symmetry
# corrected) with and without the pocket shell:
#
#     receptor          dry      +waters
#     7JXX  TTBK1      0.71      0.55      (already passing, improves)
#     4BTK  TTBK1      0.75      0.80      (already passing, native pose reaches rank 1)
#     2V5Z  MAO-B      0.56      0.73      (already passing)
#     7Q8V  TTBK1      4.21      0.70      FAIL -> PASS
#     7Q8Y  TTBK2      5.60      1.45      FAIL -> PASS
#
# 9IV binds both TTBK paralogs through a bridging water; delete it and Vina cannot place
# the ligand at all. 7Q8V needed exactly ONE water to go from 4.21 A to 0.70 A. Retaining
# the shell is applied uniformly to every receptor so that rule 1.5 (any two receptors
# whose scores are compared must be prepared by the same script in the same run) still
# holds across the whole set.
set -euo pipefail

# Add conda env paths to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

PDB=$1; CHAIN=$2; LIG=$3; COFACTORS=${4:-""}; WATER_SHELL=${5:-5.0}
D="03_receptors/$PDB"; mkdir -p "$D"

if [ ! -f "$D/raw.pdb" ]; then
    if command -v wget >/dev/null 2>&1; then
        wget -qO "$D/raw.pdb" "https://files.rcsb.org/download/$PDB.pdb"
    else
        curl -s -L -o "$D/raw.pdb" "https://files.rcsb.org/download/$PDB.pdb"
    fi
fi

# 1. native ligand out, as its own file (reference pose for redocking validation)
# Altloc filter (col 17) is REQUIRED, not cosmetic: a ligand modelled in two alternate
# conformations (7Q8Y's 9IV is A@0.65 / B@0.35) otherwise yields a file with BOTH copies,
# which obabel then perceives as one fragmented 48-atom molecule. Every RMSD measured
# against such a reference is meaningless -- this is what produced the spurious 5.29 A
# "7Q8Y redocking failure" that was previously attributed to the Vina scoring function.
# Keep blank-altloc and altloc A only, matching the --default_altloc A already used for
# the receptor in step 3, so ligand and receptor describe the same single conformer.
grep "^HETATM" "$D/raw.pdb" | awk -v l="$LIG" -v c="$CHAIN" \
    'substr($0,18,3)==l && substr($0,22,1)==c \
     && (substr($0,17,1)==" " || substr($0,17,1)=="A")' > "$D/native_$LIG.pdb"
obabel "$D/native_$LIG.pdb" -O "$D/native_$LIG.sdf" -h 2>/dev/null || obabel.exe "$D/native_$LIG.pdb" -O "$D/native_$LIG.sdf" -h

# 2. receptor: chosen chain, protein + declared cofactors + the pocket water shell.
#    Bulk solvent and cryoprotectant are still discarded; see the WATER_SHELL_A note above.
python - "$D" "$CHAIN" "$COFACTORS" "$WATER_SHELL" "$LIG" <<'PY'
import sys
import numpy as np

d, chain, cofactors, shell, lig = sys.argv[1:6]
shell = float(shell)
keep_res = {r for r in cofactors.replace(",", " ").split() if r}

lig_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in open(f"{d}/native_{lig}.pdb") if l.startswith("HETATM")])

out, nwat = [], 0
for l in open(f"{d}/raw.pdb"):
    if l.startswith("TER") or l.startswith("END"):
        out.append(l)
        continue
    if not l.startswith(("ATOM", "HETATM")) or l[21] != chain:
        continue
    # altloc: keep blank or A, consistently with --default_altloc A used below
    if l[16] not in (" ", "A"):
        continue
    if l.startswith("ATOM"):
        out.append(l)
        continue
    res = l[17:20].strip()
    if res in keep_res:
        out.append(l)
    elif res == "HOH" and shell > 0:
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if np.linalg.norm(lig_xyz - xyz, axis=1).min() < shell:
            out.append(l)
            nwat += 1

open(f"{d}/clean_noH.pdb", "w").writelines(out)
# Meeko gets the polymer + waters only. Cofactors are prepared separately by
# prep_cofactor.py and merged back in, because Meeko's residue templating on a
# covalently-tethered cofactor succeeds or fails depending on the refined bond length
# (see the note in prep_cofactor.py) -- which would make MAO-B and MAO-A, a compared
# pair, differ in preparation. clean_noH.pdb keeps the full record.
poly = [l for l in out if not (l.startswith("HETATM") and l[17:20].strip() in keep_res)]
open(f"{d}/clean_polymer.pdb", "w").writelines(poly)
print(f"{d}: kept {nwat} pocket waters (shell {shell} A), cofactors {sorted(keep_res) or 'none'}")
PY

# 3. pdbqt using Meeko
if command -v mk_prepare_receptor >/dev/null 2>&1; then
    MK_REC="mk_prepare_receptor"
elif command -v mk_prepare_receptor.exe >/dev/null 2>&1; then
    MK_REC="mk_prepare_receptor.exe"
elif command -v mk_prepare_receptor.py >/dev/null 2>&1; then
    MK_REC="mk_prepare_receptor.py"
else
    MK_REC="mk_prepare_receptor"
fi

$MK_REC -a --default_altloc A --read_pdb "$D/clean_polymer.pdb" -o "$D/receptor" -p

# 3a. cofactors: prepared through the ligand toolchain, merged back as rigid receptor atoms
for R in $COFACTORS; do
    python scripts/prep_cofactor.py "$D" "$R"
done

# 3b. GUARD: every declared cofactor must actually survive into receptor.pdbqt.
# `-a/--allow_bad_res` above tells Meeko to DROP any residue it cannot template-match
# instead of raising. That is what we want for ragged protein termini, but it also means a
# cofactor can vanish in silence. It really happened here: 2Z5X's FAD is refined with a
# 1.65 A Cys406 SG--C8M distance, so Meeko perceives the 8alpha-S-cysteinyl linkage as an
# "excess inter-residue bond", template matching fails, and FAD is discarded -- while 2V5Z's
# chemically identical FAD survives only because that crystal models the same bond at 2.31 A,
# too long to be perceived as covalent. The result was a MAO-A receptor with no flavin at all
# lining the active site, used in every MAO-A docking run for weeks without anyone noticing.
# Fail loudly instead.
for R in $COFACTORS; do
    N_IN=$(awk -v r="$R" 'substr($0,1,6)=="HETATM" && substr($0,18,3)==r' "$D/clean_noH.pdb" | wc -l)
    N_OUT=$(awk -v r="$R" '(substr($0,1,4)=="ATOM"||substr($0,1,6)=="HETATM") && substr($0,18,3)==r' "$D/receptor.pdbqt" | wc -l)
    if [ "$N_IN" -gt 0 ] && [ "$N_OUT" -eq 0 ]; then
        echo "[FAIL] $PDB: cofactor $R had $N_IN atoms in clean_noH.pdb but 0 in receptor.pdbqt." >&2
        echo "       Meeko dropped it (see --allow_bad_res). Do NOT dock against this receptor." >&2
        exit 1
    fi
    echo "[cofactor] $R: $N_IN atoms in -> $N_OUT atoms out"
done

# 4. grid box from the native ligand: centroid + bounding box + 8 A padding
python - "$D" "$LIG" <<'PY'
import sys, json, numpy as np
d, lig = sys.argv[1], sys.argv[2]
c = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
              for l in open(f"{d}/native_{lig}.pdb") if l.startswith("HETATM")])
box = {"center": [round(float(x),2) for x in c.mean(0)],
       "size":   [round(float(max(s,18.0)),1) for s in (c.max(0)-c.min(0)+8.0)]}
json.dump(box, open(f"{d}/box.json","w"), indent=2); print(d, box)
PY
echo "[ok] $PDB prepared"
