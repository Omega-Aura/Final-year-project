#!/usr/bin/env bash
# Usage: prep_receptor.sh PDBID CHAIN LIGRESNAME "COFACTORS"
# Example: prep_receptor.sh 2V5Z A SAG "FAD"
set -euo pipefail

# Add conda env paths to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi
export PATH="/c/Users/sayan/miniforge3/envs/dock/Library/bin:/c/Users/sayan/miniforge3/envs/dock/Scripts:/c/Users/sayan/miniforge3/envs/dock:$PATH"

PDB=$1; CHAIN=$2; LIG=$3; COFACTORS=${4:-""}
D="03_receptors/$PDB"; mkdir -p "$D"

if [ ! -f "$D/raw.pdb" ]; then
    if command -v wget >/dev/null 2>&1; then
        wget -qO "$D/raw.pdb" "https://files.rcsb.org/download/$PDB.pdb"
    else
        curl -s -L -o "$D/raw.pdb" "https://files.rcsb.org/download/$PDB.pdb"
    fi
fi

# 1. native ligand out, as its own file (reference pose for redocking validation)
grep "^HETATM" "$D/raw.pdb" | awk -v l="$LIG" -v c="$CHAIN" \
    'substr($0,18,3)==l && substr($0,22,1)==c' > "$D/native_$LIG.pdb"
obabel "$D/native_$LIG.pdb" -O "$D/native_$LIG.sdf" -h 2>/dev/null || obabel.exe "$D/native_$LIG.pdb" -O "$D/native_$LIG.sdf" -h

# 2. receptor: chosen chain, protein + declared cofactors only. No waters, no cryoprotectant.
KEEP="$COFACTORS"
awk -v c="$CHAIN" -v keep="$KEEP" '
  /^ATOM/   && substr($0,22,1)==c {print}
  /^HETATM/ && substr($0,22,1)==c {
      r=substr($0,18,3); gsub(/ /,"",r)
      if (keep!="" && index(keep,r)>0) print
  }
  /^TER|^END/ {print}
' "$D/raw.pdb" > "$D/clean_noH.pdb"

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

$MK_REC -a --read_pdb "$D/clean_noH.pdb" -o "$D/receptor" -p

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
