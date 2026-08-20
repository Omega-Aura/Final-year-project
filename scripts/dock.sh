#!/usr/bin/env bash
# Usage: dock.sh RECEPTOR LIGANDSET SEED
set -euo pipefail

# Add conda env Library/bin to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi
export PATH="/c/Users/sayan/miniforge3/envs/dock/Library/bin:/c/Users/sayan/miniforge3/envs/dock/Scripts:/c/Users/sayan/miniforge3/envs/dock:$PATH"

R=$1; SET=$2; SEED=$3
BOX="03_receptors/$R/box.json"
CX=$(python -c "import json;print(json.load(open('$BOX'))['center'][0])")
CY=$(python -c "import json;print(json.load(open('$BOX'))['center'][1])")
CZ=$(python -c "import json;print(json.load(open('$BOX'))['center'][2])")
SX=$(python -c "import json;print(json.load(open('$BOX'))['size'][0])")
SY=$(python -c "import json;print(json.load(open('$BOX'))['size'][1])")
SZ=$(python -c "import json;print(json.load(open('$BOX'))['size'][2])")
OUT="04_docking/${R}_${SET}_seed${SEED}"; mkdir -p "$OUT"

if [ -f "02_ligands/pdbqt/${SET}.pdbqt" ]; then 
    LIGS="02_ligands/pdbqt/${SET}.pdbqt"
else 
    LIGS=$(ls 02_ligands/pdbqt/*.pdbqt)
fi

if command -v vina.exe >/dev/null 2>&1; then
    VINA_CMD="vina.exe"
else
    VINA_CMD="vina"
fi

for L in $LIGS; do
  N=$(basename "$L" .pdbqt)
  $VINA_CMD --receptor "03_receptors/$R/receptor.pdbqt" --ligand "$L" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SX --size_y $SY --size_z $SZ \
       --exhaustiveness 32 --num_modes 9 --seed $SEED \
       --out "$OUT/${N}_out.pdbqt" > "$OUT/${N}.log" 2>&1
  obabel "$OUT/${N}_out.pdbqt" -O "$OUT/${N}_out.sdf" 2>/dev/null || obabel.exe "$OUT/${N}_out.pdbqt" -O "$OUT/${N}_out.sdf" 2>/dev/null || true
  cp "$OUT/${N}_out.sdf" "$OUT/out.sdf" 2>/dev/null || true
done
echo "[ok] $R $SET seed$SEED"
