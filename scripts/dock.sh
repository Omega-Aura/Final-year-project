#!/usr/bin/env bash
# Usage: dock.sh RECEPTOR LIGANDSET SEED
set -euo pipefail

# Add conda env Library/bin to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

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
    # single combined multi-mol file for this set
    LIGS="02_ligands/pdbqt/${SET}.pdbqt"
elif [ -f "01_smiles/${SET}.csv" ]; then
    # set membership comes from the SMILES CSV's id column (e.g. candidates_56.csv),
    # NOT a blind glob of 02_ligands/pdbqt/*.pdbqt -- that directory holds ligands
    # from every set (native_* references included), so an unscoped glob would
    # silently dock unrelated ligands alongside this set's candidates.
    LIGS=$(python -c "
import pandas as pd
df = pd.read_csv('01_smiles/${SET}.csv')
for i in df['id'].astype(str):
    print(f'02_ligands/pdbqt/{i}.pdbqt')
" | tr -d '\r')
else
    LIGS=$(ls 02_ligands/pdbqt/*.pdbqt)
fi

if command -v vina.exe >/dev/null 2>&1; then
    VINA_CMD="vina.exe"
else
    VINA_CMD="vina"
fi

if command -v mk_export >/dev/null 2>&1; then
    MK_EXPORT="mk_export"
elif command -v mk_export.exe >/dev/null 2>&1; then
    MK_EXPORT="mk_export.exe"
else
    MK_EXPORT="mk_export"
fi

for L in $LIGS; do
  N=$(basename "$L" .pdbqt)
  $VINA_CMD --receptor "03_receptors/$R/receptor.pdbqt" --ligand "$L" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SX --size_y $SY --size_z $SZ \
       --exhaustiveness 32 --num_modes 9 --seed $SEED \
       --out "$OUT/${N}_out.pdbqt" > "$OUT/${N}.log" 2>&1
  # NOTE: must use meeko's own mk_export, not obabel, to convert docked PDBQT -> SDF.
  # obabel doesn't understand meeko's dummy "glue" atoms used for flexible ring bonds,
  # and silently produces `*` wildcard atoms for any ligand needing them (breaks
  # spyrmsd's graph-isomorphism RMSD check and downstream ProLIF/analysis).
  $MK_EXPORT "$OUT/${N}_out.pdbqt" -s "$OUT/${N}_out.sdf" 2>/dev/null || true
  cp "$OUT/${N}_out.sdf" "$OUT/out.sdf" 2>/dev/null || true
done
echo "[ok] $R $SET seed$SEED"
