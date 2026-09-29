#!/usr/bin/env bash
# Usage: dock.sh RECEPTOR LIGANDSET SEED
set -euo pipefail

# Add conda env Library/bin to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

R=$1; SET=$2; SEED=$3
. "$(dirname "$0")/docking_paths.sh"
. "$(dirname "$0")/receptor_paths.sh"
RDIR=$(receptor_dir "$R") || { echo "[abort] no receptor directory for '$R'" >&2; exit 1; }

# 02_ligands/pdbqt/ is grouped (candidates/, natives/, references/). Resolve a ligand by name
# at EITHER depth: grouped, or loose at the top level as prep_ligands.py wrote them before the
# grouping. Accepting both is deliberate -- a half-migrated tree must not silently resolve to
# nothing, because the fallback below would then dock a different set entirely.
# Prints nothing and returns 1 if the name does not resolve; aborts if it resolves twice.
lig_path() {
    local n=$1 hits
    hits=$(ls "02_ligands/pdbqt/$n.pdbqt" "02_ligands/pdbqt"/*/"$n.pdbqt" 2>/dev/null)
    case $(echo "$hits" | grep -c .) in
        0) return 1 ;;
        1) echo "$hits" ;;
        *) echo "[abort] $n resolves to more than one file:" $hits >&2; exit 1 ;;
    esac
}
BOX="$RDIR/box.json"
CX=$(python -c "import json;print(json.load(open('$BOX'))['center'][0])")
CY=$(python -c "import json;print(json.load(open('$BOX'))['center'][1])")
CZ=$(python -c "import json;print(json.load(open('$BOX'))['center'][2])")
SX=$(python -c "import json;print(json.load(open('$BOX'))['size'][0])")
SY=$(python -c "import json;print(json.load(open('$BOX'))['size'][1])")
SZ=$(python -c "import json;print(json.load(open('$BOX'))['size'][2])")
OUT=$(docking_out "$R" "$SET" "$SEED"); mkdir -p "$OUT"

if SET_FILE=$(lig_path "$SET"); then
    # single combined multi-mol file for this set
    LIGS="$SET_FILE"
elif [ -f "01_smiles/${SET}.csv" ]; then
    # set membership comes from the SMILES CSV's id column (e.g. candidates_56.csv),
    # NOT a blind glob of 02_ligands/pdbqt/*.pdbqt -- that directory holds ligands
    # from every set (native_* references included), so an unscoped glob would
    # silently dock unrelated ligands alongside this set's candidates.
    IDS=$(python -c "
import pandas as pd
df = pd.read_csv('01_smiles/${SET}.csv')
for i in df['id'].astype(str):
    print(i)
" | tr -d '\r')
    LIGS=""
    for I in $IDS; do
        P=$(lig_path "$I") || { echo "[abort] no pdbqt for '$I' (set $SET)" >&2; exit 1; }
        LIGS="$LIGS $P"
    done
else
    LIGS=$(ls 02_ligands/pdbqt/*.pdbqt 02_ligands/pdbqt/*/*.pdbqt 2>/dev/null)
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
  # Per-ligand resume. run_week2_redock.sh only resumes at whole-run granularity, so an
  # interrupted 56-ligand run would otherwise re-dock everything it already finished.
  # A ligand counts as done only if its log is newer than the receptor AND its SDF ends
  # in the SDF record terminator -- the terminator is the last thing written, so a job
  # killed mid-export leaves a truncated SDF that correctly fails this test and is redone.
  if [ -f "$OUT/${N}.log" ] && [ "$OUT/${N}.log" -nt "$RDIR/receptor.pdbqt" ] \
     && [ -s "$OUT/${N}_out.sdf" ] && [ "$(tail -n 1 "$OUT/${N}_out.sdf" | tr -d '\r')" = '$$$$' ]; then
      continue
  fi
  $VINA_CMD --receptor "$RDIR/receptor.pdbqt" --ligand "$L" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SX --size_y $SY --size_z $SZ \
       --exhaustiveness 32 --num_modes 9 --seed $SEED \
       --out "$OUT/${N}_out.pdbqt" > "$OUT/${N}.log" 2>&1
  # NOTE: must use meeko's own mk_export, not obabel, to convert docked PDBQT -> SDF.
  # obabel doesn't understand meeko's dummy "glue" atoms used for flexible ring bonds,
  # and silently produces `*` wildcard atoms for any ligand needing them (breaks
  # spyrmsd's graph-isomorphism RMSD check and downstream ProLIF/analysis).
  # Do NOT silence this. A swallowed mk_export failure leaves the run directory with a
  # valid .log and .pdbqt but no .sdf, and the next rmsd_check.py globs zero poses and
  # reports "BEST None 1000000000.00 A -> FAIL" -- a preparation failure that reads as a
  # docking failure. Retry once, then fail loudly.
  if ! $MK_EXPORT "$OUT/${N}_out.pdbqt" -s "$OUT/${N}_out.sdf" 2>"$OUT/${N}_export.err"; then
      echo "[warn] mk_export failed for $N, retrying" >&2
      $MK_EXPORT "$OUT/${N}_out.pdbqt" -s "$OUT/${N}_out.sdf"
  fi
  [ -s "$OUT/${N}_out.sdf" ] || { echo "[FAIL] no SDF produced for $N" >&2; exit 1; }
  rm -f "$OUT/${N}_export.err"
  cp "$OUT/${N}_out.sdf" "$OUT/out.sdf"
done
echo "[ok] $R $SET seed$SEED"
