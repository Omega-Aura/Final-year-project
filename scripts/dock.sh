#!/usr/bin/env bash
# Usage: dock.sh RECEPTOR LIGANDSET SEED [SCORING]
#
# SCORING defaults to vina, the production scoring function. Anything else (vinardo, ad4) is an
# independent cross-check and its output goes to a SEPARATE tree that the vina collectors cannot
# see: 04_docking/crosscheck/<scoring>/<run>/ sits one level deeper than the globs in
# collect_results.py, analyze_selectivity.py and analyze_water_test.py, so a vinardo score can
# never be pooled into a vina consensus table or averaged into a vina selectivity margin. Two
# scoring functions on one axis is the same class of error as the water-shell asymmetry -- the
# numbers look comparable and are not (Vinardo runs on a different scale entirely).
set -euo pipefail

# Add conda env Library/bin to PATH
if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

R=$1; SET=$2; SEED=$3; SCORING=${4:-vina}

# DOCK_SDF=0 skips the PDBQT -> SDF conversion. Poses are needed for RMSD validation and
# interaction profiling; a scoring-function cross-check needs only the affinities in the .log, so
# it can opt out. Default is 1, because for a production run a missing SDF is a real failure and
# must stay fatal (see the mk_export note below).
#
# It is also the only way to dock at all on this machine as of 2026-09-28: mk_export.exe is
# blocked by a Device Guard policy, so any run that requires an SDF fails at the export step.
DOCK_SDF=${DOCK_SDF:-1}
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
if [ "$SCORING" = "vina" ]; then
    OUT=$(docking_out "$R" "$SET" "$SEED")
else
    OUT="04_docking/crosscheck/$SCORING/${R}_${SET}_seed${SEED}"
fi
mkdir -p "$OUT"

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

# Check the .exe FIRST. On Windows only mk_export.exe exists on disk, but `command -v mk_export`
# still resolves it and prints the name without the extension; selecting that name then fails to
# exec with a misleading "Permission denied". Checking .exe first avoids picking a name that is
# not a real file.
if command -v mk_export.exe >/dev/null 2>&1; then
    MK_EXPORT="mk_export.exe"
elif command -v mk_export >/dev/null 2>&1; then
    MK_EXPORT="mk_export"
else
    MK_EXPORT="mk_export"
fi

# Resolve the export route ONCE, here, instead of discovering the answer after the first ligand has
# already been docked. A 56-ligand run would otherwise burn ~15 minutes per seed before reporting
# it, and a 6-run batch nearly an hour and a half.
#
# Context: on 2026-09-28 Windows Smart App Control began blocking mk_export.exe (and rdchem.pyd)
# under an enforced Code Integrity policy. Vina itself is unaffected, so docking "works" and only
# the pose conversion fails -- exactly the shape of failure worth catching early. Smart App Control
# has no allowlist, so a Defender path exclusion cannot reach it; WSL2 runs ELF binaries and is not
# governed by it. See LOGBOOK entry M and scripts/wsl_run.sh.
if [ "$DOCK_SDF" != "0" ]; then
    if $MK_EXPORT --help >/dev/null 2>&1; then
        :                                   # native export works, nothing to do
    elif wsl -e bash -lc "source ~/miniconda3/etc/profile.d/conda.sh \
             && conda activate ${LIGPREP_ENV:-ligprep} \
             && command -v mk_export.py" >/dev/null 2>&1; then
        # Route the conversion through WSL. $MK_EXPORT is expanded unquoted at the call site, so a
        # multi-word command word splits correctly; wsl_run.sh cd's to the project root inside WSL,
        # and dock.sh also runs from there, so the relative paths passed to it still resolve.
        MK_EXPORT="bash scripts/wsl_run.sh mk_export"
        echo "[note] native mk_export is blocked; exporting SDFs through WSL (${LIGPREP_ENV:-ligprep})" >&2
    else
        echo "[abort] '$MK_EXPORT' cannot execute and no WSL fallback is available, so no SDF" >&2
        echo "        can be written. Vina would still dock; only the PDBQT -> SDF conversion" >&2
        echo "        fails, which breaks rmsd_check.py and ProLIF downstream." >&2
        echo "        Known cause on this machine: Smart App Control blocks mk_export.exe" >&2
        echo "        (first seen 2026-09-28). Fix: build the WSL env, then re-run --" >&2
        echo "          conda create -n ligprep --override-channels -c conda-forge \\" >&2
        echo "            python=3.11 'rdkit==2025.09.6' 'meeko==0.7.1' openbabel" >&2
        echo "          bash scripts/wsl_run.sh --check" >&2
        echo "        If this run does not need poses -- a scoring-function cross-check, say --" >&2
        echo "        re-run with DOCK_SDF=0 and only the .log affinities will be produced." >&2
        exit 1
    fi
fi

for L in $LIGS; do
  N=$(basename "$L" .pdbqt)
  # Per-ligand resume. run_week2_redock.sh only resumes at whole-run granularity, so an
  # interrupted 56-ligand run would otherwise re-dock everything it already finished.
  # A ligand counts as done only if its log is newer than the receptor AND its SDF ends
  # in the SDF record terminator -- the terminator is the last thing written, so a job
  # killed mid-export leaves a truncated SDF that correctly fails this test and is redone.
  if [ -f "$OUT/${N}.log" ] && [ "$OUT/${N}.log" -nt "$RDIR/receptor.pdbqt" ]; then
      if [ "$DOCK_SDF" = "0" ]; then
          continue          # score-only run: a fresh log is the whole product
      fi
      if [ -s "$OUT/${N}_out.sdf" ] && [ "$(tail -n 1 "$OUT/${N}_out.sdf" | tr -d '\r')" = '$$$$' ]; then
          continue
      fi
  fi
  # --scoring is passed explicitly even for vina, so the log records which function produced
  # it. Everything else is identical across scoring functions on purpose: the cross-check is only
  # interpretable if the search settings, box and seed are the same and ONLY the function differs.
  $VINA_CMD --receptor "$RDIR/receptor.pdbqt" --ligand "$L" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SX --size_y $SY --size_z $SZ \
       --scoring $SCORING \
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
  if [ "$DOCK_SDF" = "0" ]; then
      continue
  fi
  if ! $MK_EXPORT "$OUT/${N}_out.pdbqt" -s "$OUT/${N}_out.sdf" 2>"$OUT/${N}_export.err"; then
      echo "[warn] mk_export failed for $N, retrying" >&2
      $MK_EXPORT "$OUT/${N}_out.pdbqt" -s "$OUT/${N}_out.sdf"
  fi
  [ -s "$OUT/${N}_out.sdf" ] || { echo "[FAIL] no SDF produced for $N" >&2; exit 1; }
  rm -f "$OUT/${N}_export.err"
  cp "$OUT/${N}_out.sdf" "$OUT/out.sdf"
done
echo "[ok] $R $SET seed$SEED ($SCORING)"
