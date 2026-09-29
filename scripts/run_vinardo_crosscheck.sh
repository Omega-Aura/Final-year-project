#!/usr/bin/env bash
# Independent scoring-function cross-check on the VALIDATED receptors.
#
# Closes the gap stated in manuscript section 3.3 and WORKFLOW Phase 3. The only Vinardo data in
# this project is prior-phase, run on 4NFM/2V5Z -- the earlier receptor set, with an apo TTBK1
# structure that carried no passing redocking validation. So the published caveat ("docking rank
# order does not survive a change of scoring function", r = 0.56 at TTBK1 and -0.28 at MAO-B) is
# a general warning from superseded receptors, NOT a cross-check of the rankings actually
# reported. This run makes it one.
#
# Design: all 56 candidates, both validated water-symmetric on-targets, 3 seeds, at the SAME
# exhaustiveness/num_modes/box/seeds as the production vina run. Only the scoring function
# differs -- that is the whole point, and it is the same arm-symmetry rule this project applies
# everywhere else. It also improves on the prior-phase cross-check, which covered only the top 15.
#
# Output goes to 04_docking/crosscheck/vinardo/, one directory level deeper than the globs in
# collect_results.py / analyze_selectivity.py / analyze_water_test.py, so Vinardo scores can never
# be pooled into a vina consensus table or averaged into a vina selectivity margin. Vinardo runs
# on a different scale (roughly -5 to -9.5 where vina gives -7 to -12); mixing them on one axis
# would be the same class of error as the 0-vs-5 water shell.
#
# DOCK_SDF=0: this is a SCORE cross-check, so no poses are exported. That is also the only way it
# can run at present -- mk_export.exe is blocked by a Device Guard policy on this machine as of
# 2026-09-28, so any run needing an SDF fails at the export step.
#
# Resumable: dock.sh skips any ligand whose log is already newer than the receptor.
#
# Usage: bash scripts/run_vinardo_crosscheck.sh [--dry-run]
set -uo pipefail
. "$(dirname "$0")/receptor_paths.sh"
cd "$(dirname "$0")/.." || exit 1

DRY=0
for a in "$@"; do
    case "$a" in
        --dry-run) DRY=1 ;;
        *) echo "unknown option: $a" >&2; exit 2 ;;
    esac
done

RECEPTORS="7JXXdry 2V5Zdry"     # TTBK1 and MAO-B, the validated water-symmetric on-target pair
SET="candidates_56"
SEEDS="11 22 33"
SCORING="vinardo"
export DOCK_SDF=0

N_EXP=$(python -c "
import pandas as pd
print(len(pd.read_csv('01_smiles/${SET}.csv')))
" 2>/dev/null | tr -d '\r')
case "$N_EXP" in
    ''|*[!0-9]*) echo "[abort] could not count ligands in 01_smiles/${SET}.csv (got '$N_EXP')" >&2; exit 1 ;;
    0)           echo "[abort] ligand set is empty" >&2; exit 1 ;;
esac
echo "ligand set $SET: $N_EXP ligands"

total_run=0; total_skip=0
for R in $RECEPTORS; do
    RDIR=$(receptor_dir "$R") || { echo "[skip] $R: no receptor directory"; continue; }
    REC="$RDIR/receptor.pdbqt"
    [ -f "$REC" ] || { echo "[skip] $R: no receptor.pdbqt"; continue; }
    for S in $SEEDS; do
        OUT="04_docking/crosscheck/$SCORING/${R}_${SET}_seed${S}"
        N_FRESH=0
        [ -d "$OUT" ] && N_FRESH=$(find "$OUT" -name '*.log' -newer "$REC" 2>/dev/null | wc -l)
        if [ "$N_FRESH" -ge "$N_EXP" ]; then
            echo "[current] $R $SET seed$S ($N_FRESH/$N_EXP logs newer than receptor)"
            total_skip=$((total_skip+1))
            continue
        fi
        echo "[run] $R $SET seed$S  ($N_FRESH/$N_EXP current)  $(date '+%H:%M:%S')"
        total_run=$((total_run+1))
        [ "$DRY" = "1" ] && continue
        bash scripts/dock.sh "$R" "$SET" "$S" "$SCORING" \
            || { echo "[FAIL] $R $SET seed$S" >&2; exit 1; }
    done
done

echo
echo "runs needed: $total_run   already current: $total_skip"
[ "$DRY" = "1" ] && exit 0

echo "[complete] $(date '+%Y-%m-%d %H:%M:%S')"
echo
echo "Then analyse:"
echo "  conda run -n docking_project python scripts/analyze_vinardo_crosscheck.py"
