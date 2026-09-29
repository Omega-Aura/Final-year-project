#!/usr/bin/env bash
# Week 2 re-dock: every production docking run, on the validated (water-shell + cofactor-fixed)
# receptors. Safe to stop and restart -- it skips any run whose logs are already newer than the
# receptor they were docked against, so an interrupted overnight job resumes where it left off.
#
# Why everything has to be re-run: the receptors changed materially this week (5 A pocket water
# shell added; 2Z5X's FAD restored; 7Q8Y's reference ligand de-duplicated). Measured on the
# first compounds through, the water shell moves consensus scores by -1.0 to +0.7 kcal/mol
# depending on the compound, so the old table cannot be mixed with the new one.
#
# Do NOT "speed this up" by lowering --exhaustiveness or setting --cpu. GATE 1 was validated at
# exhaustiveness 32 with Vina's default CPU auto-detect, and Vina's results depend on BOTH:
# changing either decouples the production numbers from the validation that licenses them.
#
# Usage: bash scripts/run_week2_redock.sh [--dry-run]
set -uo pipefail
. "$(dirname "$0")/docking_paths.sh"
. "$(dirname "$0")/receptor_paths.sh"   # docking_out / docking_glob: 04_docking layout

DRY=0
WITH_PARALOG=0
for arg in "$@"; do
    case "$arg" in
        --dry-run)      DRY=1 ;;
        --with-paralog) WITH_PARALOG=1 ;;
    esac
done

if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

SEEDS="11 22 33"

# receptor:ligandset pairs, in priority order. Primary targets first (they replace the
# existing Results tables), then the benchmark/reference sets, then the paralog cross-check.
JOBS="
7JXX:candidates_56
2V5Z:candidates_56
7JXX:references
2V5Z:references
2Z5X:references
"

# Optional (--with-paralog): the 56 candidates on both calibration receptors. This is the
# cross-check on B4 -- if the TTBK1-minus-TTBK2 margin clusters near the +0.784 measured for
# 9IV across chemically diverse compounds, the correction is a property of the receptor pair
# rather than of one ligand. Scientifically valuable but it DOUBLES the job (another 336
# dockings, roughly 15 h), so it is off by default. Run it after the required set lands.
PARALOG_JOBS="
7Q8V:candidates_56
7Q8Y:candidates_56
"
[ "$WITH_PARALOG" = "1" ] && JOBS="$JOBS$PARALOG_JOBS"

expected_count() {
    local set=$1
    if [ -f "01_smiles/${set}.csv" ]; then
        # Deliberately NOT pandas: this count is the resume guard, and it is consulted
        # before the conda env is necessarily on PATH. A bare `python` without pandas
        # returns nothing, N_EXP goes empty, the -ge test errors out, and every run is
        # classed as stale -- i.e. an interrupted 18 h job silently restarts from zero
        # instead of resuming. Row count needs no chemistry; awk is always there.
        awk 'END{print NR-1}' "01_smiles/${set}.csv"
    else
        echo 1
    fi
}

total_run=0; total_skip=0
for J in $JOBS; do
    R="${J%%:*}"; SET="${J##*:}"
    REC=$(receptor_file "$R" receptor.pdbqt) || REC="03_receptors/$R/receptor.pdbqt"
    [ -f "$REC" ] || { echo "[skip] $R: no receptor.pdbqt"; continue; }
    N_EXP=$(expected_count "$SET")
    # Never let an unreadable count degrade into "re-run everything": that is the one
    # failure mode this script exists to prevent.
    case "$N_EXP" in
        ''|*[!0-9]*) echo "[abort] $R $SET: could not count ligands (got '$N_EXP')" >&2; exit 1 ;;
        0)           echo "[abort] $R $SET: ligand set is empty" >&2; exit 1 ;;
    esac
    for S in $SEEDS; do
        OUT=$(docking_out "$R" "$SET" "$S")
        # a run is current if it has the expected number of logs AND every one of them is
        # newer than the receptor file it was supposedly docked against
        N_FRESH=0
        if [ -d "$OUT" ]; then
            N_FRESH=$(find "$OUT" -name '*.log' -newer "$REC" 2>/dev/null | wc -l)
        fi
        if [ "$N_FRESH" -ge "$N_EXP" ]; then
            echo "[current] $R $SET seed$S ($N_FRESH/$N_EXP logs newer than receptor)"
            total_skip=$((total_skip+1))
            continue
        fi
        echo "[run] $R $SET seed$S  ($N_FRESH/$N_EXP current)"
        total_run=$((total_run+1))
        [ "$DRY" = "1" ] && continue
        bash scripts/dock.sh "$R" "$SET" "$S" || echo "[ERROR] $R $SET seed$S failed" >&2
    done
done

echo
echo "runs needed: $total_run   already current: $total_skip"
[ "$DRY" = "1" ] && exit 0

python scripts/collect_results.py 04_docking -o 08_analysis/consensus_week2.csv
echo "[done] wrote 08_analysis/consensus_week2.csv"
