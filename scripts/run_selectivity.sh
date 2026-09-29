#!/usr/bin/env bash
# Selectivity margins for manuscript Table 3.7, on minimal-water receptors.
#
# Protocol (adopted 2026-09-25 after the water-shell test): start dry, add waters back only
# until the native redock passes GATE 1. Validated across the set:
#
#     7JXX 0 waters 0.71 A     4BTK 0 waters 0.75 A     2V5Z 0 waters 0.51 A
#     2Z5X 0 waters 1.02 A     7Q8V 1 water  0.70 A     7Q8Y 5 waters 1.45 A
#
# 7Q8V keeps one water because dry fails outright there (4.21 A); it sits 3.78 A from 9IV and
# bridges nothing, so it is presumably blocking a decoy subpocket. 7Q8Y keeps the five waters
# that bridge ligand and protein -- the other five reproduce nothing (1.45 A either way) and are
# dead volume. The 5 A shell those replace imposed a polarity penalty on non-native ligands
# (score vs TPSA +0.792 wet, +0.049 dry on 7JXX) and was the sole cause of the 7JXX ranking
# change, so the candidate sets must not be scored against it.
#
# The on-target arms (7JXXdry, 2V5Zdry candidates) are already complete; this covers the
# anti-targets and the B4 cross-check pair.
#
#     2Z5Xdry  MAO-A   -> MAO-B vs MAO-A margin
#     7Q8Ybrg  TTBK2   -> TTBK1 vs TTBK2 margin
#     7Q8V     TTBK1   -> B4 cross-check against 7Q8Ybrg, and a second TTBK1 reference point
#
# Usage: bash scripts/run_selectivity.sh [--dry-run]
set -uo pipefail
. "$(dirname "$0")/docking_paths.sh"
. "$(dirname "$0")/receptor_paths.sh"   # docking_out / docking_glob: 04_docking layout

DRY=0
for arg in "$@"; do [ "$arg" = "--dry-run" ] && DRY=1; done

if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

SEEDS="11 22 33"
JOBS="
2Z5Xdry:candidates_56
7Q8Ybrg:candidates_56
7Q8V:candidates_56
"

expected_count() {
    local set=$1
    if [ -f "01_smiles/${set}.csv" ]; then
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
    case "$N_EXP" in
        ''|*[!0-9]*) echo "[abort] $R $SET: could not count ligands (got '$N_EXP')" >&2; exit 1 ;;
        0)           echo "[abort] $R $SET: ligand set is empty" >&2; exit 1 ;;
    esac
    for S in $SEEDS; do
        OUT=$(docking_out "$R" "$SET" "$S")
        N_FRESH=0
        [ -d "$OUT" ] && N_FRESH=$(find "$OUT" -name '*.log' -newer "$REC" 2>/dev/null | wc -l)
        if [ "$N_FRESH" -ge "$N_EXP" ]; then
            echo "[current] $R $SET seed$S ($N_FRESH/$N_EXP)"
            total_skip=$((total_skip+1)); continue
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
echo "[done] next: python scripts/analyze_selectivity.py"
