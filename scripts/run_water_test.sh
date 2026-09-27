#!/usr/bin/env bash
# Wet-vs-dry water-shell test.
#
# Question: does retaining crystallographic waters within 5 A of the NATIVE ligand, as rigid
# receptor atoms, improve or degrade docking for ligands that are not the native one?
#
# The shell is justified in prep_receptor.sh entirely by native-ligand redocking, which is
# circular for cross-docking: those waters are by construction the ones that coexist with the
# native ligand in its own crystal, so redocking it into a pocket moulded around it will tend
# to improve. Measured on the 56 candidates, the wet 7JXX receptor inverts Vina's normal size
# bias (score vs MW r = -0.337 dry-era -> +0.458 wet) and penalises ligands in proportion to
# polarity (delta vs TPSA r = +0.792), which is the signature of steric exclusion rather than
# chemistry.
#
# The dry receptors (03_receptors/<PDB>dry) are the wet ones with ONLY the water atoms
# deleted -- not a re-preparation -- so the waters are the single variable.
#
# Two arms, in power order:
#   1. candidates on 7JXX wet vs dry -- does removing the waters restore the normal
#      score/size relationship? Direct measurement of the artifact on the actual set. (strong)
#   2. references on 7JXX/2V5Z/2Z5X wet vs dry -- correlation against measured IC50.
#      Only 7 compounds carry use_in_correlation=yes and they split 3/2/1/1 across targets,
#      so this arm is weak on its own and is reported with that caveat, not as the decider.
#
# Usage: bash scripts/run_water_test.sh [--dry-run]
set -uo pipefail

DRY=0
for arg in "$@"; do [ "$arg" = "--dry-run" ] && DRY=1; done

if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

SEEDS="11 22 33"

# arm 1 first: it is the one with the power to settle the question
JOBS="
7JXXdry:candidates_56
2V5Zdry:candidates_56
7JXXdry:references
2V5Zdry:references
2Z5Xdry:references
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
    REC="03_receptors/$R/receptor.pdbqt"
    [ -f "$REC" ] || { echo "[skip] $R: no receptor.pdbqt -- run the dry-receptor build first"; continue; }
    N_EXP=$(expected_count "$SET")
    case "$N_EXP" in
        ''|*[!0-9]*) echo "[abort] $R $SET: could not count ligands (got '$N_EXP')" >&2; exit 1 ;;
        0)           echo "[abort] $R $SET: ligand set is empty" >&2; exit 1 ;;
    esac
    for S in $SEEDS; do
        OUT="04_docking/${R}_${SET}_seed${S}"
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
echo "[done] now run: python scripts/analyze_water_test.py"
