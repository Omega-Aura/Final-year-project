#!/usr/bin/env bash
# Minimal-water probe: GATE 1 native redocking on reduced-water receptors.
#
# The 5 A native-ligand shell was shown to impose a large polarity-dependent penalty on
# non-native ligands (delta vs TPSA r = +0.792 wet, +0.049 dry) and to be the sole cause of the
# 7JXX candidate ranking change (dry reproduces the A5 table to max |diff| 0.0000). But dry is
# not universally safe: 7Q8V goes 0.70 -> 4.21 A without its single water and 7Q8Y 1.45 -> 5.60.
#
# So the proposed rule is empirical rather than geometric: start dry and add waters back only
# until the native redock passes. Geometry alone will not do it -- 7Q8V's one required water sits
# 3.78 A from 9IV and bridges nothing, so a "keep bridging waters" rule would delete it.
#
# This probe measures the two unknowns that decide whether the rule is viable:
#   2Z5Xdry   -- does MAO-A still pass with no waters? (never tested; wet is 0.27 A)
#   7Q8Ybrg   -- do the 5 bridging waters suffice, vs 1.45 A with all 10?
# plus two confirmations on the exact files in use (7JXXdry, 2V5Zdry).
#
# 12 dockings, a few minutes. GATE 1 threshold is 2.0 A.
set -uo pipefail

if [ -n "${CONDA_PREFIX:-}" ]; then
    UNIX_PREFIX=$(cygpath -u "$CONDA_PREFIX" 2>/dev/null || echo "$CONDA_PREFIX")
    export PATH="$UNIX_PREFIX/Library/bin:$UNIX_PREFIX/Scripts:$UNIX_PREFIX:$PATH"
fi

# receptor : ligandset (pdbqt basename) : reference sdf basename
PROBES="
7JXXdry:native_VP7:native_VP7
2V5Zdry:native_SAG:native_SAG
2Z5Xdry:native_HRM_ref:native_HRM
7Q8Ybrg:native_9IV:native_9IV
"

for P in $PROBES; do
    R="${P%%:*}"; REST="${P#*:}"; SET="${REST%%:*}"; REF="${REST##*:}"
    [ -f "03_receptors/$R/receptor.pdbqt" ] || { echo "[skip] $R: no receptor"; continue; }
    for S in 11 22 33; do
        bash scripts/dock.sh "$R" "$SET" "$S" >/dev/null 2>&1 \
            || echo "[ERROR] $R $SET seed$S failed" >&2
    done
    NW=$(grep -c 'HOH' "03_receptors/$R/receptor.pdbqt" 2>/dev/null | head -1)
    echo -n "$R (${NW:-0} water atoms): "
    python scripts/rmsd_check.py --ref "03_receptors/$R/${REF}.sdf" \
        --poses "04_docking/${R}_${SET}_seed*/${SET}_out.sdf" 2>/dev/null | tail -1
done
echo "=== PROBE DONE ==="
