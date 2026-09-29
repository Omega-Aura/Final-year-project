#!/usr/bin/env bash
# Where a prepared receptor lives under 03_receptors/.
#
# The 12 receptor directories are grouped by target family (ttbk/, mao/) so that a reader sees
# which structures belong to which comparison, and so the dry/brg/noFAD variants cluster with
# the structure they vary. _cofactors/ and prior_phase/ stay at the top level.
#
# The READ side deliberately carries no mapping at all: receptor_dir resolves a PDB id by glob,
# at either depth, so regrouping a structure cannot make it silently unfindable and a
# half-migrated tree still works. The family rule exists in exactly one place -- receptor_group
# below, used only when CREATING a receptor directory (prep_receptor.sh). One rule, one writer.
#
# Sourced, not executed:  . scripts/receptor_paths.sh

# receptor_group PDB -> target-family directory. Write side only.
receptor_group() {
    case "$1" in
        7JXX*|7Q8V*|7Q8Y*|4BTK*) echo ttbk ;;   # TTBK1 / TTBK2 (4BTK is a TTBK1 structure)
        2V5Z*|2Z5X*)             echo mao  ;;   # MAO-A / MAO-B
        *)                       echo other ;;  # unknown: parked, never silently misfiled
    esac
}

# receptor_dir PDB -> the directory holding that receptor, or empty + status 1.
# Aborts if a PDB id resolves in more than one group, which would mean two copies of a
# receptor and an ambiguous choice of which one every score was computed against.
receptor_dir() {
    local r=$1 hits n
    hits=$(ls -d "03_receptors/$r" "03_receptors"/*/"$r" 2>/dev/null)
    n=$(echo "$hits" | grep -c .)
    case "$n" in
        0) return 1 ;;
        1) echo "$hits" ;;
        *) echo "[abort] receptor '$r' resolves to $n directories:" $hits >&2; exit 1 ;;
    esac
}

# receptor_file PDB RELPATH -> e.g. receptor_file 7JXX receptor.pdbqt
receptor_file() {
    local d
    d=$(receptor_dir "$1") || return 1
    echo "$d/$2"
}
