#!/usr/bin/env bash
# Where a docking run's output directory lives under 04_docking/.
#
# 04_docking/ held 78 sibling run directories named ${RECEPTOR}_${LIGANDSET}_seed${SEED},
# which is 85% of every tracked file in the project and unreadable as a listing. They are
# now grouped by purpose. The run directory name itself is UNCHANGED -- every parser still
# derives receptor/ligandset/seed from `basename $(dirname log)`, and the readers simply
# glob one extra level (04_docking/*/${R}_${SET}_seed*/). Only the write side needs to know
# which group a set belongs to, and that rule lives here and nowhere else: four scripts
# create run directories, and a copy of this mapping in each of them would drift.
#
# Sourced, not executed:  . scripts/docking_paths.sh

# docking_group LIGANDSET -> group directory name
docking_group() {
    case "$1" in
        candidates_56) echo candidates    ;;  # the production screen
        references)    echo references    ;;  # the known-inhibitor benchmark set
        native_*)      echo native_redock ;;  # redocking validation, native co-crystal ligand
        cand_*)        echo controls      ;;  # single-ligand controls (noFAD, bridging-water)
        *)             echo other         ;;  # unknown set: parked, not silently misfiled
    esac
}

# docking_out RECEPTOR LIGANDSET SEED -> full output path
docking_out() {
    echo "04_docking/$(docking_group "$2")/${1}_${2}_seed${3}"
}

# docking_glob RECEPTOR LIGANDSET -> glob matching that run's seeds in any group.
# Group-agnostic on purpose, so a set that gets reclassified is still found.
docking_glob() {
    echo "04_docking/*/${1}_${2}_seed*"
}
