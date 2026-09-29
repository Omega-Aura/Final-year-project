#!/usr/bin/env bash
# Where an MD system lives under 06_md/.
#
# The 16 system directories moved into 06_md/systems/ so that 06_md/'s own listing shows the
# things a reader chooses between -- the run scripts, params/, queue_logs/, prior_phase/ -- rather
# than 16 near-identical system names. They are NOT split further by target: the names already
# sort into target clusters (system_MAOA*, system_MAOB*, system_TTBK1*, system_TTBK2*), and every
# extra level has to be paid for again in three WSL absolute paths and the summary table.
#
# Resolution is by glob at either depth, so a half-migrated tree still works and moving a system
# cannot make it silently unfindable -- which matters more here than anywhere else in the project,
# because a "missing" system is silently dropped from the summary table rather than erroring.
#
# Sourced, not executed:  . scripts/md_paths.sh

# md_system_dir NAME -> repo-relative directory for that system, or empty + status 1.
# Aborts on an ambiguous hit: two directories for one system name would make it undefined
# which trajectory an energy came from.
md_system_dir() {
    local s=$1 hits n
    hits=$(ls -d "06_md/systems/$s" "06_md/$s" 2>/dev/null)
    n=$(echo "$hits" | grep -c .)
    case "$n" in
        0) return 1 ;;
        1) echo "$hits" ;;
        *) echo "[abort] MD system '$s' resolves to $n directories:" $hits >&2; exit 1 ;;
    esac
}
