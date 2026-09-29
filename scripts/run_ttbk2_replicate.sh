#!/usr/bin/env bash
# Velocity replicate of system_TTBK2_p2, then its MM-GBSA and cpptraj analysis.
#
# Closes the last open scientific item from the 2026-09-27 review. system_TTBK2_p2 is the ONLY
# on-pose TTBK2 trajectory in the project, and its -30.95 kcal/mol is the single number the whole
# TTBK1-vs-TTBK2 comparison rests on. It had no error bar of its own: the "no resolvable
# discrimination" verdict compared a 1.45 kcal/mol ddG against TTBK1's 2.08 replicate spread,
# i.e. it judged a two-arm difference using one arm's uncertainty. The other TTBK2 pair
# (pose 1, spread 6.89) is off-pose in both runs and cannot serve -- it describes structures that
# were never the complex. See 06_md/systems/system_TTBK2_p2_r2/PROVENANCE.md.
#
# The replicate carries byte-identical topology and coordinates (sha256 recorded in that
# PROVENANCE.md); the only difference is the random velocity seed, which OpenMM draws itself
# because run_md_system.py sets none. Same mechanism as system_TTBK1_r2.
#
# NOTE: run_md_system.py, NOT run_md_restrained.py. TTBK2 has no FAD to restrain. The MAO
# replicates used the restrained driver, and their unusually tight spreads (0.24-0.26) are partly
# an artifact of that restraint -- which is why this spread had to be measured, not borrowed.
#
# Resumable in both halves: skips MD if final_state.xml exists, skips analysis if
# mmgbsa_results.dat exists. Safe to re-run after an interruption.
#
# Usage: bash scripts/run_ttbk2_replicate.sh
set -u
. "$(dirname "$0")/md_paths.sh"

# Must go through `conda run`, NOT the env's python.exe directly. OpenMM loads its CUDA platform
# from a plugin directory and needs the environment's DLL paths set; invoking the bare interpreter
# silently yields only ['Reference', 'CPU', 'OpenCL'] and the run dies with "There is no registered
# Platform called CUDA".
CONDA_ENV="docking_project"
run_py() { conda run --no-capture-output -n "$CONDA_ENV" python "$@"; }
cd "$(dirname "$0")/.." || exit 1

S="system_TTBK2_p2_r2"
REF="system_TTBK2_p2"
D=$(md_system_dir "$S") || { echo "[abort] no system dir for $S" >&2; exit 1; }

# Refuse to run on inputs that are not the reference's, byte for byte. A replicate whose topology
# or coordinates differ is not a velocity replicate, and its spread would be meaningless.
R=$(md_system_dir "$REF") || { echo "[abort] no system dir for $REF" >&2; exit 1; }
for f in complex.prmtop complex.inpcrd; do
    a=$(sha256sum "$R/$f" | cut -d' ' -f1)
    b=$(sha256sum "$D/$f" | cut -d' ' -f1)
    if [ "$a" != "$b" ]; then
        echo "[abort] $f differs from $REF -- this would not be a velocity replicate" >&2
        exit 1
    fi
done
echo "[ok] inputs byte-identical to $REF"

if [ -f "$D/final_state.xml" ]; then
    echo "[skip md] $S already has final_state.xml"
else
    echo "[start md] $S  $(date '+%Y-%m-%d %H:%M:%S')  (~1.9 h at 131 ns/day)"
    run_py 06_md/run_md_system.py "$D" > "$D/run_md.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "[FAIL md] $S exited $rc -- see $D/run_md.log"
        tail -15 "$D/run_md.log"
        exit $rc
    fi
    echo "[done md]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(tail -1 "$D/production.log" 2>/dev/null)"
fi

if [ -f "$D/mmgbsa_results.dat" ]; then
    echo "[skip analysis] $S"
else
    echo "[start analysis] $S  $(date '+%Y-%m-%d %H:%M:%S')"
    bash scripts/run_mmgbsa.sh "$S" > "$D/mmgbsa_run.log" 2>&1 \
        || { echo "[FAIL mmgbsa] $S"; tail -15 "$D/mmgbsa_run.log"; exit 1; }
    # lig_rmsd only: no rmsf/core_rmsd here. Those exist for the MAO systems because both MAO
    # constructs end in a solvent-exposed C-terminal tail that dominates a whole-protein fit.
    # TTBK2 has no such tail, and system_TTBK2_p2 was analysed with lig_rmsd alone.
    WSLDIR="/mnt/c/Users/aritr/OneDrive/Desktop/Final year project/$D"
    wsl -e bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate mdgbsa && cd '$WSLDIR' && \
        cpptraj -i lig_rmsd.cpptraj > cpptraj_rmsd.log 2>&1" \
        || { echo "[FAIL cpptraj] $S"; exit 1; }
    echo "[done analysis]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(grep 'DELTA TOTAL' "$D/mmgbsa_results.dat")"
fi

echo
echo "[complete] $(date '+%Y-%m-%d %H:%M:%S')"
a=$(grep 'DELTA TOTAL' "$R/mmgbsa_results.dat" 2>/dev/null | awk '{print $3}')
b=$(grep 'DELTA TOTAL' "$D/mmgbsa_results.dat" 2>/dev/null | awk '{print $3}')
if [ -n "$a" ] && [ -n "$b" ]; then
    awk -v a="$a" -v b="$b" 'BEGIN{
        d=a-b; if(d<0)d=-d;
        printf "\nTTBK2 pose 2 replicate spread: %.2f vs %.2f  =  %.2f kcal/mol\n", a, b, d;
        printf "Compare: TTBK ddG is 1.45 kcal/mol (TTBK1 replicate mean -32.40 vs TTBK2 %.2f)\n", a;
        if (d > 1.45) print "  -> spread EXCEEDS the ddG: no resolvable discrimination, now with a symmetric error bar.";
        else print "  -> spread is BELOW the ddG: the verdict needs revisiting against the ~1.3 kcal/mol calibration floor.";
    }'
fi
echo
echo "Then regenerate the summary table and figures:"
echo "  conda run -n $CONDA_ENV python scripts/collect_md_summary.py"
echo "  cd 10_results && conda run -n $CONDA_ENV python make_figures.py"
