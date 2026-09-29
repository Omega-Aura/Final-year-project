#!/usr/bin/env bash
# Velocity replicates of the two best on-pose MAO systems, then their MM-GBSA and cpptraj analysis.
#
# Closes the open item left by the 2026-09-27 pose scan: the MAO ddG direction is robust across all
# six on-pose pairings, but its magnitude had no error bar because every MAO arm had exactly one
# trajectory. The only replicate spreads measured in this project are TTBK1's 2.08 and TTBK2's 6.90
# kcal/mol, both 10-30x the SEM that MMPBSA.py reports. Pose spread (0.44-1.13 kcal/mol) is a
# different and smaller source of variance and cannot stand in for it -- see 07_mmgbsa/README.md.
#
# The replicate systems carry byte-identical topology and coordinates (verified by sha256, see each
# directory's PROVENANCE.md); the only difference is the random velocity seed, which OpenMM draws
# itself because run_md_restrained.py sets none. Same mechanism as system_TTBK1_r2.
#
# One GPU job at a time, and resumable in both halves: skips MD for any system that already has
# final_state.xml, and skips analysis for any that already has mmgbsa_results.dat. Safe to re-run
# after an interruption -- which is not hypothetical, a machine restart killed the pose-scan queue
# mid-run on 2026-09-26.
#
# Usage: bash scripts/run_mao_replicate_queue.sh
set -u
. "$(dirname "$0")/md_paths.sh"

# Must go through `conda run`, NOT the env's python.exe directly. OpenMM loads its CUDA platform
# from a plugin directory and needs the environment's DLL paths set; invoking the bare interpreter
# silently yields only ['Reference', 'CPU', 'OpenCL'] and the run dies with "There is no registered
# Platform called CUDA" after the restraints are already set up.
CONDA_ENV="docking_project"
run_py() { conda run --no-capture-output -n "$CONDA_ENV" python "$@"; }
cd "$(dirname "$0")/.." || exit 1

SYSTEMS="system_MAOA_p3_r2 system_MAOB_p2_r2"

for S in $SYSTEMS; do
  if [ -f "$(md_system_dir "$S")/final_state.xml" ]; then
    echo "[skip md] $S already has final_state.xml"
    continue
  fi
  echo "[start md] $S  $(date '+%Y-%m-%d %H:%M:%S')"
  run_py 06_md/run_md_restrained.py "$(md_system_dir "$S")" FAD > "$(md_system_dir "$S")/run_md.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "[FAIL md] $S exited $rc -- see $(md_system_dir "$S")/run_md.log; stopping queue"
    tail -15 "$(md_system_dir "$S")/run_md.log"
    exit $rc
  fi
  echo "[done md]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(tail -1 "$(md_system_dir "$S")/production.log" 2>/dev/null)"
done

for S in $SYSTEMS; do
  if [ -f "$(md_system_dir "$S")/mmgbsa_results.dat" ]; then
    echo "[skip analysis] $S"
    continue
  fi
  echo "[start analysis] $S  $(date '+%Y-%m-%d %H:%M:%S')"
  bash scripts/run_mmgbsa.sh "$S" > "$(md_system_dir "$S")/mmgbsa_run.log" 2>&1 \
    || { echo "[FAIL mmgbsa] $S"; tail -15 "$(md_system_dir "$S")/mmgbsa_run.log"; exit 1; }
  WSLDIR="/mnt/c/Users/aritr/OneDrive/Desktop/Final year project/$(md_system_dir "$S")"
  wsl -e bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate mdgbsa && cd '$WSLDIR' && \
    cpptraj -i lig_rmsd.cpptraj > cpptraj_rmsd.log 2>&1 && \
    cpptraj -i rmsf.cpptraj > rmsf.log 2>&1 && \
    cpptraj -i core_rmsd.cpptraj > core_rmsd.log 2>&1" \
    || { echo "[FAIL cpptraj] $S"; exit 1; }
  echo "[done analysis]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(grep 'DELTA TOTAL' "$(md_system_dir "$S")/mmgbsa_results.dat")"
done

echo "[replicate queue complete] $(date '+%Y-%m-%d %H:%M:%S')"
echo
echo "Replicate spreads (run 1 vs run 2):"
for pair in "system_MAOA_p3 system_MAOA_p3_r2" "system_MAOB_p2 system_MAOB_p2_r2"; do
  set -- $pair
  a=$(grep 'DELTA TOTAL' "$(md_system_dir "$1")/mmgbsa_results.dat" 2>/dev/null | awk '{print $3}')
  b=$(grep 'DELTA TOTAL' "$(md_system_dir "$2")/mmgbsa_results.dat" 2>/dev/null | awk '{print $3}')
  [ -n "$a" ] && [ -n "$b" ] && awk -v a="$a" -v b="$b" -v n="$1" \
    'BEGIN{d=a-b; if(d<0)d=-d; printf "  %-18s %8.2f vs %8.2f  spread %.2f kcal/mol\n", n, a, b, d}'
done
echo
echo "Then regenerate the summary table:"
echo "  conda run -n $CONDA_ENV python scripts/collect_md_summary.py"
