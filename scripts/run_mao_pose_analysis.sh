#!/usr/bin/env bash
# Post-MD analysis for the MAO pose scan, identical to what was run by hand on system_MAOA_p2:
# strip + MM-GBSA (run_mmgbsa.sh), then ligand/FAD RMSD, per-residue RMSF and core-fit RMSD.
# One system at a time; skips any system that already has mmgbsa_results.dat.
#
# Usage: bash scripts/run_mao_pose_analysis.sh
set -u
. "$(dirname "$0")/md_paths.sh"
cd "$(dirname "$0")/.." || exit 1

for S in system_MAOA_p3 system_MAOB_p2 system_MAOB_p3; do
  if [ -f "$(md_system_dir "$S")/mmgbsa_results.dat" ]; then echo "[skip] $S"; continue; fi
  echo "[start] $S  $(date '+%Y-%m-%d %H:%M:%S')"
  bash scripts/run_mmgbsa.sh "$S" > "$(md_system_dir "$S")/mmgbsa_run.log" 2>&1 || { echo "[FAIL] $S mmgbsa"; tail -15 "$(md_system_dir "$S")/mmgbsa_run.log"; exit 1; }
  WSLDIR="/mnt/c/Users/aritr/OneDrive/Desktop/Final year project/$(md_system_dir "$S")"
  wsl -e bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate mdgbsa && cd '$WSLDIR' && \
    cpptraj -i lig_rmsd.cpptraj > cpptraj_rmsd.log 2>&1 && \
    cpptraj -i rmsf.cpptraj > rmsf.log 2>&1 && \
    cpptraj -i core_rmsd.cpptraj > core_rmsd.log 2>&1" || { echo "[FAIL] $S cpptraj"; exit 1; }
  echo "[done]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(grep 'DELTA TOTAL' "$(md_system_dir "$S")/mmgbsa_results.dat")"
done
echo "[analysis complete] $(date '+%Y-%m-%d %H:%M:%S')"
