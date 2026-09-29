#!/usr/bin/env bash
# Sequential MD queue for the MAO pose scan (MAO-A and MAO-B, docked poses 2 and 3).
#
# Runs one system at a time on purpose. The machine was overwhelmed by concurrent MD earlier in
# this project and restarted mid-run; nothing was lost, but a single GPU job at a time keeps the
# queue survivable and the throughput numbers comparable between systems.
#
# MAO-A order first: the open scientific item is an ON-POSE MAO-A trajectory, without which the
# MAO-A/MAO-B ddG cannot be interpreted (the pose-1 MAO-A run drifted 3.24 A off the docked pose
# while MAO-B held at 2.28 A). MAO-B poses follow to keep the pose scan rank-symmetric.
#
# Usage: bash scripts/run_mao_pose_queue.sh
set -u
. "$(dirname "$0")/md_paths.sh"

# Must go through `conda run`, NOT the env's python.exe directly. OpenMM loads its CUDA platform
# from a plugin directory and needs the environment's DLL paths set; invoking the bare
# interpreter silently yields only ['Reference', 'CPU', 'OpenCL'] and the run dies with
# "There is no registered Platform called CUDA" after the restraints are already set up.
CONDA_ENV="docking_project"
run_py() { conda run --no-capture-output -n "$CONDA_ENV" python "$@"; }
cd "$(dirname "$0")/.." || exit 1

for S in system_MAOA_p2 system_MAOA_p3 system_MAOB_p2 system_MAOB_p3; do
  if [ -f "$(md_system_dir "$S")/final_state.xml" ]; then
    echo "[skip] $S already has final_state.xml -- already completed"
    continue
  fi
  echo "[start] $S  $(date '+%Y-%m-%d %H:%M:%S')"
  run_py 06_md/run_md_restrained.py "$(md_system_dir "$S")" FAD > "$(md_system_dir "$S")/run_md.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "[FAIL] $S exited $rc -- see $(md_system_dir "$S")/run_md.log; stopping queue"
    tail -15 "$(md_system_dir "$S")/run_md.log"
    exit $rc
  fi
  echo "[done]  $S  $(date '+%Y-%m-%d %H:%M:%S')  $(tail -1 "$(md_system_dir "$S")/production.log" 2>/dev/null)"
done
echo "[queue complete] $(date '+%Y-%m-%d %H:%M:%S')"
