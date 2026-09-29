#!/usr/bin/env bash
# Single-trajectory GB MM-GBSA for a prepared complex, run through the WSL mdgbsa env
# (AmberTools lives there; OpenMM/CUDA runs natively on Windows -- see LOGBOOK 2026-08-19).
#
# Same protocol as the TTBK1 baseline in 06_md/system: igb=5, saltcon=0.150, every 5th frame.
# Kept identical on purpose -- the whole point of the TTBK2 run is a like-for-like delta-G
# comparison against TTBK1, so the energy model and frame stride must not drift between them.
#
# Usage: bash scripts/run_mmgbsa.sh <system_dir_name>      e.g. system_TTBK2
set -euo pipefail
. "$(dirname "$0")/md_paths.sh"

SYS="${1:?usage: run_mmgbsa.sh <system_dir_name under 06_md/>}"
SYSDIR=$(md_system_dir "$SYS") || { echo "[abort] no MD system dir for '$SYS'" >&2; exit 1; }
WSLDIR="/mnt/c/Users/aritr/OneDrive/Desktop/Final year project/$SYSDIR"

wsl -e bash -lc "
source ~/miniconda3/etc/profile.d/conda.sh && conda activate mdgbsa
cd '$WSLDIR' || exit 1
set -e

echo '--- stripping solvent from trajectory ---'
cpptraj -i strip_traj.cpptraj

echo '--- building GB topologies ---'
ante-MMPBSA.py -p complex.prmtop -c gb_complex.prmtop -r gb_receptor.prmtop \
               -l gb_ligand.prmtop -s ':WAT,Na+,Cl-' -n ':MOL' --radii=mbondi2

echo '--- MM-GBSA ---'
# No -sp here. -sp is the SOLVATED topology and is only correct when the trajectory still has
# solvent in it. This trajectory has already been stripped by cpptraj above, so it matches -cp
# (4833 atoms) and passing the 48158-atom solvated topology just aborts with an atom-count
# mismatch.
MMPBSA.py -O -i mmpbsa.in -o mmgbsa_results.dat \
          -cp gb_complex.prmtop -rp gb_receptor.prmtop -lp gb_ligand.prmtop \
          -y production_full_stripped.nc

echo '--- result ---'
grep -A4 'DELTA TOTAL' mmgbsa_results.dat | head -6
"
