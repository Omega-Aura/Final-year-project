#!/usr/bin/env bash
# Strict receptor preparation for the reference-inhibitor redock, from the RCSB raw files, one receptor at
# a time. Same recipe as the cand_003 redock (see scripts/prep_receptor_strict.py):
#   1 remove water  2 remove every heteroatom  3 rebuild missing atoms (missing loops NOT invented)
#   4 add hydrogens at pH 7.4   5 Meeko: polar-H merge, AutoDock atom types, partial charges -> PDBQT
# Runs in WSL from the project root (envs: recprep = pdbfixer/openmm, ligprep = meeko).
#   bash scripts/refs_redock_prep_receptors.sh
# Validation is a separate step (scripts/refs_validate_receptors.py) and docking refuses to start until
# it passes.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
OUT=04_docking/cand003_redock/refs/receptors
for spec in 2V5Z:mao 2Z5X:mao 7JXX:ttbk 7Q8Y:ttbk; do
  PDB=${spec%%:*}; FAM=${spec##*:}
  mkdir -p "$OUT/$PDB"
  conda activate recprep
  python scripts/prep_receptor_strict.py "$PDB" A "03_receptors/$FAM/$PDB/raw.pdb" "$OUT/$PDB" \
      > "$OUT/$PDB/prep_stdout.txt"
  conda activate ligprep
  mk_prepare_receptor.py -i "$OUT/$PDB/${PDB}_fixed.pdb" -o "$OUT/$PDB/$PDB" -p --default_altloc A \
      > "$OUT/$PDB/meeko_stdout.txt" 2>&1
  echo "prepared $PDB"; ls "$OUT/$PDB"
done
