#!/usr/bin/env bash
# Reference-inhibitor docking, one run per inhibitor on its OWN target receptor, prepared the same way as the
# cand_003 redock. Receptors: refs/receptors/<PDB>/<PDB>.pdbqt; ligands: refs/ligands/<id>/<id>.pdbqt; both
# prepared and validated separately (scripts/refs_redock_prep_receptors.sh, refs_validate_receptors.py,
# refs_prep_ligands.py). The script refuses to start unless the validators pass.
#
# Active site known (native-ligand box, 03_receptors/*/box.json, 18 A cube) -> targeted docking,
# exhaustiveness 32, 9 modes, ONE run per pair (--seed 11).
cd "$(dirname "$0")"
PY="/c/Users/aritr/.conda/envs/docking_project/python.exe"
"$PY" ../../scripts/refs_validate_receptors.py > refs/receptors/validation_stdout.txt || { echo "receptor validation failed"; exit 1; }
for L in refs/ligands/*/report.json; do
  "$PY" -c "import json,sys;sys.exit(1 if json.load(open('$L'))['failures'] else 0)" || { echo "ligand validation failed: $L"; exit 1; }
done
VINA="/c/Program Files (x86)/The Scripps Research Institute/Vina/vina.exe"
mkdir -p refs/docking
declare -A DIR=( [2V5Z]=mao [2Z5X]=mao [7JXX]=ttbk [7Q8Y]=ttbk )
declare -A LIGS=(
  [2V5Z]="safinamide selegiline rasagiline lazabemide isatin"   # MAO-B inhibitors
  [2Z5X]="kaempferol quercetin clorgyline harmine"               # MAO-A (measured / native)
  [7JXX]="VP7 DTQ 9IV"                                           # TTBK1
  [7Q8Y]="9IV" )                                                 # TTBK2
for R in 2V5Z 2Z5X 7JXX 7Q8Y; do
  B=../../03_receptors/${DIR[$R]}/$R/box.json
  read -r CX CY CZ SX SY SZ < <("$PY" -c "import json;b=json.load(open('$B'));print(*b['center'],*b['size'])" | tr -d '\r')
  for L in ${LIGS[$R]}; do
    "$VINA" --receptor refs/receptors/$R/$R.pdbqt --ligand refs/ligands/$L/$L.pdbqt \
      --center_x $CX --center_y $CY --center_z $CZ --size_x $SX --size_y $SY --size_z $SZ \
      --scoring vina --exhaustiveness 32 --num_modes 9 --seed 11 \
      --out refs/docking/${R}_${L}_seed11_out.pdbqt > refs/docking/${R}_${L}_seed11.log 2>&1
    echo "$R $L: $(grep -m1 -E '^ +1 ' refs/docking/${R}_${L}_seed11.log)"
  done
done
echo ALLDONE
