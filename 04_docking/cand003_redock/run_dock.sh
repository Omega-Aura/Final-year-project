#!/usr/bin/env bash
# AutoDock Vina 1.2.7, cand_003 vs four strictly-cleaned receptors. The active site is known (native
# ligand box from 03_receptors/*/box.json), so this is targeted docking, not blind: exhaustiveness 32.
cd "$(dirname "$0")"
VINA="/c/Program Files (x86)/The Scripps Research Institute/Vina/vina.exe"
PY="/c/Users/aritr/.conda/envs/docking_project/python.exe"
declare -A DIR=( [2V5Z]=mao [2Z5X]=mao [7JXX]=ttbk [7Q8Y]=ttbk )
for R in 2V5Z 2Z5X 7JXX 7Q8Y; do
  B=../../03_receptors/${DIR[$R]}/$R/box.json
  # Windows python prints CRLF; the trailing CR would poison the last argument, so strip it.
  read -r CX CY CZ SX SY SZ < <("$PY" -c "import json;b=json.load(open('$B'));print(*b['center'],*b['size'])" | tr -d '\r')
  for S in 11 22 33; do
    "$VINA" --receptor ${R}_meeko.pdbqt --ligand cand_003.pdbqt \
      --center_x $CX --center_y $CY --center_z $CZ --size_x $SX --size_y $SY --size_z $SZ \
      --scoring vina --exhaustiveness 32 --num_modes 9 --seed $S \
      --out ${R}_seed${S}_out.pdbqt > ${R}_seed${S}.log 2>&1
    echo "$R seed$S: $(grep -m1 -E '^ +1 ' ${R}_seed${S}.log)"
  done
done
echo ALLDONE
