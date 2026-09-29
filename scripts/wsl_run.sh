#!/usr/bin/env bash
# Run an RDKit/meeko-dependent project command inside WSL2.
#
# WHY THIS EXISTS
# On 2026-09-28 Windows Smart App Control (VerifiedAndReputablePolicyState = 1, usermode Code
# Integrity enforced) began blocking two unsigned conda-forge binaries in the native
# docking_project env:
#
#     Lib\site-packages\rdkit\Chem\rdchem.pyd   -> "An Application Control policy has blocked this file"
#     Scripts\mk_export.exe                     -> "blocked by your organization's Device Guard policy"
#
# Smart App Control has NO exclusion list, so no Defender path exclusion can reach them:
# Add-MpPreference configures the antivirus, a different subsystem. Confirmed twice -- exclusions
# were added and both files stayed blocked. The only alternatives are turning Smart App Control
# off, which cannot be undone without reinstalling Windows, or running the affected code somewhere
# the policy does not apply. WSL2 runs ELF binaries and is not governed by it, which is also why
# AmberTools in the `mdgbsa` env kept working throughout.
#
# rdchem.pyd is only the first of 63 .pyd files in that tree, so unblocking it individually would
# likely just move the error to the next one. Moving the whole RDKit-dependent step is the fix that
# holds.
#
# WHAT IT AFFECTS
# Blocked:     prep_ligands.py, prep_receptor.sh, filter_cascade.py, bbb_score.py,
#              fix_native_bondorders.py, and PDBQT -> SDF export.
#
#              NOT run_prolif_v2.py: that one is prior-phase and cannot run at all, because its
#              inputs (prep/, ligands/, docking/) were removed in the cleanup. ProLIF and
#              MDAnalysis are installed in `ligprep` anyway, so it has a home if ever revived.
# Unaffected:  Vina scoring, OpenMM MD (CUDA), WSL MM-GBSA, every analysis and figure script.
#
# DESIGN
# The project's scripts take relative paths and are run from the project root, so the whole script
# can be executed inside WSL unchanged -- no per-binary shim, no argument rewriting, nothing for a
# caller to get subtly wrong. `shutil.which` inside those scripts then finds the Linux entry point
# names (mk_export.py, not mk_export.exe) on its own.
#
# The `ligprep` env is pinned to the SAME rdkit (2025.09.6) and meeko (0.7.1) as the Windows env,
# so files produced here are comparable with the ones already in the repo. It is a drop-in
# replacement for the blocked binaries, not a second toolchain. It is kept separate from `mdgbsa`
# because a dependency solve that broke that AmberTools install would cost far more than it saves.
#
# Usage:
#   bash scripts/wsl_run.sh python scripts/prep_ligands.py --csv 01_smiles/candidates_56.csv -o 02_ligands
#   bash scripts/wsl_run.sh mk_export.py in.pdbqt -s out.sdf
#   bash scripts/wsl_run.sh --check          # verify the env can do real work, then exit
set -uo pipefail

ENV_NAME="${LIGPREP_ENV:-ligprep}"
PROJ_WSL="/mnt/c/Users/aritr/OneDrive/Desktop/Final year project"

# Run a command inside the ligprep env, from the project root, and resolve meeko's entry points to
# their Linux names. The entry point is mk_export.py / mk_prepare_ligand.py on conda-forge Linux
# and mk_export.exe / mk_prepare_ligand.exe on Windows -- same meeko version, different script
# name -- so a caller that hardcodes either one breaks on the other platform. Aliases let the
# scripts and the command line use whichever name they already use.
in_wsl() {
    wsl -e bash -lc "
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh || exit 1
conda activate '$ENV_NAME' || { echo '[abort] no conda env $ENV_NAME in WSL' >&2; exit 1; }
cd '$PROJ_WSL' || { echo '[abort] cannot cd to project root in WSL' >&2; exit 1; }
for t in mk_export mk_prepare_ligand mk_prepare_receptor; do
    command -v \$t >/dev/null 2>&1 || eval \"\$t() { \$t.py \\\"\\\$@\\\"; }\"
done
$1"
}

# --check does a REAL operation, never a --help probe. A --help probe already produced one false
# "OK" for mk_export.exe this session: bash printed "Permission denied" and the probe read the
# exit status as success. So the check embeds a molecule, optimises it, writes a PDBQT and reads it
# back -- the actual work the prep scripts do.
if [ "${1:-}" = "--check" ]; then
    echo "checking WSL env '$ENV_NAME' ..."
    in_wsl '
python - <<'"'"'EOF'"'"'
import subprocess, shutil, sys, tempfile, os
import rdkit, meeko
from rdkit import Chem
from rdkit.Chem import AllChem
print("  rdkit  ", rdkit.__version__, "(Windows env has 2025.09.6)")
print("  meeko  ", meeko.__version__, "(Windows env has 0.7.1)")
m = Chem.AddHs(Chem.MolFromSmiles("O=c1c(O)c(-c2ccccc2)oc2ccccc12"))   # a flavonol, like the real set
ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
assert AllChem.EmbedMultipleConfs(m, numConfs=3, params=ps), "embed failed"
AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
d = tempfile.mkdtemp()
sdf = os.path.join(d, "t.sdf"); pq = os.path.join(d, "t.pdbqt")
Chem.SDWriter(sdf).write(m, confId=0)
for tool, args in (("obabel", [sdf, "-O", sdf, "-p", "7.4"]),
                   ("mk_prepare_ligand.py", ["-i", sdf, "-o", pq])):
    exe = shutil.which(tool) or sys.exit("  MISSING: %s" % tool)
    r = subprocess.run([exe] + args, capture_output=True)
    if r.returncode != 0:
        sys.exit("  FAILED: %s -> %s" % (tool, r.stderr.decode()[:300]))
n = sum(1 for L in open(pq) if L.startswith(("ATOM", "HETATM")))
assert n > 10, "pdbqt has only %d atoms" % n
print("  embed + MMFF + obabel + mk_prepare_ligand: OK (%d atoms written)" % n)
EOF'
    rc=$?
    [ $rc -eq 0 ] && echo "[ok] '$ENV_NAME' can do real RDKit/meeko work" \
                  || echo "[FAIL] '$ENV_NAME' check failed (exit $rc)" >&2
    exit $rc
fi

[ $# -ge 1 ] || { echo "usage: wsl_run.sh <command...>   |   wsl_run.sh --check" >&2; exit 2; }

# Quote each argument for the remote shell so paths containing spaces survive. The project root
# itself sits under "Final year project", so this is not hypothetical.
CMD=""
for arg in "$@"; do
    CMD="$CMD '$(printf '%s' "$arg" | sed "s/'/'\\\\''/g")'"
done

in_wsl "$CMD"
