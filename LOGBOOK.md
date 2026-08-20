# Logbook

## 2026-08-20 — A

Day 1 environment setup (§Part 2, Day 1).

**Native Windows, conda env `docking_project`:** was rdkit-only at start. Installed
`openmm mdanalysis prolif parmed` via conda-forge, and `spyrmsd meeko dimorphite-dl`
via pip. `python -m openmm.testInstallation` confirms a working CUDA platform
(RTX 4050, driver CUDA UMD 13.3) — MD runs (§A6–A8) can run natively, no WSL needed
for the simulation itself.

**AmberTools has no Windows build** (confirmed: zero conda-forge win-64 matches for
`ambertools`). This blocks `acpype` and `MMPBSA.py` natively. Machine already had
WSL2 Ubuntu with GPU passthrough working (`nvidia-smi` visible inside WSL). Installed
Miniconda in WSL2 and created a minimal env `mdgbsa` (`ambertools acpype`, conda-forge
only — had to use `--override-channels` to skip the `defaults` channel's ToS gate).
Both `acpype` (v2023.10.27) and `MMPBSA.py` (v14.0) confirmed callable inside WSL.

**Net environment split for Machine A:**
- Native Windows `docking_project`: rdkit, openmm(CUDA), MDAnalysis, prolif, parmed,
  spyrmsd, meeko, dimorphite-dl, vina, obabel — docking + MD + analysis.
- WSL2 `mdgbsa`: acpype (GAFF2 parameterization), MMPBSA.py (MM-GBSA) only — used
  because AmberTools has no Windows port, not a general duplicate stack.

**NOTE:** meeko's CLI installs as `mk_prepare_receptor.exe` / `mk_prepare_ligand.exe`
on Windows, not `mk_prepare_receptor.py` as written in the doc's Part 6 scripts. Both
`mk_prepare_receptor` (no extension) and the `.exe` form work from bash. Scripts in
`scripts/` need this accounted for per machine, per the doc's own warning about
Meeko's CLI changing across versions.

**STILL OPEN:** the lead compound's acpype/GAFF2 parameter files (named as one of the
three Day-1-critical files) were not found anywhere in the recovered `artifacts/` —
no `.itp`/`.top`/`.frcmod`/`.mol2`. Will need to be regenerated via `acpype` (WSL)
before Week 2's MD system builds (§A6), which assume these already exist and are
just reused across all four systems.

Files: none yet (environment-only session, no data outputs).

## 2026-08-20 — B

- Initialized environment and scripts in `scripts/` (`prep_receptor.sh`, `prep_ligands.py`, `dock.sh`, `rmsd_check.py`, `collect_results.py`).
- Prepared MAO-B receptor (2V5Z Chain A) with FAD cofactor using `scripts/prep_receptor.sh`.
- Extracted and prepared native safinamide ligand (SAG) with `scripts/prep_ligands.py`.
- Executed calibration redocking with `scripts/dock.sh 2V5Z native_SAG 11` (exhaustiveness=32, seed=11).
- Computed symmetry-corrected RMSD with `scripts/rmsd_check.py`:
  - Pose 1 RMSD: 1.38 Å
  - Pose 3 (best) RMSD: 0.56 Å -> **PASS** (< 2.0 Å threshold).
- **GATE 0 (Calibration exercise) passed** on B's machine.

NOTE (A, reviewing): §Day 2 requires *both* machines to reproduce the ~1.57 Å
safinamide redock and agree with each other to within ~0.3 Å before proceeding.
B's best pose (0.56 Å) is a clean pass in isolation but hasn't yet been
cross-checked against an independent run on Machine A — that comparison is still
outstanding before GATE 0 can be called fully closed per the doc's own rule.

## 2026-08-20 — A

Day 3: structure download and inspection (§Part 2, Day 3). A's task in full.

**Downloaded all six receptors** into `03_receptors/<PDB>/raw.pdb` (2V5Z already
present from B's Day 2 work, skipped re-download). NOTE: the doc's download loop
literally lists eight IDs (`7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X 7JXY 4BTM`) under a "six
receptors" heading — `7JXY` and `4BTM` appear nowhere else in the doc (not in the
table, not in any later step) and look like a copy/typo artifact (near-duplicates of
7JXX/4BTK). Went with the doc's own 6-receptor table as authoritative and did not
fetch those two. Flag if that's wrong.

**Inspection table** (HETATM species, missing-residue count from REMARK 465,
resolution — all pulled directly from the downloaded headers):

| PDB | Protein | Chains | Res. | HETATM species | REMARK 465 lines | Native ligand |
|---|---|---|---|---|---|---|
| 7JXX | TTBK1 | A | 1.56 Å | NA, HOH, VP7 | 51 | VP7 |
| 4BTK | TTBK1 | A | 2.00 Å | DMS, DTQ, HOH | 57 | DTQ |
| 7Q8V | TTBK1 | A | 2.13 Å | 9IV, HOH, PO4 | 24 | 9IV |
| 7Q8Y | TTBK2 | A, B | 1.60 Å | 9IV, HOH, PO4 | 20 | 9IV (present in both chains; using chain A for consistency with 2V5Z convention) |
| 2V5Z | MAO-B | A, B | 1.60 Å | FAD, HOH, SAG | 54 | SAG |
| 2Z5X | MAO-A | A | 2.20 Å | DCX, FAD, GOL, HOH, HRM | 0 | HRM |

All resolutions and native ligands match the doc's table exactly — no surprises there.

**Check 1 — does 2Z5X contain FAD?** Yes, confirmed: `HET FAD A 600 53` (full
53-atom FAD, not partial), plus a REMARK 500 close-contact record between Cys406 SG
and the FAD C8M atom (1.65 Å) — consistent with the expected 8α-S-cysteinyl
covalent linkage at Cys406. **The custom FAD-transplant step described in the
manuscript's Methods was unnecessary** — 2Z5X ships with a complete, correctly
linked FAD. This deletes a caveat from Limitations.

**Check 2 — do 7JXX pocket residues match 4NFM numbering?** Superposed 7JXX chain A
onto `artifacts/4NFM.pdb` chain A in PyMOL (`cmd.align`, 0.43 Å RMSD over 1815 atoms
before refinement). Checked all seven pocket residues named in the manuscript
(ILE40, ILE48, GLN89, GLN110, GLY111, ASN113, LEU175) directly against both raw PDB
files (not just the alignment tool, since the alignment's raw-pair extraction gave a
false "not aligned" for residue 40 that direct inspection contradicted):

| 4NFM | 7JXX | Match |
|---|---|---|
| ILE40 | ILE40 | identical number + identity |
| ILE48 | ILE48 | identical number + identity |
| **GLN89** | HIS89 (4NFM is *also* HIS89, not GLN) | **manuscript residue name is wrong in both structures** |
| GLN110 | GLN110 | identical number + identity |
| GLY111 | GLY111 | identical number + identity |
| ASN113 | ASN113 | identical number + identity |
| LEU175 | LEU175 | identical number + identity |

**7JXX uses the exact same numbering as 4NFM** — no residue-numbering lookup table
is needed for §2.2. The one required fix is a residue-identity correction:
**GLN89 → HIS89** in the manuscript's pocket-residue list; it was never a numbering
issue, and 4NFM (the structure the manuscript numbers are drawn from) has HIS at
that position too, so this looks like a plain transcription error, not something
that changed between structures.

**GATE 0 status:** still open pending A's own independent calibration run to cross-
check against B's 0.56 Å (see note above) — did not run it today, scope was Day 3
only per instruction.

Files: `03_receptors/{7JXX,4BTK,7Q8V,7Q8Y,2Z5X}/raw.pdb` (new), `03_receptors/2V5Z/raw.pdb` (pre-existing, B).

## 2026-08-20 — A

Track A, Week 1 (§Part 3). A's tasks A1-A5.

**Fixed two portability bugs in shared scripts before running anything**, both
harmless-on-this-machine but violating §1.5 (same script, same run):
`scripts/prep_receptor.sh` and `scripts/dock.sh` both hardcoded B's personal path
(`/c/Users/sayan/miniforge3/envs/dock/...`) into `$PATH`. Removed; both scripts
already autodetect the active conda env correctly without it.

**A1 — prepare all six receptors.** Two (7JXX, 7Q8Y) initially failed: alternate-
location (altloc) residues that `-a/--allow_bad_res` alone doesn't resolve. Added
`--default_altloc A` to the `mk_prepare_receptor` call in `prep_receptor.sh` and
re-ran **all six** uniformly (not just the two failures) per §1.5. All six now
prepare cleanly: 7JXX, 4BTK, 7Q8V, 7Q8Y, 2V5Z, 2Z5X.

**A2 — GATE 1 redocking validation (7JXX, 4BTK).** Found a real pipeline bug while
running this: 7JXX's redock crashed spyrmsd with `NonIsomorphicGraphs` instead of
producing a number. Root cause: `dock.sh` used plain `obabel` to convert Vina's
docked PDBQT back to SDF. Meeko inserts dummy "glue" atoms into the PDBQT for
ligands needing flexible-ring handling (VP7 apparently needs this, DTQ apparently
doesn't -- explains why 4BTK "passed" and 7JXX crashed instead of just scoring
badly). `obabel` doesn't understand Meeko's own convention and silently emits `*`
wildcard atoms, corrupting the molecular graph. Fix: use Meeko's own `mk_export`
to convert docked poses, not `obabel` (verified: `mk_export` reconstructs VP7 as
45 atoms / 48 bonds, exact match to the reference; `obabel` gave 28 atoms with two
`*` atoms). Re-exported existing poses with `mk_export` (no need to re-dock,
deterministic given the same seed) and re-ran the RMSD check:

  - **7JXX (native VP7): best pose 0.71 Å -> PASS.**
  - **4BTK (native DTQ): best pose 0.75 Å -> PASS** (unchanged by the fix; DTQ
    didn't hit the dummy-atom issue, confirms the fix doesn't regress a working case).

**GATE 1 PASSES on 7JXX**, the doc's designated primary TTBK1 receptor (1.56 Å,
holo, VP7-bound) -- a dramatic improvement over the manuscript's original 5.29 Å
failure on the apo 4NFM structure. Per §Part 3: **Limitations #2 and #3 are deleted
from the manuscript.** No fallback to 4BTK needed, though 4BTK remains available as
a validated cross-check receptor.

This `mk_export` fix matters beyond GATE 1 -- `dock.sh` is used for every docking
run in the project, including A5's 56-candidate consensus docking below. Fixed
before running A5, not after.

**A3 — rebuild the 56-candidate ligand set from SMILES.** `01_smiles/candidates_56.csv`
(the `cand_001`...`cand_056` reformatted copy of `rl_v2_shortlist_56.csv`) rebuilt
from scratch via `prep_ligands.py`. 56/56 succeeded, 0 failures.

**Found and fixed a third `dock.sh` bug while testing the multi-ligand path** (never
exercised before today -- B's GATE 0 calibration and my GATE 1 validation were both
single-ligand runs): the set-resolution logic fell back to globbing *every*
`.pdbqt` file in `02_ligands/pdbqt/` when there's no single combined
`<SET>.pdbqt` file, which would have silently swept the native reference ligands
(native_VP7, native_DTQ, native_SAG, native_9IV) into any "candidates_56" docking
run. Fixed: when `01_smiles/<SET>.csv` exists, the ligand list now comes from that
CSV's `id` column (the authoritative set membership), not a blind glob.

**Then hit a fourth bug from that same fix**: the python-generated ligand-path list
inside `dock.sh` picked up trailing `\r` on every line (Python's stdout on Windows
does universal-newline translation to `\r\n`; bash's `$()` doesn't strip embedded
`\r`). This corrupted every path, so `basename "$L" .pdbqt` didn't strip the
suffix, vina couldn't find the file, and `set -euo pipefail` killed each of the 6
candidates_56 docking runs after only the first ligand. First A5 attempt silently
produced only 6 log files total (1 per seed/receptor combo) instead of 336.
Caught this by checking `collect_results.py`'s output table (only 3 rows -- the
native ligands -- with zero candidate entries) rather than trusting a clean exit
code. Fixed with `| tr -d '\r'` on the python output; verified the fix directly
(56/56 correct paths, no `\r`, file existence confirmed) before re-running the full
batch.

**A4 — re-run the filtering cascade.** Hard blocker, not a version-drift issue: the
actual filter code was never recovered at all, only its output CSVs. `bbb_score.py`
is a different thing (REINVENT4's RL-time BBB reward hook, not the offline
Lipinski/GI/BBB/alerts cascade), and it itself imports a `filtering/boiled_egg_coords.py`
that also didn't exist anywhere in the project.

Asked how to proceed; reconstructed from published/standard definitions rather than
holding indefinitely:
  - **Lipinski Ro5**: standard, pass = <=1 violation of {MW>500, WLogP>5, HBD>5, HBA>10}.
  - **GI absorption / BBB penetration**: point-in-polygon against the digitized
    BOILED-Egg ellipses from Daina & Zoete (ChemMedChem 2016), coordinates pulled
    from the open-source reimplementation PyBOILEDegg (github.com/bfmilne/PyBOILEDegg,
    GPLv3), which cites the same source paper. Saved to the (recreated)
    `filtering/boiled_egg_coords.py` that `bbb_score.py` already expected to exist.
    TPSA must use `includeSandP=True` -- confirmed against the existing TPSA column
    in the recovered data (0/50 mismatch with S/P included vs 1/50 without).
  - **Structural alerts**: RDKit's combined PAINS + BRENK `FilterCatalog` -- alert
    names in the recovered data (`Michael_acceptor_1`, `catechol_A(92)`, etc.)
    matched this combination's naming convention.

**Verified against `00_library/reinvent4_output/campaign2_v2/rl_v2_filtered_full.csv`
(405 molecules, already carrying the original pass/fail columns) -- zero mismatches
on all four filters, individually, per molecule** (not just matching aggregate
counts, which could be coincidental -- every single one of 405 x 4 boolean labels
matched). Funnel: 369 Lipinski / 81 BBB / 334 GI / 253 alert-free / 56 passing all
-- **exact match to the doc's target numbers.**

Wrote this up as a real, re-runnable script (`scripts/filter_cascade.py`), not just
inline verification code, since that's what "re-run the filtering cascade" actually
requires. Ran it independently on the 405-molecule pool (369/81/334/253/56, exact)
and on our own `01_smiles/candidates_56.csv` (56/56 pass all four -- confirms A3's
ligand set is internally consistent with A4's filter).

**Side finding, not fixed (out of scope -- the RL run itself isn't being redone):**
`bbb_score.py`'s TPSA call is missing `includeSandP=True`, so the RL-time BBB reward
and the actual offline BBB filter are not perfectly consistent for S/P-containing
molecules, contradicting that script's own docstring claim. Doesn't affect anything
in this project's remaining phases since the generative campaigns are being kept
as-is, not rerun.

**A5 — consensus re-dock all 56 candidates on 7JXX and 2V5Z.** 3 seeds x 56
ligands x 2 receptors = 336 dockings. First attempt silently broke on the `\r` bug
above (caught before trusting it). Re-ran with the fixed `dock.sh` -- **complete,
336/336 logs, `collect_results.py` -> `08_analysis/consensus_new.csv`.**

**Ranking on the validated 7JXX changed, as the doc predicted.** The original
pipeline's top-ranked dual-target candidate by RL score, `cand_001`
(`rl_v2_shortlist_56.csv` row 1, `Score`=0.7497, the highest in the set), is now
**#3/56 on 7JXX** (consensus -8.79 kcal/mol) -- `cand_013` (-8.83) and `cand_002`
(-8.80) both edge it out on the holo, validated receptor. Consistent with the doc's
"a holo pocket at 1.56 Å is a different shape from an apo pocket" expectation.

**More striking: `cand_001` ranks only #49/56 on 2V5Z (MAO-B)**, consensus -10.32
kcal/mol, well off the top cluster (best on 2V5Z: cand_043 at -11.63). 2V5Z itself
wasn't revalidated today (it's the doc's "unchanged" primary MAO-B receptor, already
established), so this isn't a validation artifact -- it's either real chemistry
(cand_001 may simply bind TTBK1 much better than MAO-B despite being scored as a
top dual-target hit during RL) or an artifact of how the RL reward combined
`TTBK1_dock` and `MAOB_dock` into one scalar score (a compound can rank #1 on a
sum/product of two scores while being mediocre on one of them). Worth raising at
sync -- this bears directly on which compound A should advance as "the lead" for
MM-GBSA in Week 2-3, and on the dual-engagement claim generally.

New top TTBK1 hit: **cand_013** (-8.83 kcal/mol, seed SD 0.005 -- very stable
across seeds). New top MAO-B hit: **cand_043** (-11.63 kcal/mol).

Files: `03_receptors/*` (all six, receptor.pdbqt/box.json), `02_ligands/{sdf,pdbqt}/cand_*`
(56 each), `filtering/boiled_egg_coords.py` (new), `scripts/filter_cascade.py` (new),
`08_analysis/filter_cascade_A4_check.csv`, `08_analysis/filter_cascade_candidates_56.csv`,
`08_analysis/consensus_new.csv`, `05_validation/7JXX_redock.txt`, `05_validation/4BTK_redock.txt`,
`04_docking/{7JXX,2V5Z}_candidates_56_seed{11,22,33}/` (336 logs total).

**Track A Week 1 (A1-A5) is now complete.** Per §Part 3 checklist: all six receptors
prepared by one script [x]; GATE 1 (7JXX < 2.0 Å) [x]; 56 candidates rebuilt [x];
filtering cascade re-run, numbers reproduce exactly [x]; 56 candidates re-docked on
7JXX + 2V5Z [x]. Waiting on B's Track B Week 1 (reference set + benchmark docking)
before SYNC POINT 1.
