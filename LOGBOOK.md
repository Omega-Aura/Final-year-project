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
