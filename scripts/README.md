# Scripts — Pipeline Drivers

Every step of the workflow is re-runnable from here. Scripts are grouped by the step they drive.

## Preparation

| Script | Drives | Notes |
|---|---|---|
| [`prep_receptor.sh`](prep_receptor.sh) | [step 3](../03_receptors/) | receptor → `receptor.pdbqt` + `box.json` |
| [`prep_ligands.py`](prep_ligands.py) | [step 2](../02_ligands/) | SMILES → SDF + PDBQT; `--ph 7.4 --nconf 20` |
| [`prep_cofactor.py`](prep_cofactor.py) | [step 3](../03_receptors/), [step 6](../06_md/) | FAD from RCSB ideal chemistry → GAFF2 residue |
| [`fix_native_bondorders.py`](fix_native_bondorders.py) | [step 5](../05_validation/) | repairs PDB-derived bond orders before RMSD |
| [`pdbqt_utils.py`](pdbqt_utils.py) | shared | PDBQT parsing helpers |

## Docking

| Script | Drives | Notes |
|---|---|---|
| [`dock.sh`](dock.sh) | [step 4](../04_docking/) | the Vina call; `--exhaustiveness 32 --num_modes 9 --seed N` |
| [`docking_paths.sh`](docking_paths.sh) | [step 4](../04_docking/) | sourced, not run: the **only** definition of which `04_docking/` group a ligand set writes into |
| [`receptor_paths.sh`](receptor_paths.sh) | [step 3](../03_receptors/) | sourced, not run: resolves a receptor by PDB id under `03_receptors/<family>/`, and holds the **only** family rule |
| [`project_paths.py`](project_paths.py) | [step 3](../03_receptors/) | the same receptor resolution for Python readers — glob only, **no** family rule, so nothing can drift against the shell copy |
| [`md_paths.sh`](md_paths.sh) | [step 6](../06_md/) | sourced, not run: resolves an MD system name under `06_md/systems/` |
| [`run_week2_redock.sh`](run_week2_redock.sh) | [step 4](../04_docking/) | production re-dock, **resumable** |
| [`run_water_test.sh`](run_water_test.sh) | [step 4](../04_docking/) | paired wet/dry docking |
| [`run_minwater_probe.sh`](run_minwater_probe.sh) | [step 4](../04_docking/) | minimal-water probe |
| [`run_selectivity.sh`](run_selectivity.sh) | [step 4](../04_docking/) | anti-target counter-screen |
| [`collect_results.py`](collect_results.py) | [step 8](../08_analysis/) | Vina logs → consensus score table |

## Validation and analysis

| Script | Drives | Notes |
|---|---|---|
| [`rmsd_check.py`](rmsd_check.py) | [step 5](../05_validation/) | redocking RMSD, symmetry-corrected via `spyrmsd` |
| [`calibration_9iv.py`](calibration_9iv.py) | [step 5](../05_validation/) | 9IV margin vs experiment → systematic bias |
| [`filter_cascade.py`](filter_cascade.py) | [step 8](../08_analysis/) | Lipinski / GI / BBB / alerts cascade |
| [`bbb_score.py`](bbb_score.py) | [step 8](../08_analysis/) | BOILED-Egg scoring |
| [`analyze_selectivity.py`](analyze_selectivity.py) | [step 8](../08_analysis/) | per-candidate selectivity margins |
| [`analyze_water_test.py`](analyze_water_test.py) | [step 8](../08_analysis/) | wet-minus-dry deltas |
| [`run_prolif_v2.py`](run_prolif_v2.py) | [step 4](../04_docking/) | protein–ligand interaction fingerprints |

## MD and free energy

| Script | Drives | Notes |
|---|---|---|
| [`../06_md/run_md_system.py`](../06_md/) | [step 6](../06_md/) | builds a solvated system |
| [`../06_md/run_md.py`](../06_md/) | [step 6](../06_md/) | minimise → equilibrate → produce |
| [`../06_md/run_md_extend.py`](../06_md/) | [step 6](../06_md/) | continues from `final_state.xml` |
| [`../06_md/run_md_restrained.py`](../06_md/) | [step 6](../06_md/) | FAD-restrained variant for MAO systems |
| [`run_ttbk2_replicate.sh`](run_ttbk2_replicate.sh) | [step 6](../06_md/) | velocity replicate of `system_TTBK2_p2`, the only on-pose TTBK2 trajectory; **aborts** unless its inputs are byte-identical to run 1 |
| [`run_vinardo_crosscheck.sh`](run_vinardo_crosscheck.sh) | [step 4](../04_docking/) | independent scoring-function cross-check on the **validated** receptors; writes to `04_docking/crosscheck/vinardo/`, deeper than every vina glob |
| [`analyze_vinardo_crosscheck.py`](analyze_vinardo_crosscheck.py) | [step 8](../08_analysis/) | Vina vs Vinardo correlation and rank agreement; never pools or averages the two functions |
| [`run_mmgbsa.sh`](run_mmgbsa.sh) | [step 7](../07_mmgbsa/) | strip → `ante-MMPBSA.py` → `MMPBSA.py`, inside WSL2 |

## Three environments, one pipeline

AmberTools has no native Windows build, so parameterisation and MM-GBSA run in a **WSL2** conda
environment (`mdgbsa`) while OpenMM runs **natively on the GPU** in the Windows `docking_project`
environment. Scripts that cross the boundary call `wsl -e bash -lc` and translate paths to
`/mnt/c/...`; [`run_mmgbsa.sh`](run_mmgbsa.sh) is the reference example.

Since 2026-09-28 there is a third: **`ligprep`**, also in WSL2, for everything that needs RDKit.

### Why RDKit moved out of Windows

Windows Smart App Control blocks two unsigned conda-forge binaries in the native env:
`rdkit\Chem\rdchem.pyd` and `Scripts\mk_export.exe`. Smart App Control has **no allowlist**, so a
Defender `-ExclusionPath` cannot reach them -- that configures the antivirus, a different subsystem.
Disabling Smart App Control is irreversible without reinstalling Windows. WSL2 runs ELF binaries and
is not governed by it, which is why `mdgbsa` was unaffected throughout. See LOGBOOK entry M.

**Run these through the bridge**, because they import RDKit:

```bash
bash scripts/wsl_run.sh python scripts/prep_ligands.py --csv 01_smiles/candidates_56.csv -o 02_ligands
bash scripts/wsl_run.sh python scripts/filter_cascade.py ...
bash scripts/wsl_run.sh python scripts/bbb_score.py ...
bash scripts/wsl_run.sh python scripts/fix_native_bondorders.py ...
bash scripts/wsl_run.sh --check      # verify the env by doing real work, not a --help probe
```

`prep_receptor.sh` is a shell script that calls `obabel` directly rather than importing RDKit, so
route it the same way if it is ever re-run.

`run_prolif_v2.py` is **not** in this list. It is prior-phase and cannot run at all -- its inputs
(`prep/`, `ligands/`, `docking/`) were removed in the cleanup, so it would fail with or without
RDKit. ProLIF 2.2.1 and MDAnalysis 2.9.0 are installed in `ligprep` regardless, so it has a home if
it is ever brought forward to the numbered tree.

Unaffected, still native: Vina scoring, OpenMM MD on CUDA, and every analysis and figure script.

`ligprep` is pinned to the **same rdkit (2025.09.6) and meeko (0.7.1)** as the Windows env so it is
a drop-in replacement, not a divergent toolchain. `openbabel` could not be matched (3.2.1 vs 3.1.1;
3.1.1 conflicts with python 3.11 + rdkit 2025.09.6 via libxml2). That was verified rather than
assumed harmless: WSL-prepped conformers are internally identical to the repo's (`GetBestRMS` =
0.000, the difference being rigid-body frame only), and re-docking one at the identical box and seed
reproduced the repo scores to **<=0.02 kcal/mol**, well under the inter-seed SD.

Env spec: [`env_ligprep_wsl.yml`](env_ligprep_wsl.yml).

**Existing ligands do not need re-prepping.** `ligprep` is for new work. Re-prepping is chemically
safe but perturbs the absolute frame for no gain.

## Path dependencies — read before reorganising

These scripts hardcode the numbered step directories in **16+ places**: `03_receptors/`,
`01_smiles/`, `04_docking/`, `02_ligands/pdbqt/`, `06_md/systems/`, `05_validation/`,
`08_analysis/`, `00_library/reinvent4_output/campaign2_v2/`.

**Renaming or moving a numbered directory breaks reproducibility.** Add to the structure rather than
rearranging it. This is why the workspace keeps its `00_`–`09_` layout and documents each step in
place with a README instead of regrouping files into new folders.

## Traps these scripts already handle

Worth knowing about before writing new ones:

- **Silent docking failure.** Vina reports a preparation failure as `BEST None 1000000000.00` — a
  failure that reads like a score. `dock.sh` catches it.
- **Truncated SDF.** Detected by checking for the record terminator, which is the last thing
  written, so a partially written file is distinguishable from a complete one.
- **No blind globbing of `02_ligands/pdbqt/`.** It is now split into `candidates/`, `natives/` and
  `references/`, and `dock.sh` resolves a ligand by name at either depth via its `lig_path` helper,
  so a half-migrated tree cannot silently resolve to nothing. The directory still holds more ligands than any single run
  should use; ligand lists are explicit.
- **`MMPBSA.py -sp` is wrong for a stripped trajectory.** `-sp` expects the solvated topology and
  aborts on an atom-count mismatch. `run_mmgbsa.sh` omits it deliberately, with the reasoning in a
  comment.
- **Bond orders from PDB are unreliable**, and a wrong bond order changes the symmetry graph an RMSD
  is computed over — hence `fix_native_bondorders.py`.
