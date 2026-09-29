# Workspace Inventory and Classification

Pass 1 of the cleanup: every item inventoried, classified, dependency-checked, and either kept in
place, relocated into its workflow step, or quarantined. **Nothing has been permanently deleted.**

Generated 26 Sep 2026.

## Classification summary

| Category | Meaning | Outcome |
|---|---|---|
| **FINAL** | directly used in the final analysis / proposal | kept in place |
| **SUPPORTING** | needed to explain, reproduce or validate a final result | kept in place |
| **REFERENCE** | papers, protocols, parameter sources, characterisation docs | kept in place |
| **INTERMEDIATE** | temporary byproducts, regenerable | quarantined (1.83 GB) |
| **OUTDATED** | superseded and carrying no remaining scientific meaning | quarantined (190 MB) |
| **DUPLICATE** | byte-identical copy of a file that stays live | quarantined (198 MB) |
| **REVIEW_REQUIRED** | could not be classified from content | resolved 27 Sep: the two project plans moved to `provenance/` |

## Where everything stands

| Directory | Size | Files | Classification |
|---|---|---|---|
| `00_library/` | 1.5 MB | 26 | REFERENCE + SUPPORTING — curated library with PMIDs; RL campaign history |
| `01_smiles/` | 56 KB | 8 | **FINAL** — the frozen candidate and reference sets |
| `02_ligands/` | 792 KB | 165 | SUPPORTING — prepared 3D inputs, regenerable but cheap to keep |
| `03_receptors/` | 15 MB | 91 | **FINAL** + REFERENCE — docking targets, controls, characterisation |
| `04_docking/` | 101 MB | 5,160 | **FINAL** — all production docking, incl. wet/dry and cofactor controls |
| `05_validation/` | 126 KB | 11 | **FINAL** — redocking and calibration; the protocol's credibility rests here |
| `06_md/` | 18 GB | 675 | **FINAL** + SUPPORTING — 17 MD systems, raw + stripped trajectories |
| `10_results/` | <1 MB | 5 | **FINAL** — consolidated results, interpretation, 3 figures + their generator |
| `07_mmgbsa/` | 296 KB | 7 | **FINAL** — collected energy summary + prior-phase results |
| `08_analysis/` | 700 KB | 33 | **FINAL** — consensus, cascade, selectivity, water test |
| `09_manuscript/` | 5.6 MB | 6 | **FINAL** — the deliverable |
| `filtering/` | 20 KB | 3 | SUPPORTING — the BBB model, used by both cascade and RL reward |
| `generation/` | 95 MB | 12 | **FINAL** + SUPPORTING — v3 campaign, config and checkpoint |
| `scripts/` | 116 KB | 20 | **FINAL** — every step is re-runnable from here |
| `REINVENT4/` | 1.1 GB | 3,679 | third-party clone, untracked; required to re-run generation |
| `provenance/` | 47 KB | 5 | original project plans + the record of what the cleanup removed |

Workspace total: **11 GB**, of which 2.1 GB is quarantined.

## Dependency check (rule 5) — what it changed

Before moving anything, every candidate was checked for being referenced by a script, needed to
reproduce a result, or the only copy of something. Two findings redirected the whole approach:

### The numbered directories are load-bearing — so they were not renamed

`scripts/` and `generation/` hardcode the step directories in **16+ places**: `03_receptors/`,
`01_smiles/`, `04_docking/`, `02_ligands/pdbqt/`, `06_md/systems/`, `05_validation/`,
`08_analysis/`, `00_library/reinvent4_output/campaign2_v2/`.

Renaming or regrouping them would have broken reproducibility for no scientific gain. The
restructure therefore **adds a documentation layer in place** rather than rearranging files: a
`README.md` per step answering *what / why / reference / inputs / analysis / result / files*, plus
[`README.md`](README.md) at the root as the workflow index.

### `artifacts/` was not a duplicate mirror — it was the missing documentation

The initial read was that `artifacts/` (202 MB) duplicated other directories. Hash comparison showed
**144 of its 174 files were unique**, including the entire `phase0`–`phase8` report series that
[`LOGBOOK.md`](LOGBOOK.md) cites constantly, plus the only copies of the earlier receptor set and the
prior-phase MM-GBSA results.

Deleting it on the filename-level assumption would have destroyed the project's rationale trail. Its
30 verified duplicates were quarantined; its 144 unique files were distributed into the workflow step
each belongs to, under `<step>/prior_phase/`. `artifacts/` no longer exists.

`prior_phase/` is kept as a subfolder rather than merged into each step's root because the two passes
used **different receptor structures** (`4NFM`/`6U0K`/`2V60` vs `7JXX`/`7Q8Y`/`2V5Z`/`2Z5X`). Pooling
them would invite exactly the kind of silent mismatch this project was bitten by twice.

## Quarantined — and why each is safe

| Item | Size | Category | Evidence |
|---|---|---|---|
| `06_md/*/reference.frc` (7) | 1.83 GB | INTERMEDIATE | `MMPBSA.py` forces dump. Nothing references `.frc`. Present in exactly the 7 dirs where MM-GBSA ran, absent from the 3 where it did not. Timestamps match each `mmgbsa_results.dat`. Regenerated on re-run. |
| `artifacts/*` (30) | 198 MB | DUPLICATE | md5-identical to a copy that stays live; the manifest names the survivor for each |
| 2 RL checkpoints + their duplicates (4) | 379 MB | OUTDATED | superseded campaigns; v3 checkpoint retained. ⚠ **no live copy remains** — by decision |
| `filtering/__pycache__/*.pyc` (1) | 8 KB | INTERMEDIATE | bytecode |

**Deletion carried out 27 Sep 2026**: 2.2 GB removed — 7 `reference.frc` byproducts (1.83 GB),
28 verified duplicates (~11 MB, each re-hashed against its surviving copy at deletion time) and
4 RL checkpoints (379 MB, the one irreversible loss: campaigns 1 and 2 can no longer be resumed).
Full record with per-file reasons:
[`provenance/deleted_files_manifest.csv`](provenance/deleted_files_manifest.csv).
Relocation log: [`provenance/restructure_moves.csv`](provenance/restructure_moves.csv).
Rationale: [`provenance/README.md`](provenance/README.md).

## Deliberately kept, against a "minimal" reading

The brief said to remove outdated, duplicate, intermediate and irrelevant files. These look like they
qualify and do not — each is load-bearing for a conclusion:

| Kept | Why removing it would break something |
|---|---|
| `system_TTBK2m/`, `system_TTBK2_p3/`, `system_TTBK1_p2/`, `system_TTBK1_p3/` | **Dissociated and drifted MD runs.** The headline result — *TTBK1 holds the ligand in 1 of 3 poses, TTBK2 in 1 of 3, so there is no stability difference* — exists **only** because these failures were run and retained. Delete them and the central negative finding becomes unsupported. |
| `00_library/.../campaign1_first_failed/` CSVs, plots | A docking-only reward produced **zero** BBB-passing molecules. That failure is what justifies the BBB-aware reward in the final campaign. |
| `06_md/params/cand_003.acpype/` | Superseded free-conformer parameters. Kept so the charge-mismatch incident in the logbook remains demonstrable against the docked-pose set actually used. |
| `06_md/system/mmgbsa_results.dat` | The 1 ns estimate (−34.55) alongside the 10.1 ns result (−31.36); kept to show convergence. |
| `06_md/*/production*.dcd` | **13.3 GB** of raw solvated trajectories across 17 files (was 6.75 GB before the MAO pose scan and velocity replicates; `06_md/` totals 17.5 GB, 91% of the workspace). Retained by decision — all reported numbers come from the stripped `.nc`, but re-stripping with different masks needs these. Gitignored, so this is a disk cost only. |
| `04_docking/*dry*`, `*noFAD*`, `*brg*` | Not redundant variants — they are the **controls** that quantified the water artifact and the cofactor dependency. |
| `01_smiles/references.csv` rows flagged "CONFIRM primary source" | Flagged rather than removed, so the gap stays visible before submission. |

The principle: **negative and superseded results were kept wherever they carry scientific meaning.**
Only true byproducts and verified duplicates were quarantined.

## Verification (rule 8)

- [x] Every final result has a documented method — 16 `README.md` files, one per step
- [x] Every methodological choice has a rationale — each README's *Why we did it* section, cross-referenced to `LOGBOOK.md`
- [x] Required inputs and scripts present — all 17 MD systems retain trajectory, topology, RMSD and energy files; the 3 without MM-GBSA are exactly the 3 that intentionally have none (verified 2026-09-28, after the TTBK2 pose-2 replicate)
- [x] No final result depends on a quarantined file — the only quarantined data-like items are `reference.frc` (referenced by nothing) and verified duplicates
- [x] Structure is navigable by someone who did not do the work — [`README.md`](README.md) walks method → rationale → reference → result
- [x] Scripts still resolve — no numbered directory was renamed or moved

## Structure

`04_docking/` held 77 sibling run directories — 5,160 files, 85% of everything tracked — and is
now grouped by purpose into `candidates/` (24 runs), `references/` (18), `native_redock/` (30) and
`controls/` (5), alongside the untouched `prior_phase/`. The **run directory names are unchanged**,
so every parser still recovers receptor, ligand set and seed from the run name; readers glob one
extra level and the write-side mapping lives only in `scripts/docking_paths.sh`. Verified by
re-running `analyze_selectivity.py`, `analyze_water_test.py`, `calibration_9iv.py` and
`collect_results.py` before and after the move: **stdout, stderr and every output file
byte-identical**, and all three docking drivers report `runs needed: 0` on `--dry-run`, so the
resume logic still recognises the completed work rather than re-docking it.

The numbered step directories were deliberately **not** touched: 127 hardcoded path references
live in `scripts/`, `generation/` and `10_results/`, 48 of them to `06_md` alone.

## Remaining decisions

1. ~~Approve permanent deletion of the quarantine~~ — **done 27 Sep 2026**, 2.2 GB removed.
2. ~~Decide on `_REVIEW_REQUIRED/`~~ — **done**: both project plans kept, moved to `provenance/`.
3. **The manuscript body is still the 20 Aug prior-phase draft** — `4NFM`/`6U0K` receptors, 20 ns,
   the −26.32/−36.85 MM-GBSA pair quoted as SEM, a 15-candidate shortlist, §3.7 unsplit, and a
   resolved limitation (§4.5 item 5, "no replicate MD") still listed. `09_manuscript/README.md`
   states what each claim must become; none of it has been applied to
   `manuscript_dualtarget_flavonol.md`, and the `.html`/`.docx` renders and figure index follow
   from it. This is a scientific task, not a cleanup one, and it is the largest remaining item in
   the project.
4. **Propagate `no_pose` into `08_analysis/selectivity_margins.csv`.**
   `scripts/analyze_selectivity.py` applies its `FAIL_THRESHOLD` when printing but drops the
   `no_pose` column before writing, so the CSV ships all 56 rows unflagged and a plain mean of its
   MAO column returns the withdrawn −3.89 rather than −3.28. That trap had already reached two
   READMEs (logbook entry H, 27 Sep). Flagging the 8 rows closes it at the source.
5. **Vinardo has never been re-run on the validated receptors** — the only cross-check is
   prior-phase, on `4NFM` (`WORKFLOW.md` Phase 3; manuscript §3.3).
6. ~~TTBK2 pose 2 has no velocity replicate~~ — **done 2026-09-28**
   (`06_md/systems/system_TTBK2_p2_r2`, via `scripts/run_ttbk2_replicate.sh`). Every arm in the
   project now has a measured replicate spread. TTBK2's is 3.97 kcal/mol, the ΔΔG fell from 1.45
   to 0.53, and the pose stability reproduced (1.89 Å against 1.99 Å).
