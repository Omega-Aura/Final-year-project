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
| **REVIEW_REQUIRED** | could not be classified from content | parked in `_REVIEW_REQUIRED/` |

## Where everything stands

| Directory | Size | Files | Classification |
|---|---|---|---|
| `00_library/` | 1.5 MB | 26 | REFERENCE + SUPPORTING — curated library with PMIDs; RL campaign history |
| `01_smiles/` | 56 KB | 8 | **FINAL** — the frozen candidate and reference sets |
| `02_ligands/` | 792 KB | 165 | SUPPORTING — prepared 3D inputs, regenerable but cheap to keep |
| `03_receptors/` | 15 MB | 91 | **FINAL** + REFERENCE — docking targets, controls, characterisation |
| `04_docking/` | 101 MB | 5,160 | **FINAL** — all production docking, incl. wet/dry and cofactor controls |
| `05_validation/` | 126 KB | 11 | **FINAL** — redocking and calibration; the protocol's credibility rests here |
| `06_md/` | 14 GB | 582 | **FINAL** + SUPPORTING — 14 MD systems, raw + stripped trajectories |
| `10_results/` | <1 MB | 5 | **FINAL** — consolidated results, interpretation, 3 figures + their generator |
| `07_mmgbsa/` | 296 KB | 7 | **FINAL** — collected energy summary + prior-phase results |
| `08_analysis/` | 700 KB | 33 | **FINAL** — consensus, cascade, selectivity, water test |
| `09_manuscript/` | 5.6 MB | 6 | **FINAL** — the deliverable |
| `filtering/` | 20 KB | 3 | SUPPORTING — the BBB model, used by both cascade and RL reward |
| `generation/` | 95 MB | 12 | **FINAL** + SUPPORTING — v3 campaign, config and checkpoint |
| `scripts/` | 116 KB | 20 | **FINAL** — every step is re-runnable from here |
| `REINVENT4/` | 1.1 GB | 3,679 | third-party clone, untracked; required to re-run generation |
| `_ARCHIVE_TO_DELETE/` | 2.1 GB | 43 | quarantined, awaiting deletion approval |
| `_REVIEW_REQUIRED/` | 40 KB | 3 | parked pending decision |

Workspace total: **11 GB**, of which 2.1 GB is quarantined.

## Dependency check (rule 5) — what it changed

Before moving anything, every candidate was checked for being referenced by a script, needed to
reproduce a result, or the only copy of something. Two findings redirected the whole approach:

### The numbered directories are load-bearing — so they were not renamed

`scripts/` and `generation/` hardcode the step directories in **16+ places**: `03_receptors/`,
`01_smiles/`, `04_docking/`, `02_ligands/pdbqt/`, `06_md/system/`, `05_validation/`,
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

Full record with per-file reasons: [`_ARCHIVE_TO_DELETE/MANIFEST.csv`](_ARCHIVE_TO_DELETE/MANIFEST.csv).
Relocation log: [`_ARCHIVE_TO_DELETE/RESTRUCTURE_MOVES.csv`](_ARCHIVE_TO_DELETE/RESTRUCTURE_MOVES.csv).

## Deliberately kept, against a "minimal" reading

The brief said to remove outdated, duplicate, intermediate and irrelevant files. These look like they
qualify and do not — each is load-bearing for a conclusion:

| Kept | Why removing it would break something |
|---|---|
| `system_TTBK2m/`, `system_TTBK2_p3/`, `system_TTBK1_p2/`, `system_TTBK1_p3/` | **Dissociated and drifted MD runs.** The headline result — *TTBK1 holds the ligand in 1 of 3 poses, TTBK2 in 1 of 3, so there is no stability difference* — exists **only** because these failures were run and retained. Delete them and the central negative finding becomes unsupported. |
| `00_library/.../campaign1_first_failed/` CSVs, plots | A docking-only reward produced **zero** BBB-passing molecules. That failure is what justifies the BBB-aware reward in the final campaign. |
| `06_md/params/cand_003.acpype/` | Superseded free-conformer parameters. Kept so the charge-mismatch incident in the logbook remains demonstrable against the docked-pose set actually used. |
| `06_md/system/mmgbsa_results.dat` | The 1 ns estimate (−34.55) alongside the 10.1 ns result (−31.36); kept to show convergence. |
| `06_md/*/production*.dcd` | 6.75 GB of raw solvated trajectories. Retained by decision — all reported numbers come from the stripped `.nc`, but re-stripping with different masks needs these. |
| `04_docking/*dry*`, `*noFAD*`, `*brg*` | Not redundant variants — they are the **controls** that quantified the water artifact and the cofactor dependency. |
| `01_smiles/references.csv` rows flagged "CONFIRM primary source" | Flagged rather than removed, so the gap stays visible before submission. |

The principle: **negative and superseded results were kept wherever they carry scientific meaning.**
Only true byproducts and verified duplicates were quarantined.

## Verification (rule 8)

- [x] Every final result has a documented method — 16 `README.md` files, one per step
- [x] Every methodological choice has a rationale — each README's *Why we did it* section, cross-referenced to `LOGBOOK.md`
- [x] Required inputs and scripts present — all 14 MD systems retain trajectory, topology, RMSD and energy files; the 3 without MM-GBSA are exactly the 3 that intentionally have none (verified 2026-09-27, after the MAO pose scan added 4 systems)
- [x] No final result depends on a quarantined file — the only quarantined data-like items are `reference.frc` (referenced by nothing) and verified duplicates
- [x] Structure is navigable by someone who did not do the work — [`README.md`](README.md) walks method → rationale → reference → result
- [x] Scripts still resolve — no numbered directory was renamed or moved

## Remaining decisions

1. **Approve permanent deletion** of `_ARCHIVE_TO_DELETE/` (2.1 GB) — see its README first.
2. **Decide on `_REVIEW_REQUIRED/`** — two 36 KB original project-plan JSONs.
3. **Three manuscript claims need revision** against the Week 1–2 data (§3.5, §3.6, §3.7) — see
   [`09_manuscript/README.md`](09_manuscript/README.md). This is a scientific task, not a cleanup one,
   but it is the most consequential item on this list.
