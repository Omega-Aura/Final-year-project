# Quarantine — Awaiting Deletion Approval

**Nothing here has been deleted.** These files were *moved* out of the workspace, keeping their
original directory layout, so anything can be put back. Review, then approve permanent deletion.

**42 files, 2.1 GB.**

[`MANIFEST.csv`](MANIFEST.csv) records every move: `original_path`, `archive_path`, `category`,
`bytes`, and the `reason` it was quarantined. [`RESTRUCTURE_MOVES.csv`](RESTRUCTURE_MOVES.csv) is a
separate log of the `artifacts/` reorganisation (those files were *relocated into the workflow*, not
quarantined).

## What is in here, and why

### `06_md/*/reference.frc` — 7 files, 1.83 GB (INTERMEDIATE)

Raw forces dumps left behind by `MMPBSA.py`. Each line reads
`Forces on <n> at <x> <y> <z> - <fx> <fy> <fz>`.

Evidence they are leftovers, not data:

- **Nothing reads them.** No script, cpptraj input, `mmpbsa.in`, or notebook references `.frc`.
- **They exist in exactly the 7 directories where MM-GBSA ran** and in none of the 3 where it did
  not (`system_TTBK1_p2`, `system_TTBK1_p3`, `system_TTBK2_p3`).
- **Each file's timestamp matches its `mmgbsa_results.dat`** to the minute.
- They are regenerated on any re-run of `scripts/run_mmgbsa.sh`.

No reported number depends on them. This is the bulk of the reclaimed space.

### `artifacts/*` — 30 files, 198 MB (DUPLICATE)

Byte-identical copies, verified by md5, of files that remain live elsewhere. The `reason` column in
the manifest names the surviving canonical copy for each one, for example:

| Quarantined | Canonical copy still in place |
|---|---|
| `artifacts/manuscript_dualtarget_flavonol.{md,html,docx}` | `09_manuscript/` |
| `artifacts/phase0_flavonoid_library.csv` | `00_library/flavonoid_library.csv` |
| `artifacts/2Z5X.pdb`, `artifacts/7Q8Y.pdb` | `03_receptors/<PDB>/raw.pdb` |
| `artifacts/pdbqt_utils.py`, `artifacts/run_prolif_v2.py` | `scripts/` |
| `artifacts/rl_*chunk*.csv`, `rl_v2_shortlist_56.csv`, `rl_v2_top20_shortlist.csv` | `00_library/reinvent4_output/campaign*/` |

`artifacts/` was a staging mirror. Its 144 *unique* files were **not** quarantined — they were moved
into the workflow steps they belong to (see `RESTRUCTURE_MOVES.csv` and each step's
`prior_phase/`).

### RL checkpoints — 4 files, 379 MB (OUTDATED + DUPLICATE)

`campaign1_first_failed/rl_stage1_chunk1.chkpt` and `campaign2_v2/rl_v2_stage1_chunk1.chkpt`, plus
their `artifacts/` duplicates.

⚠ **For these two checkpoints, no copy remains in the live tree.** That was the explicit decision:
the campaign that produced the final candidate set is v3, whose checkpoint
(`generation/rl_stage1_v3.chkpt`) is retained.

What stays live is everything needed to *read* those campaigns' history — score-progression plots,
per-chunk CSVs, filtered output, combined run tables and the `.toml` configs, all still in
`00_library/reinvent4_output/`. Only the 95 MB model weights are gone, so the campaigns remain
documented but cannot be resumed or re-sampled without re-running them.

### `filtering/__pycache__/*.pyc` — 1 file (INTERMEDIATE)

Python bytecode, regenerated on import.

## What was deliberately NOT quarantined

| Kept | Size | Why |
|---|---|---|
| `06_md/*/production*.dcd` | 6.75 GB | raw solvated trajectories — primary data; retained by decision so trajectories can be re-stripped with different masks later |
| `REINVENT4/` | 1.1 GB | third-party clone; `generation/`'s config references `REINVENT4/priors/libinvent.prior`, so it is needed to re-run generation |
| `06_md/*_stripped.nc` | 665 MB | **every reported MD and MM-GBSA number comes from these** |
| Failed / dissociated MD runs | — | `system_TTBK2m`, `system_TTBK2_p3`, `system_TTBK1_p2/p3` are **negative results the main conclusion depends on**. The "1 of 3 poses" finding is only possible because these were run and kept. |
| `campaign1_first_failed/` CSVs and plots | — | the docking-only reward producing zero BBB-passing molecules is what motivated the BBB-aware reward. A kept failure. |
| `06_md/params/cand_003.acpype/` | — | superseded free-conformer parameters, kept for provenance against the docked-pose set actually used |
| `06_md/system/mmgbsa_results.dat` | — | the 1 ns estimate, kept to show the 10.1 ns result converged |

The pattern: **negative and superseded results were kept wherever they carry scientific meaning.**
Only true byproducts and verified duplicates were quarantined.

## To restore something

Move it back to the path in the manifest's `original_path` column, e.g.

```bash
mv "_ARCHIVE_TO_DELETE/06_md/system_MAOA/reference.frc" "06_md/system_MAOA/reference.frc"
```

## To delete permanently

Only after review, and only with explicit approval:

```bash
rm -rf "_ARCHIVE_TO_DELETE"
```

Keep `MANIFEST.csv` somewhere first if you want a record of what was removed.
