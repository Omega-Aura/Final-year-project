# Provenance

The paper trail for *how* this workspace came to look the way it does. Nothing here is input to
the science; it exists so that any question of the form "where did this come from?" or "what used
to be here?" has a documented answer.

| File | What it is |
|---|---|
| [`original_plan_first_2026-08-17.json`](original_plan_first_2026-08-17.json) | the project plan as first written, 17 Aug 2026 12:07 |
| [`original_plan_revised_2026-08-17.json`](original_plan_revised_2026-08-17.json) | the revision made 18 minutes later, 12:25 |
| [`deleted_files_manifest.csv`](deleted_files_manifest.csv) | every file removed in the 27 Sep cleanup, with its original path, size and the reason |
| [`restructure_moves.csv`](restructure_moves.csv) | the 26 Sep reorganisation: which `artifacts/` file went to which workflow step |

## The two plans

Both predate any results. They are kept because they show the pipeline was specified before it was
run, and because comparing them against [`../LOGBOOK.md`](../LOGBOOK.md) shows where reality
diverged from the plan — which is most of the interesting part of this project. The plan did not
anticipate the water-shell asymmetry, the missing MAO-A cofactor, the off-pose MM-GBSA sign
inversion, or that docking rank would fail to predict pose stability in three of four targets.

## What was deleted, and why it was safe

`deleted_files_manifest.csv` is the full record. In summary, **2.2 GB was removed on 27 Sep 2026**
in three groups:

- **7 × `reference.frc`, 1.83 GB.** Raw force dumps left behind by `MMPBSA.py`. Nothing reads
  them — the only mentions anywhere are `.gitignore` and a note in `06_md/README.md` explaining
  that they are a byproduct — and `scripts/run_mmgbsa.sh` regenerates them.
- **28 duplicate files, ~11 MB.** Copies left in a staging directory. **Each was re-hashed against
  its surviving canonical copy immediately before deletion**, and the manifest names that copy, so
  no file was removed on the strength of an earlier session's claim. Three had been mislabelled as
  having no surviving twin because the canonical copy has a different filename
  (`artifacts/2Z5X.pdb` → `03_receptors/mao/2Z5X/raw.pdb`, and similarly for `7Q8Y.pdb` and
  `phase0_flavonoid_library.csv`); all three were confirmed byte-identical before removal.
- **4 × RL checkpoints, 379 MB.** The stage-1 model weights for generative campaigns 1 and 2, held
  as two unique files plus two exact duplicates. **This is the one deletion that lost something
  irreversible:** those two campaigns can no longer be resumed or re-sampled. Everything needed to
  *read* them survives in [`../00_library/reinvent4_output/`](../00_library/) — per-chunk CSVs,
  score-progression plots, filtered output, combined run tables and the `.toml` configs — so the
  campaigns remain documented. Campaign 3, the one that produced the lead compound `cand_003`,
  keeps its checkpoint at `generation/rl_stage1_v3.chkpt`.

## What was deliberately kept

Negative and superseded results were kept wherever they carry scientific meaning, because in this
project they carry most of it:

| Kept | Why |
|---|---|
| the failed and dissociated MD runs | the "1 of 3 poses holds" finding exists *because* these were run and kept |
| `campaign1_first_failed/` CSVs and plots | a docking-only reward yielding zero BBB-passing molecules is what motivated the BBB-aware reward |
| `06_md/params/cand_003.acpype/` | superseded free-conformer parameters, kept as provenance against the docked-pose set actually used |
| `06_md/system/mmgbsa_results.dat` | the 1 ns estimate, kept to show the 10.1 ns result had converged |
| `system_MAOA/` pose 1 and its −44.22 | the off-pose trajectory that inverted the MAO conclusion — the project's cautionary case |

Trajectories (`*.dcd`, `*.nc`) and topologies (`*.prmtop`) are excluded from version control by
`.gitignore` for size, not deleted; they are on local disk. See
[`../06_md/README.md`](../06_md/README.md) for what that means for reproducing the reported numbers.
