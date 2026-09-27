# Step 5 — Protocol Validation

## What we did

Before trusting any candidate score, tested whether the docking protocol can reproduce answers
that are already known — by redocking each receptor's own crystallographic ligand back into its
own site, and by checking the scoring function against measured MAO affinities.

## Why we did it

A docking score is only meaningful if the protocol can recover a pose that was determined
experimentally. Redocking is the cheapest available falsification test: if Vina cannot place a
native ligand back where crystallography found it, nothing it says about a novel candidate in
that site is worth reporting.

The MAO benchmark exists for the same reason at the level of *ranking* rather than *geometry* —
a protocol can place poses well and still rank affinities backwards.

## Reference

- **Redocking criterion** — heavy-atom RMSD to the crystallographic pose, ≤ 2.0 Å is the
  conventional success threshold
- **RMSD engine** — `spyrmsd` (symmetry-corrected; a plain atom-order RMSD inflates values for
  symmetric groups)
- **Docking protocol** — identical to [step 4](../04_docking/): AutoDock Vina,
  `--exhaustiveness 32 --num_modes 9`, three seeds (11, 22, 33)
- **Experimental MAO data** — IC50 values sourced via PubMed for kaempferol, quercetin,
  lazabemide and isatin; see [`../LOGBOOK.md`](../LOGBOOK.md) and
  [`../01_smiles/references.csv`](../01_smiles/references.csv)

## Inputs and parameters

Driven by [`../scripts/rmsd_check.py`](../scripts/rmsd_check.py) and
[`../scripts/calibration_9iv.py`](../scripts/calibration_9iv.py). Native ligand references come
from each receptor directory (`03_receptors/<PDB>/native_*_ref.sdf`), extracted from the
deposited structure.

Bond orders in native ligands extracted from PDB files are unreliable; they are repaired by
[`../scripts/fix_native_bondorders.py`](../scripts/fix_native_bondorders.py) before RMSD, since a
wrong bond order changes the symmetry graph and therefore the RMSD.

## Analysis performed

**Redocking, all six receptors** — best pose across three seeds:

| Receptor | Protein | Native ligand | Best RMSD | Result |
|---|---|---|---|---|
| 2Z5X | MAO-A | HRM (harmine) | **0.27 Å** | PASS |
| 7JXX | TTBK1 | VP7 | **0.55 Å** | PASS |
| 2V5Z | MAO-B | SAG (safinamide) | **0.66 Å** | PASS |
| 7Q8V | TTBK1 | 9IV | **0.68 Å** | PASS |
| 4BTK | — | DTQ | **0.79 Å** | PASS |
| 7Q8Y | TTBK2 | 9IV | **1.45 Å** | PASS |

**Scoring-function calibration on the 9IV pair** — 9IV is a ligand with measured affinity for
both TTBK1 and TTBK2, so the docking ΔΔG can be compared directly against experiment:

| Quantity | Value |
|---|---|
| Docking consensus, 7Q8V (TTBK1) | −8.115 ± 0.015 |
| Docking consensus, 7Q8Y (TTBK2) | −8.899 ± 0.022 |
| **Docking margin** | **0.784 kcal/mol favouring TTBK2** |
| **Experimental ΔΔG** | **−0.234 kcal/mol** (TTBK1 IC50 330 nM vs TTBK2 490 nM, same assay) |
| **Systematic bias** | **1.018 kcal/mol** |

## Final result

**The geometry is trustworthy. The selectivity margin is not.**

All six receptors pass redocking, four of them under 0.7 Å. The protocol reliably reproduces
crystallographic poses, so pose-level conclusions rest on solid ground.

The 9IV calibration tells a different story. Experiment says 9IV has **essentially no preference**
between TTBK1 and TTBK2 (ΔΔG −0.234, a slight preference for TTBK1). The docking protocol
reports a 0.784 kcal/mol preference for TTBK2 — a **1.018 kcal/mol systematic bias toward
TTBK2** on a case
where the true answer is known to be "no difference".

This matters directly for the study's central claim. Any TTBK2-over-TTBK1 margin smaller than
~1.0 kcal/mol is within the protocol's demonstrated bias and cannot be read as selectivity. The
MAO benchmark in [`benchmark_mao.csv`](benchmark_mao.csv) plays the same role for the
flavoenzymes.

## Relevant files

| File | Role |
|---|---|
| `*_redock.txt` | per-receptor redocking, all seeds and poses, with the PASS/FAIL verdict |
| [`calibration_9IV.csv`](calibration_9IV.csv) | per-seed consensus scores for the 9IV pair |
| [`calibration_9IV_margin.json`](calibration_9IV_margin.json) | the margin, the experimental ΔΔG, and the derived systematic bias |
| [`benchmark_mao.csv`](benchmark_mao.csv) | MAO reference-ligand scores vs measured affinity |
| [`benchmark_mao_scatter.png`](benchmark_mao_scatter.png) | predicted vs experimental, MAO benchmark |
| [`../scripts/rmsd_check.py`](../scripts/rmsd_check.py) | redocking RMSD driver |
| [`../scripts/calibration_9iv.py`](../scripts/calibration_9iv.py) | calibration driver |
| [`../scripts/fix_native_bondorders.py`](../scripts/fix_native_bondorders.py) | native-ligand bond-order repair |
