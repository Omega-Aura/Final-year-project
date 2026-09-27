# Dual-Target Flavonol CADD Study — Workflow Index

A computer-aided drug design study of a flavonol scaffold against two target families:
**TTBK1/TTBK2** (casein kinase 1 superfamily, tau kinases) and **MAO-A/MAO-B** (monoamine
oxidases). The lead compound carried through the full pipeline is `cand_003`:

```
Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F
```

Read this file first. Every numbered directory below is one step of the workflow and carries
its own `README.md` answering: **what we did → why → reference → inputs/parameters →
analysis → result → relevant files**.

---

## How to read the workflow

| Step | Directory | Question it answers |
|------|-----------|--------------------|
| 0 | [`00_library/`](00_library/) | Which plant source and which flavonoid chemical space? |
| 1 | [`01_smiles/`](01_smiles/) | Which molecules go into the pipeline, and which references anchor them? |
| 2 | [`02_ligands/`](02_ligands/) | How were those molecules turned into dockable 3D structures? |
| 3 | [`03_receptors/`](03_receptors/) | Which protein structures, and are they fit for docking? |
| 4 | [`04_docking/`](04_docking/) | How do the ligands score and pose in each binding site? |
| 5 | [`05_validation/`](05_validation/) | Does the docking protocol reproduce known answers? |
| 6 | [`06_md/`](06_md/) | Do the docked poses survive unrestrained dynamics? |
| 7 | [`07_mmgbsa/`](07_mmgbsa/) | What is the binding free energy, and is it resolvable? |
| 8 | [`08_analysis/`](08_analysis/) | Consensus ranking, filtering cascade, selectivity margins |
| 9 | [`09_manuscript/`](09_manuscript/) | Write-up and figure/table plan |
| 10 | [`10_results/`](10_results/) | **What the study concludes** — results, figures, and what can and cannot be claimed |

Supporting directories:

| Directory | Role |
|-----------|------|
| [`generation/`](generation/) | REINVENT4 reinforcement-learning ligand generation (campaign v3) |
| [`filtering/`](filtering/) | BOILED-Egg / BBB permeability gate |
| [`scripts/`](scripts/) | All driver scripts; every step is re-runnable from here |
| `REINVENT4/` | Third-party clone, installed locally (not vendored, untracked) |
| `_ARCHIVE_TO_DELETE/` | Quarantined leftovers awaiting deletion approval — see its README |
| `_REVIEW_REQUIRED/` | Items whose fate needs a human decision |

**[`LOGBOOK.md`](LOGBOOK.md)** is the chronological record and the authoritative account of
*why* each decision was made, including every defect found and corrected. Where a step README
and the logbook disagree, the logbook is the primary source.

---

## Two phases, one workflow

This workspace contains two passes over the same scientific question:

- **Phase 0–8 (17–21 Aug)** — the original pipeline, docked against `4NFM`, `6U0K`, `2V60`
  alongside `7Q8Y`/`2V5Z`/`2Z5X`. Produced the manuscript and the `phase*_report.md` series.
- **Week 1–2 (24–26 Sep)** — the pipeline rebuilt on properly validated structures
  (`7JXX`, `7Q8Y`, `2V5Z`, `2Z5X`), with redocking validation, a water-shell control, MD and
  MM-GBSA.

They are organised **by workflow step, not by phase**. Each step directory holds its current
work at the top level and the same step's earlier outputs in a `prior_phase/` subfolder.

`prior_phase/` is kept separate rather than merged because the two passes used *different
receptor structures*. Mixing results from `4NFM` and `7JXX` in one directory is precisely the
kind of silent mismatch this project has already been bitten by more than once (see the
water-asymmetry and ligand-charge-mismatch entries in the logbook). The separation is a
safety rail, not an archive.

---

## Headline results

Summarised below. [`10_results/`](10_results/) is the fuller account: the same numbers
organised by **claim**, each with its status, plus the figures and the limitations that bound
them. Read that directory before writing any of this up.

### Molecular dynamics — pose stability

Ligand RMSD from the docked pose, mean over the final 100 frames of 10.1 ns:

| System | Pose | Last-100 RMSD | Verdict |
|---|---|---|---|
| TTBK1 (7JXX) | 1, run 1 | 1.88 Å | stable |
| TTBK1 (7JXX) | 1, run 2 | 1.38 Å | stable |
| TTBK1 (7JXX) | 2 | 4.18 Å | drifts |
| TTBK1 (7JXX) | 3 | 3.42 Å | drifts |
| TTBK2 (7Q8Y) | 1, run 1 | 3.72 Å | drifts |
| TTBK2 (7Q8Y) | 1, run 2 | 5.69 Å | dissociates |
| TTBK2 (7Q8Y) | 2 | 1.99 Å | stable |
| TTBK2 (7Q8Y) | 3 | 5.79 Å | dissociates |
| MAO-A (2Z5X) | 1 | 3.24 Å | drifts off docked pose |
| MAO-A (2Z5X) | 2 | 2.43 Å (1.79 core-fit) | holds, loose |
| MAO-A (2Z5X) | 3 | 1.43 Å | stable |
| MAO-B (2V5Z) | 1 | 2.28 Å | stable |
| MAO-B (2V5Z) | 2 | 1.10 Å | stable |
| MAO-B (2V5Z) | 3 | 2.57 Å | holds, loose |

**TTBK1 holds the ligand in 1 of 3 docked poses; TTBK2 also in 1 of 3. There is no
pose-stability difference between the two kinases.** Protein backbones were stable in every
run (1.12–1.81 Å, and 1.06–1.63 Å core-fit for MAO), so the ligand motion is ligand motion,
not a collapsing binding site.

**The MAO complexes are the more stable pair.** No MAO run left the site — the worst is 3.24 Å,
against two TTBK runs past 5 Å. Counting strictly, MAO-A holds 1 of 3 and MAO-B 2 of 3.

**Docking rank does not predict pose stability**, and with 3 poses × 4 targets this is now the
project's best-supported methodological result: the top-ranked docked pose was the most stable
one in **one of four** targets (TTBK1, the exception). The best-holding MAO-A pose is its
*worst*-ranked.

### MM-GBSA — binding free energy

| System | ΔG (kcal/mol) | SD | Usable for selectivity? |
|---|---|---|---|
| TTBK1 pose 1 run 1 | −31.36 | 3.82 | yes |
| TTBK1 pose 1 run 2 | −33.44 | 2.63 | yes |
| TTBK2 pose 2 | −30.95 | 3.81 | yes |
| MAO-A pose 1 | −44.22 | 3.15 | **no** — ligand left the docked pose |
| MAO-A pose 2 | −35.54 | 2.73 | direction only |
| **MAO-A pose 3** | **−36.67** | 2.40 | yes — best on-pose MAO-A |
| MAO-B pose 1 | −39.48 | 2.20 | yes |
| **MAO-B pose 2** | **−39.20** | 2.13 | yes — best on-pose MAO-B |
| MAO-B pose 3 | −39.04 | 2.86 | direction only |

No MM-GBSA was computed for the runs that drifted or dissociated: an energy averaged over a
trajectory that has left the docked pose describes a structure that was never docked.

**The −44.22 is the cautionary number of this project.** It is the most favourable MAO-A value,
it is the only off-pose one, and taking it at face value inverted the MAO selectivity conclusion
for two days. Every on-pose MAO-A trajectory lands near −36. Corrected, all six on-pose pairings
favour **MAO-B** by 2.37–3.94 kcal/mol, which agrees in sign with docking (3.03–3.22 pose for
pose). The direction is robust; the magnitude sits inside the replicate spread and is not.

### What the dynamics does and does not establish

| Evidence | TTBK1 vs TTBK2 | MAO-A vs MAO-B |
|---|---|---|
| Docking, water-symmetric | TTBK2 favoured by ~1.6 kcal/mol; 0/56 candidates favourable | MAO-B favoured, −3.89 mean, 48/48 favourable |
| Docking, validated pair | +0.12, confounded by a 0-vs-5 water asymmetry | — |
| Protocol bias (known-answer calibration) | 0.861 kcal/mol toward TTBK2 where experiment says none | — |
| MD pose stability | **no difference** — 1/3 poses each, symmetric test | MAO-A 1/3, MAO-B 2/3; no run left the site |
| MM-GBSA | **no resolvable difference** (1.45 < 2.08 replicate spread) | MAO-B favoured 2.37–3.94, **direction only** |

The two pairs now read differently, and the write-up must not treat them as one finding.

**TTBK1 vs TTBK2 — no discrimination.** The manuscript's TTBK2 liability claim rests on docking
scores alone. Dynamics neither confirms nor refutes it; it lacks the resolution to speak, and that
is the honest outcome.

**MAO-A vs MAO-B — a direction, not a number.** Docking and MM-GBSA independently agree that
MAO-B is favoured. But the margin is comparable to the 2.08 kcal/mol replicate spread measured on
TTBK1, so "no *significant* discrimination" still holds even though "no discrimination" no longer
does. Velocity replicates of the two best on-pose MAO systems are running to measure the MAO arms'
own spread directly.

### Two methodological findings worth more than the numbers

1. **The SEM is not the error bar.** `MMPBSA.py` computes the standard error as if 200 frames
   sampled 10 ps apart were independent draws. Reported SEMs run 0.15–0.27 kcal/mol while the
   true replicate-to-replicate spread is ~2.08 — about 10× larger. Quoting the SEM on a ΔΔG is
   the easiest way to manufacture a significant selectivity result from this pipeline.
2. **Both arms of a comparison must be symmetric.** Three separate false positives this
   project produced all had the same root cause: an asymmetry between the two things being
   compared (a 0-vs-5 water shell; a pose scan against a replicate scan; an on-pose run against
   an off-pose run). Checking arm symmetry before interpreting is now a standing pre-flight step.

---

## Reproducing the pipeline

Environments (see `LOGBOOK.md` for provenance):

| Environment | Holds | Used for |
|---|---|---|
| `docking_project` (conda, Windows) | vina, meeko, rdkit, openmm+CUDA, spyrmsd, obabel | docking, MD, analysis |
| `mdgbsa` (conda, WSL2) | AmberTools, acpype | GAFF2 parameterisation, `MMPBSA.py`, `cpptraj`, `tleap` |
| `reinvent4` (conda) | REINVENT4 | RL ligand generation |

AmberTools has no native Windows build, so parameterisation and MM-GBSA run in WSL2 while
OpenMM runs natively on the GPU. Scripts that cross this boundary (`scripts/run_mmgbsa.sh`)
invoke `wsl -e bash -lc` and translate paths to `/mnt/c/...`.

**The numbered directory names are load-bearing.** Over sixteen hardcoded references to
`03_receptors/`, `01_smiles/`, `04_docking/`, `02_ligands/pdbqt/`, `06_md/system/`,
`05_validation/` and `08_analysis/` live in `scripts/` and `generation/`. Renaming a step
directory breaks reproducibility; add to the structure rather than rearranging it.
