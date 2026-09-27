# Step 6 — Molecular Dynamics

## What we did

Ran 10.1 ns of unrestrained all-atom MD on fourteen explicitly solvated protein–ligand complexes:
`cand_003` docked into TTBK1, TTBK2, MAO-A and MAO-B, covering **three independent docked poses
per target** and two velocity-seed replicates of the TTBK top poses.

| Directory | Target | Structure | Start | Ligand RMSD (last 100 frames) | Verdict |
|---|---|---|---|---|---|
| `system/` | TTBK1 | 7JXX | pose 1, run 1 | 1.88 Å | stable |
| `system_TTBK1_r2/` | TTBK1 | 7JXX | pose 1, run 2 | 1.38 Å | stable |
| `system_TTBK1_p2/` | TTBK1 | 7JXX | pose 2 | 4.18 Å | drifts |
| `system_TTBK1_p3/` | TTBK1 | 7JXX | pose 3 | 3.42 Å | drifts |
| `system_TTBK2/` | TTBK2 | 7Q8Y | pose 1, run 1 | 3.72 Å | drifts |
| `system_TTBK2m/` | TTBK2 | 7Q8Y | pose 1, run 2 | 5.69 Å | dissociates |
| `system_TTBK2_p2/` | TTBK2 | 7Q8Y | pose 2 | 1.99 Å | stable |
| `system_TTBK2_p3/` | TTBK2 | 7Q8Y | pose 3 | 5.79 Å | dissociates |
| `system_MAOA/` | MAO-A | 2Z5X | pose 1 | 3.24 Å | drifts off pose |
| `system_MAOA_p2/` | MAO-A | 2Z5X | pose 2 | 2.43 Å (1.79 core-fit) | holds, loose |
| `system_MAOA_p3/` | MAO-A | 2Z5X | pose 3 | 1.43 Å | stable |
| `system_MAOB/` | MAO-B | 2V5Z | pose 1 | 2.28 Å | stable |
| `system_MAOB_p2/` | MAO-B | 2V5Z | pose 2 | 1.10 Å | stable |
| `system_MAOB_p3/` | MAO-B | 2V5Z | pose 3 | 2.57 Å | holds, loose |

Every directory is self-contained: its own protein PDB, its own ligand parameters, its own
`build.leap`, topology, trajectory and analysis. That redundancy is deliberate — see *Why*.

### MAO pose scan (complete)

The MAO pair was **symmetric in design** — both arms started from pose 1 of the *dry*
(water-symmetric) docking at seed 11, confirmed by exact coordinate match — but **asymmetric in
outcome**: MAO-B held its pose at 2.28 Å while MAO-A drifted to 3.24 Å. Because MM-GBSA on an
off-pose trajectory describes a structure that was never docked, the MAO ΔΔG was uninterpretable
until each arm had an on-pose trajectory.

Four runs extended both arms to three top-ranked docked poses each, rank-symmetric from the start
rather than after a false positive. **Both arms now hold**: MAO-A pose 3 at 1.43 Å and MAO-B
pose 2 at 1.10 Å, tighter than either pose-1 run. The energies and what they license are in
[`../07_mmgbsa/README.md`](../07_mmgbsa/README.md); in short, the sign of the MAO margin flipped
to agree with docking once the off-pose MAO-A pose-1 trajectory was set aside.

Pose spacing from pose 1 (all-atom):

| | pose 2 | pose 3 | Vina affinity (pose 1 / 2 / 3) |
|---|---|---|---|
| MAO-A (2Z5Xdry) | 2.90 Å | 8.46 Å | −8.236 / −8.228 / −7.399 |
| MAO-B (2V5Zdry) | 2.41 Å | 8.07 Å | −11.42 / −11.26 / −10.62 |

The two scans are well matched in spacing, so "pose 2" and "pose 3" mean comparable things on
both proteins — which is what makes the arms comparable.

All six MAO systems are neutral at +0.0010 e with 0 tleap errors and per-protein-identical
warning counts (MAO-A 3, MAO-B 2). They differ in **starting ligand pose alone** — same protein
file, same FAD lib/frcmod, same ligand parameters, same build recipe, all copied byte-for-byte
from the pose-1 system.

Queue drivers, both one job at a time and both resumable by skipping finished work:
[`../scripts/run_mao_pose_queue.sh`](../scripts/run_mao_pose_queue.sh) runs the MD and skips any
system that already has `final_state.xml`;
[`../scripts/run_mao_pose_analysis.sh`](../scripts/run_mao_pose_analysis.sh) runs the MM-GBSA and
cpptraj chain that follows and skips any system that already has `mmgbsa_results.dat`.

## Why we did it

**Docking scores rank poses; they do not tell you whether a pose is physically stable.** The
docking stage produced a TTBK2-over-TTBK1 margin that the manuscript interpreted as a
selectivity liability. MD was run to test whether that margin survives dynamics, and whether
the docked geometry is even a real minimum.

Three specific design choices each exist because of a defect found earlier in this project:

1. **Three poses per protein, not one.** An initial run tested TTBK2 from three poses but TTBK1
   from only one pose (twice). That produced an apparent "TTBK2 is less stable" result which was
   purely an artifact of comparing a *pose scan* against a *replicate scan*. Running TTBK1 poses
   2 and 3 made the test symmetric — and the difference vanished.
2. **Two velocity replicates of the top pose.** One replicate per target would have made a
   1.45 kcal/mol ΔΔG look like a five-sigma result against the reported SEM. The replicate spread
   (~2.08 kcal/mol) is the real error bar, and it is only knowable from replicates.
3. **Parameters generated from the docked pose, not a free conformer.** `params/cand_003.acpype/`
   was built from a relaxed free conformer; `params/docked/cand_003.acpype/` was rebuilt from the
   actual docked pose. Only the latter is used, so the MD starting geometry is the geometry the
   docking produced.
4. **Alternative poses transferred by coordinate swap, never re-embedded.**
   [`../scripts/make_pose_pdb.py`](../scripts/make_pose_pdb.py) writes pose *N* onto the acpype
   template positionally, because the docking SDF and the template share atom count and atom
   order. It refuses to proceed unless pose 1 round-trips onto the template (it does, at
   0.0002 Å) and the element sequences agree — so a silent atom-order mismatch cannot produce a
   scrambled ligand that still builds cleanly.

## Reference

Methods and parameter sources, all recorded with dates in [`../LOGBOOK.md`](../LOGBOOK.md):

- **Protein force field** — AMBER `ff14SB`
- **Ligand & cofactor force field** — `GAFF2`, charges via `acpype` (AmberTools)
- **Water** — TIP3P, explicit solvent, neutralising Cl⁻/Na⁺ to 0 net charge
- **Engine** — OpenMM with the CUDA platform (native Windows, RTX 4050)
- **System building** — `tleap` (AmberTools, WSL2 `mdgbsa` env)
- **Protein preparation** — `pdb4amber` for heavy-atom cleanup and gap detection, HIS→HIE,
  hydrogens added with OpenMM `Modeller.addHydrogens`
- **FAD cofactor** — parameterised as an 84-atom GAFF2 residue from RCSB ideal chemistry via
  [`../scripts/prep_cofactor.py`](../scripts/prep_cofactor.py), positionally restrained to its
  crystallographic coordinates. Zero guessed parameters.

Prior-phase MD on a different (apo) structure is in `prior_phase/`; its conclusions are
superseded by the runs above, which use validated holo structures.

## Inputs and parameters

Built by [`run_md_system.py`](run_md_system.py); run by [`run_md.py`](run_md.py),
[`run_md_extend.py`](run_md_extend.py) (continuation from `final_state.xml`) and
[`run_md_restrained.py`](run_md_restrained.py) (FAD-restrained MAO systems).

```
nonbonded method     PME
nonbonded cutoff     1.0 nm
constraints          HBonds, rigid water
integrator           LangevinMiddle, 300 K, 1 ps-1 friction, 2 fs timestep
pressure control     MonteCarloBarostat (active through equilibration and production)
minimisation         2000 steps
equilibration        100 ps (50,000 steps)
production           10.1 ns (5,050,000 steps), frame every 10 ps -> 1000 frames
```

System sizes and throughput:

| System | Atoms | Throughput |
|---|---|---|
| TTBK1 / TTBK2 | ~40,000 | 156–177 ns/day |
| MAO-A (2Z5X) | 118,896–118,923 | 59.8–63.1 ns/day (3 runs) |
| MAO-B (2V5Z) | 90,129 | 80.9–82.0 ns/day (3 runs) |

All fourteen runs completed to step 5,050,000 with temperature stable at ~300 K and potential
energy stable. `tleap` reported zero errors on every build; all systems neutral.

One interruption, no lost science: a machine restart at 16:06 on 2026-09-26 killed the queue
0.77 ns into `system_MAOA_p3`. `run_md_restrained.py` writes no intermediate checkpoint, so that
partial run could not be resumed and the system was simply re-run from the top; the queue skipped
`system_MAOA_p2`, which had already written `final_state.xml`. The completed run is a normal full
10.1 ns trajectory with nothing stitched together.

## Analysis performed

Per system, via `cpptraj` (WSL2 `mdgbsa` env):

- **`strip_traj.cpptraj`** — `autoimage`, then strip `:WAT,Na+,Cl-` → `production_full_stripped.nc`
  plus `stripped.complex.prmtop`. All downstream analysis uses these.
- **`lig_rmsd.cpptraj`** — backbone fit on protein `@CA,C,N`, then mass-weighted `nofit` RMSD of
  the ligand (`:MOL`) and, for the MAO systems, of `:FAD`. Outputs `lig_rmsd.dat`,
  `prot_rmsd.dat`, `fad_rmsd.dat` (1000 frames each).

MAO systems only — both constructs end in a solvent-exposed C-terminal tail that has no membrane
to sit in here, and it dominates any whole-protein fit:

- **`rmsf.cpptraj`** — per-residue `atomicfluct` on `@CA`, which is what identified the tail as
  the culprit: the maximum is always the last residues of the construct (MAO-A 510–513 of 513,
  MAO-B 496–499 of 499) at 5–8 Å.
- **`core_rmsd.cpptraj`** — refits on the ordered core `:1-496@CA,C,N` and re-measures ligand and
  FAD against that. Outputs `prot_core_rmsd.dat`, `lig_corefit_rmsd.dat`, `fad_corefit_rmsd.dat`.
  Without this, `system_MAOA_p2` reads 3.44 Å backbone RMSD — the highest in the project — and
  looks like an unstable fold; its core is 1.28 Å and its ligand 1.79 Å rather than 2.43 Å.

Stability was judged on the **mean of the final 100 frames**, not the whole-run mean: a ligand
that leaves late still shows a low whole-run average. The collected table is built from these
files by [`../scripts/collect_md_summary.py`](../scripts/collect_md_summary.py) rather than
transcribed by hand.

## Final result

**No pose-stability difference between the kinases.** TTBK1 holds `cand_003` in 1 of 3 docked
poses; TTBK2 also in 1 of 3. Which pose survives differs (TTBK1 pose 1, TTBK2 pose 2), but that
is not a selectivity signal — it says the docking pose ranking does not predict dynamic
stability, equally for both proteins. Protein backbones were stable throughout every run
(1.12–1.81 Å, and 1.06–1.63 Å core-fit for the MAO systems).

**Docking rank does not predict pose stability in the MAO pair either.** The pose that holds best
is MAO-A's *worst*-ranked pose 3 (1.43 Å, Vina −7.399 against pose 1's −8.236) and MAO-B's pose 2
(1.10 Å, Vina −11.26 against −11.42). Across the four targets, the top-ranked docked pose was the
most stable one in **one of four** cases — TTBK1, the exception. With 3 poses × 4 targets now run,
this is the project's best-supported methodological result, and it is the reason every MM-GBSA
number is tied to a named pose rather than to a target.

**The FAD restraint scheme is validated six times.** FAD stayed within 0.46–0.73 Å of its
crystallographic position across all six MAO systems, using the same lib and frcmod with only
coordinates changed. The flavin wall of the cavity was present and rigid as intended, which
removes the blocker that made MAO MM-GBSA look like multi-day work.

**The two MAO caveats that blocked a selectivity claim are now resolved** — both arms have an
on-pose trajectory, and MM-GBSA agrees with docking on the sign once the off-pose MAO-A pose-1
trajectory is excluded. What replaces them is a narrower caveat about magnitude, since neither MAO
arm has a velocity replicate. See [`../07_mmgbsa/README.md`](../07_mmgbsa/README.md), where the
energies live.

## Relevant files

Per system directory:

| File | Role |
|---|---|
| `build.leap`, `leap.log` | system construction and its log |
| `<PDB>_prot.pdb` | cleaned protein used for the build |
| `cand_003_NEW.pdb`, `cand_003_AC.{lib,frcmod}` | docked-pose ligand and its GAFF2 parameters |
| `FAD_NEW.pdb`, `FAD_AC.{lib,frcmod}` | cofactor (MAO systems only) |
| `complex.{prmtop,inpcrd}` | solvated topology and coordinates |
| `complex_dry_check.pdb` | dry complex, for visual verification of the build |
| `production*.dcd` | raw solvated trajectory (primary data) |
| `production*.log` | per-frame state: step, time, energy, temperature, volume, speed |
| `production_full_stripped.nc` + `stripped.complex.prmtop` | **analysis trajectory** — all reported numbers come from here |
| `final_state.xml`, `final.pdb` | restart state and final frame |
| `lig_rmsd.dat`, `prot_rmsd.dat`, `fad_rmsd.dat` | RMSD traces |
| `rmsf_byres.dat` | per-residue CA fluctuation (MAO only) — identifies the flailing C-terminal tail |
| `prot_core_rmsd.dat`, `lig_corefit_rmsd.dat`, `fad_corefit_rmsd.dat` | the same RMSDs refit on the ordered core `:1-496` (MAO only) |
| `mmgbsa_results*.dat` | MM-GBSA output (see step 7) |
| `gb_{complex,receptor,ligand}.prmtop` | GB topologies from `ante-MMPBSA.py` |

Shared:

- [`params/docked/cand_003.acpype/`](params/) — **the ligand parameters actually used**
- `params/cand_003.acpype/` — free-conformer parameters, superseded, kept for provenance
- `params/fad/FAD.acpype/` — FAD parameters, shared by both MAO systems
- `params/cand_003_docked_pose.sdf` — the docked pose the parameters were built from

`system/` additionally holds `production2.dcd` and `mmgbsa_results.dat`: this run was first
analysed at 1 ns (−34.55) and then extended to 10.1 ns. `mmgbsa_results_full.dat` (−31.36) is
the figure to quote; the 1 ns file is kept to show the estimate converged.
