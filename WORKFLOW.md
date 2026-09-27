# Workflow — TTBK1/MAO-B Dual-Target CADD

**6 weeks · RTX 4050 (6 GB)**

The protocol this project follows: conventions, phase-by-phase steps, go/no-go gates, and failure
modes. Self-contained, so any step can be picked up from this document alone.

> **This is the plan as specified before execution.** [`LOGBOOK.md`](LOGBOOK.md) is the record of
> what actually happened, and it is the primary source wherever the two disagree. Execution
> diverged from this plan in several places worth knowing about up front:
>
> | Planned here | What was actually done |
> |---|---|
> | 4 MD systems, 3 × 20 ns replicates at 310 K | **16 systems**, 10.1 ns each at 300 K — three docked poses per target, plus velocity replicates of the best pose |
> | replicates as the only error estimate | pose scan *and* velocity replicates; the two measure different variance and must not be conflated |
> | FAD stated as covalently attached in Methods | FAD parameterised as a **positionally restrained GAFF2 residue**; the 8α-S-cysteinyl bond is not modelled |
> | GATE 1 on 7JXX alone | all six receptors passed redocking, 0.27–1.45 Å |
>
> The 9IV calibration below was carried out and returned a **0.861 kcal/mol systematic bias toward
> TTBK2**, which is why margins below ~0.9 kcal/mol are not interpretable as selectivity.

**Design principle throughout:** *rebuild everything that is cheap and deterministic; keep
everything that is expensive and stochastic.* Ligands, receptors, boxes, filters and docking get
rebuilt from scratch. The REINVENT4 generative output and the flavonoid library are kept as-is.

---

## Part 1 — Conventions

Settle these on day one. Every hour spent here saves a day in week four.

### 1.1 Folder structure

```
project/
├── scripts/            all driver scripts (Part 6) — SINGLE source of truth
├── 00_library/         flavonoid library CSV (recovered, unchanged)
├── 01_smiles/          every SMILES set as CSV: id, smiles, source
│   ├── candidates_56.csv
│   ├── references.csv
│   └── native_ligands.csv      ← extracted from crystal structures
├── 02_ligands/         PREPARED 3D ligands (rebuilt from scratch)
│   ├── sdf/  pdbqt/
├── 03_receptors/       PREPARED receptors (rebuilt from scratch)
│   ├── 7JXX/ 4BTK/ 2V5Z/ 7Q8V/ 7Q8Y/ 2Z5X/
│   │   ├── raw.pdb  clean.pdb  receptor.pdbqt  box.json  native_ligand.sdf
├── 04_docking/         one subfolder per (receptor × ligandset × seed)
├── 05_validation/      redocking RMSDs, calibration, benchmark plots
├── 06_md/              one subfolder per system per replicate
├── 07_mmgbsa/
├── 08_analysis/        final tables and figures
├── 09_manuscript/
├── 10_results/         consolidated results and interpretation
└── LOGBOOK.md          ← see 1.3
```

The numbered directory names are load-bearing: scripts hardcode them. Add to the structure rather
than rearranging it.

### 1.2 Naming rules

- Receptor folders: **PDB ID in caps**, nothing else
- Ligand files: `<setname>_<id>.sdf` — e.g. `cand_012.sdf`, `ref_safinamide.sdf`, `native_VP7.sdf`
- Docking runs: `<RECEPTOR>_<ligandset>_seed<NN>/`
- MD runs: `<RECEPTOR>_<ligand>_rep<N>/`

No spaces in any filename, ever.

### 1.3 The logbook

[`LOGBOOK.md`](LOGBOOK.md), one entry per work session, appended never rewritten:

```
## 2026-08-22
Prepared 7Q8V and 7Q8Y with scripts/prep_receptor.sh (commit a3f9c1).
Redock 9IV → 7Q8V: 1.12 Å (pass). → 7Q8Y: 4.87 Å (FAIL).
Retried 7Q8Y with dimorphite protonation → 1.44 Å (pass).
NOTE: the original 5.29 Å failure was a protonation problem, not a scoring problem.
Files: 05_validation/calibration_9IV/
```

Two reasons this is not bureaucracy. First, the Methods section gets written from it directly.
Second, when a number looks wrong in week five, this is how you find out which run produced it.

### 1.4 Git, from hour one

```bash
cd project && git init
printf '*.xtc\n*.trr\n*.dcd\n*.prmtop\n*.rst7\n*.nc\n*.frc\n' > .gitignore
git add -A && git commit -m "baseline: recovered artifacts + conventions"
```

Commit after every completed step. Trajectories and MM-GBSA byproducts are gitignored — they are
far too large for a repository, and `reference.frc` alone runs 194–339 MB per system, over GitHub's
100 MB hard per-file limit. Back them up to external storage separately, and note that a fresh
clone therefore cannot reproduce the reported numbers without re-running the MD.

### 1.5 The rule that protects the central claim

> **Any two receptors whose scores will be compared must be prepared by the same script, in the
> same run, with only the filename changed.**

The headline selectivity result is a *difference* between two receptors. If TTBK1 and TTBK2 were
prepared with even slightly different protonation or cleanup, that difference is an artifact of
preparation, not chemistry. This is the most common silent error in comparative docking papers.
Never hand-prepare a receptor.

This generalises into the rule this project kept relearning the hard way: **both arms of any
comparison must be symmetric.** Every false positive recorded in the logbook had the same root
cause — a 0-vs-5 water shell, a pose scan compared against a replicate scan, an on-pose run
compared against an off-pose one. Check arm symmetry before interpreting anything.

---

## Part 2 — Phase 0: Setup (days 1–3)

### Day 1 — Recovery and environment

1. Assemble every recovered artifact into the tree above. The critical files are the 56-candidate
   SMILES list, the REINVENT4 config, and the lead compound's acpype/GAFF2 parameters.
2. `git init`, first commit.
3. Create `scripts/` and populate it (Part 6).
4. Build the environments.

Two environments are needed, because AmberTools has no native Windows build: OpenMM runs natively
on the GPU while parameterisation and MM-GBSA run in WSL2.

```bash
# docking + MD (native, CUDA)
conda activate docking_project
python -c "import rdkit, MDAnalysis, prolif, parmed, openmm; print('ok')"
python -m openmm.testInstallation      # MUST show a CUDA platform
vina --version
pip install spyrmsd meeko dimorphite-dl

# parameterisation + MM-GBSA (WSL2)
conda activate mdgbsa
which MMPBSA.py acpype obabel tleap cpptraj
```

A trap worth recording once: the `python` first on PATH may be MGLTools' Python 2.7, present for
the AutoDock prep scripts. Project scripts must be invoked through `conda run -n docking_project
python`, never a bare `python`. OpenMM likewise must go through `conda run` — invoking the bare
interpreter silently yields only `['Reference', 'CPU', 'OpenCL']` and the run dies with *"There is
no registered Platform called CUDA"* after the restraints are already set up.

### Day 2 — The calibration exercise

Before any real work, reproduce a known number end-to-end. This proves the toolchain agrees with
itself before five weeks of results depend on it.

**Task:** prepare MAO-B 2V5Z, extract safinamide, dock it back, compute symmetry-corrected RMSD.

```bash
bash scripts/prep_receptor.sh 2V5Z A SAG "FAD"
python scripts/prep_ligands.py --native 03_receptors/2V5Z/native_SAG.sdf -o 02_ligands
bash scripts/dock.sh 2V5Z native_SAG 11
python scripts/rmsd_check.py \
    --ref 03_receptors/2V5Z/native_SAG.sdf \
    --poses 04_docking/2V5Z_native_SAG_seed11/out.sdf
```

**Expected: ≈1.57 Å.** If it does not reproduce, stop and fix it today.

### Day 3 — Structure download and inspection

```bash
cd 03_receptors
for p in 7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X; do
  mkdir -p $p && wget -O $p/raw.pdb https://files.rcsb.org/download/$p.pdb
done
for p in 7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X; do
  echo "=== $p ==="
  grep "^HETATM" $p/raw.pdb | cut -c18-20 | sort -u | tr '\n' ' '; echo
  grep -c "REMARK 465" $p/raw.pdb
done
```

Record for each: chain to keep, residue range, HETATM species, and whether any missing residues
fall inside or near the binding pocket.

| PDB | Protein | Res. | Native ligand | Role |
|---|---|---|---|---|
| **7JXX** | TTBK1 | 1.56 Å | VP7 | **Primary TTBK1** |
| 4BTK | TTBK1 | 2.00 Å | DTQ (240 nM) | TTBK1 cross-check |
| 7Q8V | TTBK1 | 2.13 Å | 9IV (330–530 nM) | Calibration pair |
| 7Q8Y | TTBK2 | 1.60 Å | 9IV (490 nM) | Calibration pair + anti-target |
| 2V5Z | MAO-B | 1.60 Å | SAG safinamide | Primary MAO-B |
| 2Z5X | MAO-A | 2.20 Å | HRM harmine | Anti-target |

Two specific checks:

1. **Does 2Z5X actually contain FAD?** `grep " FAD " 03_receptors/2Z5X/raw.pdb | head` — the entry
   lists FAD in chain B. *(This mattered: the prepared MAO-A receptor was initially built without
   its FAD. The flavin forms one wall of the substrate cavity, so docking into a site lacking it is
   docking into the wrong site, and every MAO-A docking before this was caught is invalid.)*
2. **Do the 7JXX pocket residues match the 4NFM numbering used in the manuscript?** Superpose and
   build a lookup table. ILE40, ILE48, GLN89, GLN110, GLY111, ASN113, LEU175 in 4NFM numbering must
   be reported correctly in 7JXX numbering.

**GATE 0:** toolchain reproduces 1.57 Å; all six structures downloaded and inspected;
conventions agreed. → proceed.

---

## Part 3 — Phase 1: Rebuild the foundation (week 1)

Everything from here is built fresh — not patching the old pipeline, but rebuilding its cheap
layers so every step is owned and describable.

**1. Prepare all six receptors, in one batch, one command per receptor.**

```bash
bash scripts/prep_receptor.sh 7JXX A VP7 ""
bash scripts/prep_receptor.sh 4BTK A DTQ ""
bash scripts/prep_receptor.sh 7Q8V A 9IV ""
bash scripts/prep_receptor.sh 7Q8Y A 9IV ""
bash scripts/prep_receptor.sh 2V5Z A SAG "FAD"
bash scripts/prep_receptor.sh 2Z5X A HRM "FAD"
git add 03_receptors && git commit -m "receptors prepared, single protocol"
```

Same script, same flags, only the arguments differ — §1.5 enforced mechanically.

**2. Redocking validation on both TTBK1 receptors.**

```bash
for R in 7JXX 4BTK; do
  L=$(basename 03_receptors/$R/native_*.sdf .sdf)
  python scripts/prep_ligands.py --native 03_receptors/$R/$L.sdf -o 02_ligands
  for S in 11 22 33; do bash scripts/dock.sh $R $L $S; done
  python scripts/rmsd_check.py --ref 03_receptors/$R/$L.sdf \
      --poses "04_docking/${R}_${L}_seed*/out.sdf" | tee 05_validation/${R}_redock.txt
done
```

**GATE 1 — the single most important checkpoint in this project.**

- **7JXX redocks < 2.0 Å** → a validated on-target TTBK1 receptor. Proceed.
- **7JXX fails but 4BTK passes** → use 4BTK as primary. Proceed, note it.
- **Both fail** → do not proceed. Work Part 7 in order before concluding anything.

**3. Rebuild the ligand set from scratch.**

```bash
python scripts/prep_ligands.py --csv 01_smiles/candidates_56.csv -o 02_ligands --ph 7.4
```

Regenerates 3D coordinates, protonation states and conformers for all 56 candidates. Minutes of
compute, and it means every downstream number comes from a describable protocol.

**4. Re-run the filtering cascade.** A deterministic script over SMILES. Re-run it so the Results
numbers are reproducible from current code. Confirm 369 Lipinski / 81 BBB / 334 GI / 253
alert-free / 56 passing all. **If any number differs, find out why before proceeding** — a
discrepancy means the filter definitions drifted.

**5. Build the reference compound set.** `01_smiles/references.csv`, columns
`id, smiles, target, measured_value, measured_unit, source`. Every SMILES from PubChem.

| id | Target | Measured | Note |
|---|---|---|---|
| safinamide | MAO-B | 7.67 nM | also the native ligand of 2V5Z |
| selegiline | MAO-B | 7.0 nM | approved drug |
| rasagiline | MAO-B | 4.4 nM | approved drug |
| kaempferol | MAO-B | — | natural flavonol seed |
| quercetin | MAO-B | — | natural flavonol seed |
| harmine | MAO-A | — | native ligand of 2Z5X |
| clorgyline | MAO-A | — | reference MAO-A inhibitor |
| VP7 | TTBK1 | — | native ligand of 7JXX |
| DTQ | TTBK1 | 240 nM | native ligand of 4BTK |
| 9IV | TTBK1/TTBK2 | 330–530 / 490 nM | native ligand of 7Q8V and 7Q8Y |

**Every measured value must carry a primary source.** A value "as cited in the manuscript
introduction" is not a source, and neither is a Wikipedia chembox. Chase each one to the paper that
reports the assay before submission.

```bash
python scripts/prep_ligands.py --csv 01_smiles/references.csv -o 02_ligands --ph 7.4
```

**Checkpoint before docking anything:** open every generated 3D structure and look at it. Aromatic
rings flat? No overlapping atoms? Secondary amines protonated at pH 7.4? Ten minutes of looking
catches errors that would otherwise survive to week five.

**6. Consensus re-dock all 56 candidates** on the validated TTBK1, and re-run on 2V5Z with the
rebuilt ligands so both targets use identical ligand preparation.

```bash
for R in 7JXX 2V5Z; do
  for S in 11 22 33; do bash scripts/dock.sh $R candidates_56 $S; done
done
python scripts/collect_results.py 04_docking -o 08_analysis/consensus_new.csv
```

**Expect the ranking to change.** A holo pocket at 1.56 Å is a different shape from an apo pocket
at 2.12 Å. If the lead drops in rank, that is the result — the old ranking was an apo-structure
artifact, which is what this exercise set out to test. Report it openly.

**7. Reference benchmark docking.**

```bash
for S in 11 22 33; do
  bash scripts/dock.sh 2V5Z references $S
  bash scripts/dock.sh 2Z5X references $S
  bash scripts/dock.sh 7JXX references $S
done
python scripts/collect_results.py 04_docking -o 05_validation/benchmark_mao.csv
```

**Deliverable:** a CSV plus a scatter of consensus Vina score against pIC50 for every compound with
a measured value, Pearson r in the corner. This is what makes it possible to say whether
−11.4 kcal/mol is good.

**GATE 2:** safinamide into 2V5Z reproduces ≈1.57 Å and scores in the expected range. If not,
something in the rebuilt ligand prep changed, and the candidate re-dock uses the same prep.

---

## Part 4 — Phase 2: Validation and calibration (week 2)

**8. The 9IV matched-pair calibration.** *The most valuable single experiment in the project.*

**Why it matters, in one sentence:** the paper claims a compound binds TTBK2 better than TTBK1 by
1.7 kcal/mol, but nobody has checked whether this protocol can measure such a difference correctly
at all.

The test case: 9IV is crystallised in **both** proteins with nearly identical measured potency
(TTBK1 330–530 nM, TTBK2 490 nM). The true ΔΔG is close to zero — roughly 0.0–0.25 kcal/mol.

```bash
python scripts/prep_ligands.py --native 03_receptors/7Q8V/native_9IV.sdf -o 02_ligands
for R in 7Q8V 7Q8Y; do
  for S in 11 22 33; do bash scripts/dock.sh $R native_9IV $S; done
  python scripts/rmsd_check.py --ref 03_receptors/$R/native_9IV.sdf \
      --poses "04_docking/${R}_native_9IV_seed*/out.sdf"
done
python scripts/collect_results.py 04_docking -o 05_validation/calibration_9IV.csv
```

Then `margin = consensus_best(7Q8V) − consensus_best(7Q8Y)`.

- **Margin near zero** → the protocol is unbiased; selectivity numbers stand as measured.
- **Margin clearly non-zero** → that is the protocol's **systematic bias**. Report it and subtract
  it from every selectivity margin.

*Outcome: 0.861 kcal/mol toward TTBK2, where experiment says there is no preference. Any margin
below ~0.9 kcal/mol is therefore inside the bias and cannot be read as selectivity.*

**9. Re-examine the old 7Q8Y redocking failure.** The manuscript attributes a 5.29 Å redock to a
Vina scoring limitation. For a 1.60 Å structure of an ATP-competitive hinge binder that is
unlikely. Test in this order, changing one thing at a time:

1. Was the original RMSD symmetry-corrected? A naive atom-order RMSD on symmetric aromatic rings
   can inflate a correct pose by several Å. `rmsd_check.py` uses `spyrmsd`, which corrects for it.
2. **Protonation and tautomer of the pyrrolopyrimidine** — the most likely culprit. Open Babel
   routinely assigns the wrong tautomer for aminopyrimidines; with the hinge N–H on the wrong
   nitrogen the key hydrogen bonds cannot form. Compare `obabel -p 7.4` against `dimorphite-dl`.
3. Box size — enlarge by 4 Å per dimension.
4. Structural waters bridging ligand to hinge — try retaining them.

If any of these fixes it, **that is a finding**: the anti-target validation failure was a
preparation problem, now corrected.

**10. Build the MD systems.** Three independent docked poses per target, so pose stability is
tested symmetrically rather than assumed. **Never run two systems at once on 6 GB** — concurrent MD
overwhelmed the machine and forced a restart mid-run.

Reuse the lead's GAFF2/acpype parameters, built **from the docked pose, not a relaxed free
conformer**, so the MD starting geometry is the geometry docking produced.

For MAO-B and MAO-A, **FAD belongs to the receptor**, not the ligand. It is covalently attached
through an 8α-S-cysteinyl bond (Cys397 in MAO-B, Cys406 in MAO-A). That bond is *not* modelled
here: FAD is parameterised as an 84-atom GAFF2 residue and positionally restrained to its
crystallographic coordinates at k = 10 kcal/mol/Å². State this explicitly in Methods, with its
consequence — FAD cannot relax in response to the ligand, so induced fit involving the flavin is
suppressed, and the tight replicate spread this produces is a property of the protocol rather than
of MM-GBSA.

**11. Replicates differ only in the initial velocity seed.** Same starting structure, same
everything else. `setVelocitiesToTemperature` without an explicit seed draws its own, so re-running
the same inputs produces an independent trajectory.

**Replicates are not optional.** They are the only source of a real error bar, and the reported SEM
is not one — see §Part 5.

---

## Part 5 — Phase 3–5: Analysis and writing (weeks 3–6)

### Week 3

- MD runs continue overnight. As each completes, run MM-GBSA (igb=5, 0.150 M salt, every 5th frame,
  **per replicate, never pooled**). Re-run Vinardo on the top-15 poses from the validated TTBK1.
- **Strip solvent before MM-GBSA and do not pass `-sp`.** That flag expects the solvated topology
  and aborts on an atom-count mismatch against an already-stripped trajectory.
- ADMET re-screen: all 56 candidates plus kaempferol and quercetin through **ADMETlab 3.0** and
  **pkCSM** for BBB, P-gp substrate likelihood, hERG and hepatotoxicity. Then run the **666
  candidates from the first failed generative campaign** through the same BBB predictor: if a second
  model also finds none permeant, the "absolute BBB ceiling" claim is confirmed independently.
- Reference chasing. Retrieve full bibliographic details for **PMID 17473466, 19748554, 23357036**
  and confirm each genuinely reports isolation of the flavonoid it is cited for. Add Xue 2013
  (*ChemMedChem* 8:1846), Nozal 2022 (*J Med Chem* 65:1585), Bashore 2023 (*Sci Rep* 13:6118),
  Ahamad 2024 (*Pharmaceuticals* 17:952), Jo 2014 (*Nat Med* 20:886) and the two MAO-B astrocyte
  papers. Confirm the "Sharma 2020" → **Chaurasiya et al. 2020, *Molecules* 25:5358** correction.

### Week 4 — Analysis

RMSD and RMSF per replicate; ProLIF interaction fingerprints; the corrected selectivity table; and
**honest error bars**.

> **The SEM `MMPBSA.py` reports is not the error bar.** It treats ~200 frames sampled 10 ps apart as
> independent draws. Across this project's four replicate pairs the reported SEM stayed in a narrow
> 0.15–0.27 kcal/mol band while the measured spread between two runs differing *only* in velocity
> seed ranged from 0.24 to 6.89 — a 29-fold range.
>
> **How wrong it is tracks pose stability, monotonically:** about 2× for two tightly held MAO poses
> (1.2 Å ligand RMSD), 8× for TTBK1 pose 1 (1.6 Å), **33× for the TTBK2 pair that drifted and
> dissociated** (4.7 Å). A trajectory leaving the site samples structures that were never the
> complex, so its energy swings between seeds; a tightly held pose returns nearly the same number
> twice.
>
> The rule: **the SEM is only as good as the pose is stable, and you cannot tell which case you are
> in without running a replicate.** Quoting the SEM on a ΔΔG is the easiest way to manufacture a
> significant selectivity result from this pipeline.

Two further rules earned the hard way:

- **Never compute MM-GBSA on a trajectory that has left its docked pose.** An energy averaged over
  such a trajectory describes a structure that was never docked. In this project one off-pose
  trajectory produced not merely a noisy number but an **inverted selectivity conclusion**.
- **Pose spread is not replicate spread.** Agreeing from different starting geometries says the
  basin is well defined; it says nothing about how far the energy wanders under resampling. The
  smaller variance cannot bound the larger.

Also compute the two decisive numbers, `ΔΔG(TTBK1 − TTBK2)` and `ΔΔG(MAO-B − MAO-A)`, and regenerate
every figure from the new data at 300 dpi minimum. Optionally rescore the top-15 poses with gnina —
two disagreeing empirical scoring functions cannot adjudicate anything, but a third built on a
different principle can break the tie.

**Absolute MM-GBSA values are not comparable between different proteins.** They carry
protein-specific desolvation and surface terms that do not cancel. A ΔΔG *within* a paralog pair is
quotable; setting a flavoenzyme's number against a kinase's is not.

### Week 5 — Writing

Manuscript checklist:

- [ ] **Title** — lead with protocol and constraints, not with a lead compound
- [ ] **Abstract** — BBB ceiling → redesign → dual engagement → *calibrated* selectivity ceiling
- [ ] **Intro ¶1–2** — reframe around Alzheimer's: neuronal TTBK1/tau + astrocytic MAO-B
- [ ] **Intro ¶3** — fix the overstated claim that Biogen's is the only characterised TTBK1 series
- [ ] **§2.1** — compress species ranking; reframe *Evolvulus* as the source of the scaffold
      hypothesis, not of compounds. Be explicit that the lead is a synthetic trifluoromethylated
      derivative, not a natural product
- [ ] **§2.2** — 7JXX replaces 4NFM; add the passing redock table
- [ ] **§2.6** — new subsection: the 9IV matched-pair calibration
- [ ] **§2.7** — replicates, seeds, and why the reported SEM is not used
- [ ] **§3.1–3.2** — add the second-BBB-predictor confirmation
- [ ] **§3.3** — updated Vina/Vinardo correlation on a validated receptor
- [ ] **§3.4–3.6** — regenerated; between-replicate error bars
- [ ] **New subsection** — reference benchmark and score-vs-pIC50 plot
- [ ] **§3.5** — stability is **pose-dependent**; "both complexes are stable" is true only of
      specific poses, and which pose differs between proteins
- [ ] **§3.6** — the MAO-B preference, with its measured error bar and the restraint caveat
- [ ] **§3.7** — **split by pair.** TTBK1/TTBK2 shows no resolvable discrimination; the MAO pair
      does. They now reach opposite conclusions and must not be written as one finding
- [ ] **§4.2** — a field-wide structural limitation, not a defect unique to this series. **Keep the
      SCA11 safety point** — it is correct and well made
- [ ] **§4.4** — Limitations that the rebuild deleted should be *deleted*, not softened
- [ ] **§4.5** — anti-target MM-GBSA moves from "next step" into Results
- [ ] **Data availability** — deposit on Zenodo, cite the DOI

### Week 6 — Buffer and submission

Read-through, plagiarism check, journal formatting, cover letter leading with the counter-screen
and the negative results. Submit.

**Target journals:** *Journal of Cheminformatics* · *Molecules* · *Frontiers in Chemistry* ·
*Journal of Biomolecular Structure and Dynamics* · *Molecular Diversity* · *Scientific Reports*.

---

## Part 6 — Scripts

All drivers live in [`scripts/`](scripts/) and are the single source of truth. **Never edit one
mid-project without a commit and a logbook entry.**

Earlier revisions of this document inlined full script bodies. They are deliberately not reproduced
here: a second copy drifts from the real one, and a reader following the stale copy gets different
numbers. Read the scripts themselves.

| Script | Step | Job |
|---|---|---|
| [`prep_receptor.sh`](scripts/prep_receptor.sh) | 3 | native ligand out, chain + declared cofactors only, protonate at pH 7.4, pdbqt, grid box from the native ligand + 8 Å padding |
| [`prep_ligands.py`](scripts/prep_ligands.py) | 2 | SMILES or native SDF → 3D, protonated, conformers, pdbqt |
| [`dock.sh`](scripts/dock.sh) | 4 | one Vina run per (receptor × ligandset × seed) |
| [`rmsd_check.py`](scripts/rmsd_check.py) | 5 | symmetry-corrected RMSD via `spyrmsd` |
| [`collect_results.py`](scripts/collect_results.py) | 8 | Vina logs → consensus table; **drops logs older than the receptor they name**, so a stale result cannot be averaged in silently |
| [`make_pose_pdb.py`](scripts/make_pose_pdb.py) | 6 | writes pose *N* onto the acpype template positionally; refuses unless pose 1 round-trips and element sequences agree |
| [`prep_cofactor.py`](scripts/prep_cofactor.py) | 6 | FAD as an 84-atom GAFF2 residue from RCSB ideal chemistry |
| [`run_mmgbsa.sh`](scripts/run_mmgbsa.sh) | 7 | strip → GB topologies → MM-GBSA, inside the WSL2 env |
| [`collect_md_summary.py`](scripts/collect_md_summary.py) | 7 | rebuilds the summary table from primary `.dat` files; `--check` fails if it is stale |
| [`run_mao_pose_queue.sh`](scripts/run_mao_pose_queue.sh) · [`run_mao_replicate_queue.sh`](scripts/run_mao_replicate_queue.sh) | 6 | sequential GPU queues, resumable by skipping finished work |

Meeko's CLI has changed across versions — check `mk_prepare_receptor.py --help` on the install and
adjust the flags once, at the top.

**Tables that get read must be generated, never hand-edited.** A mistyped number in a transcribed
table is invisible and survives review. `collect_md_summary.py --check` exists for exactly this, and
the same discipline applies to prose: a claim in a README is as much a data artifact as a row in a
CSV.

---

## Part 7 — When things fail

**7JXX redock fails (> 2.0 Å).** Change one thing at a time, in order: (a) confirm
symmetry-corrected RMSD; (b) enlarge box by 4 Å per dimension; (c) check ligand protonation and
tautomer — try `dimorphite-dl` instead of `obabel -p`; (d) retain any structural water bridging
ligand to hinge; (e) switch primary to 4BTK. Apply the same list to the 7Q8Y failure; it is very
likely (c).

**The lead compound drops in rank on the validated receptor.** Good, and expected. Advance whatever
is now top-ranked and robust across scoring functions. "The validated receptor changed the ranking"
is a genuine methodological finding.

**MM-GBSA disagrees in sign with docking.** First check whether the trajectory behind the
surprising number stayed on its docked pose. In this project the one sign inversion came entirely
from an off-pose run, and every on-pose trajectory agreed with docking. If both arms are on-pose
and they still disagree, report both and do not suppress either.

**A number looks wrong.** Diff the logbook entries, then the script versions (`git log scripts/`).
Almost always a different Meeko or Open Babel version. Pin versions and re-run.

**The machine restarts mid-run.** `run_md_restrained.py` writes no intermediate checkpoint, so a
partial run cannot be resumed — re-run that system from the start. The queue scripts skip any system
that already has `final_state.xml`, so simply re-running the queue is safe.

**Falling behind schedule.** Cut in this order: gnina → the third MD replicate → the 4BTK
cross-check → the ADMET re-screen. **Never cut:** the receptor rebuild, the 9IV calibration, or the
reference benchmark. Those three are what make the work publishable.

---

## Part 8 — Master checklist

**Phase 0 — Setup**
- [ ] Artifacts recovered; manifest checked; git initialised
- [ ] Environments installed; CUDA platform confirmed
- [ ] **GATE 0:** safinamide redock reproduces ≈1.57 Å
- [ ] Six structures downloaded and inspected; 2Z5X FAD checked; 7JXX residue map built

**Phase 1 — Rebuild (week 1)**
- [ ] All six receptors prepared by one script
- [ ] **GATE 1:** 7JXX (or 4BTK) redocks < 2.0 Å
- [ ] All 56 candidate ligands rebuilt from SMILES
- [ ] Filtering cascade re-run; numbers reproduce
- [ ] 56 candidates re-docked on 7JXX + 2V5Z
- [ ] Reference set built and docked; benchmark plot produced
- [ ] Every reference value traced to a primary source
- [ ] **GATE 2**

**Phase 2 — Validation (week 2)**
- [ ] 9IV calibration complete; systematic bias known
- [ ] Old 7Q8Y failure diagnosed
- [ ] MD systems built — three docked poses per target
- [ ] Replicate queue started, one system at a time

**Phase 3 — Compute (week 3)**
- [ ] All MD runs finished
- [ ] MM-GBSA per run, on-pose trajectories only
- [ ] Vinardo re-run on validated receptor
- [ ] ADMET + second BBB predictor + failed-campaign confirmation
- [ ] All references chased and formatted

**Phase 4 — Analysis (week 4)**
- [ ] RMSD/RMSF per run; ProLIF on new poses
- [ ] Between-replicate spread replaces the single-trajectory SEM everywhere
- [ ] ΔΔG(TTBK1−TTBK2) and ΔΔG(MAO-B−MAO-A) computed, each against its own measured spread
- [ ] All figures regenerated from the collected table, not by hand

**Phase 5 — Writing (week 5)**
- [ ] Manuscript checklist completed
- [ ] §3.7 split by pair
- [ ] Resolved limitations deleted, not softened

**Phase 6 — Submit (week 6)**
- [ ] Read-through, plagiarism check, Zenodo DOI, cover letter, submitted

---

## The one-paragraph version

Rebuild the cheap layers from scratch so you own them — ligands, receptors, boxes, filters,
docking — using one script per job so that any two receptors you compare are treated identically.
Replace the empty TTBK1 structure with a drug-bound one and validate it by redocking. Calibrate the
selectivity margin against a compound crystallised in both paralogs with known, nearly equal
potency. Dock real drugs so the scores mean something. Run three docked poses per target, because
docking rank does not predict which pose survives dynamics, and replicate the ones that hold,
because the reported SEM is not an error bar. Never compute a binding energy on a trajectory that
has left its docked pose. Then write the selectivity result as what it is — a resolvable preference
in one paralog pair and an unresolvable one in the other. Keep the generative run and the library;
those were never the problem.
