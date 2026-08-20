# Complete Workflow — TTBK1/MAO-B Manuscript Rescue

**Two people · 6 weeks · RTX 4050 (6 GB) + RTX 2050 (4 GB)**

This document replaces the earlier plan. It is self-contained: conventions, scripts, day-by-day steps for both people, sync points, and go/no-go gates.

- **A = Aritra** (Machine A, RTX 4050) — receptor science, MD, MM-GBSA, interpretation, Results/Discussion
- **B = collaborator** (Machine B, RTX 2050) — ligand prep, validation docking, benchmarks, ADMET, references, figures

**Design principle throughout:** *rebuild everything that is cheap and deterministic; keep everything that is expensive and stochastic.* Ligands, receptors, boxes, filters and docking get rebuilt from scratch. The REINVENT4 generative output and the flavonoid library are kept as-is.

---

## Part 1 — Shared conventions

Agree on these on Day 1. Every hour spent here saves a day in Week 4.

### 1.1 Folder structure

Both machines use the identical tree. This is what makes handover trivial.

```
project/
├── scripts/            all shared scripts (see Part 6) — SINGLE source of truth
├── 00_library/         flavonoid library CSV (recovered, unchanged)
├── 01_smiles/          every SMILES set as CSV: id, smiles, source
│   ├── candidates_56.csv
│   ├── references.csv          ← B builds this
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
└── LOGBOOK.md          ← see 1.3
```

### 1.2 Naming rules

- Receptor folders: **PDB ID in caps**, nothing else
- Ligand files: `<setname>_<id>.sdf` — e.g. `cand_012.sdf`, `ref_safinamide.sdf`, `native_VP7.sdf`
- Docking runs: `<RECEPTOR>_<ligandset>_seed<NN>/`
- MD runs: `<RECEPTOR>_<ligand>_rep<N>/`

No spaces in any filename, ever.

### 1.3 The logbook

`LOGBOOK.md`, one entry per work session, both people appending:

```
## 2026-08-22 — B
Prepared 7Q8V and 7Q8Y with scripts/prep_receptor.sh (commit a3f9c1).
Redock 9IV → 7Q8V: 1.12 Å (pass). → 7Q8Y: 4.87 Å (FAIL).
Retried 7Q8Y with dimorphite protonation → 1.44 Å (pass).
NOTE: the original 5.29 Å failure was a protonation problem, not a scoring problem.
Files: 05_validation/calibration_9IV/
```

Two reasons this is not bureaucracy. First, your Methods section gets written from it directly. Second, when a number looks wrong in Week 5, this is how you find out which run produced it.

### 1.4 Git, from hour one

```bash
cd project && git init
printf '*.xtc\n*.trr\n*.dcd\n*.prmtop\n*.rst7\n06_md/*/traj*\n' > .gitignore
git add -A && git commit -m "baseline: recovered artifacts + conventions"
```

Commit after every completed step. Trajectories are gitignored (too large) — back those up to an external drive separately.

### 1.5 The rule that protects your central claim

> **Any two receptors whose scores will be compared must be prepared by the same script, in the same run, with only the filename changed.**

Your headline selectivity result is a *difference* between two receptors. If TTBK1 and TTBK2 were prepared with even slightly different protonation or cleanup, that difference is an artifact of your preparation, not chemistry. This is the most common silent error in comparative docking papers. Never hand-prepare a receptor.

### 1.6 Sync protocol

- **Weekly call, 30 minutes, same day each week.** B reports numbers, A decides what happens next.
- **Handover format:** B never hands over raw output. Every task ends with a CSV plus one paragraph in `LOGBOOK.md` saying what the numbers mean.
- **Escalate immediately, don't wait for the weekly call, if:** a redocking validation fails, a script errors in a way you don't understand, or a number looks wrong. Losing half a day to a question is cheap; losing a week to a silent error is not.

---

## Part 2 — Phase 0: Setup (Days 1–3, both together)

Do this phase sitting together, in one room or one call. B learns the conventions by watching A use them.

### Day 1 — Recovery and environments

**Both:**

1. A downloads every artifact from the Claude Science project into the tree above. Check the manifest — the critical files are the 56-candidate SMILES list, the REINVENT4 config, and the lead compound's acpype/GAFF2 parameters.
2. `git init`, first commit.
3. Create `scripts/` and copy in every script from Part 6.
4. Both set up environments.

**A's environment** — verify the existing one rather than rebuilding:

```bash
conda activate <existing-env>
python -c "import rdkit, MDAnalysis, prolif, parmed, openmm; print('ok')"
python -m openmm.testInstallation      # MUST show a CUDA platform
vina --version
which MMPBSA.py acpype obabel
pip install spyrmsd meeko dimorphite-dl
```

**B's environment** — new, docking only, no MD stack:

```bash
conda create -n dock python=3.11 -y && conda activate dock
conda install -c conda-forge rdkit openbabel numpy pandas matplotlib scipy -y
pip install meeko spyrmsd dimorphite-dl
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64 && sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina
vina --version
```

### Day 2 — The calibration exercise (both, together)

Before any real work, B reproduces a number A already has. This proves the two environments agree.

**Task:** prepare MAO-B 2V5Z, extract safinamide, dock it back, compute symmetry-corrected RMSD.

```bash
cd project
bash scripts/prep_receptor.sh 2V5Z A SAG "FAD"
python scripts/prep_ligands.py --native 03_receptors/2V5Z/native_SAG.sdf -o 02_ligands
bash scripts/dock.sh 2V5Z native_SAG 11
python scripts/rmsd_check.py \
    --ref 03_receptors/2V5Z/native_SAG.sdf \
    --poses 04_docking/2V5Z_native_SAG_seed11/out.sdf
```

**Expected: ≈1.57 Å.** Both machines must agree to within ~0.3 Å.

If they don't agree, stop and fix it today. Two people producing subtly different numbers for five weeks is the worst possible outcome of this project, and it is entirely preventable on Day 2.

### Day 3 — Structure download and inspection

**A downloads all six receptors:**

```bash
cd 03_receptors
for p in 7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X 7JXY 4BTM; do
  mkdir -p $p && wget -O $p/raw.pdb https://files.rcsb.org/download/$p.pdb
done
```

**A inspects each one and records in `LOGBOOK.md`:**

```bash
for p in 7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X; do
  echo "=== $p ==="
  grep "^HETATM" $p/raw.pdb | cut -c18-20 | sort -u | tr '\n' ' '; echo
  grep -c "REMARK 465" $p/raw.pdb
done
```

Record for each: chain to keep, residue range, HETATM species, and whether any missing residues fall inside or near the binding pocket.

| PDB | Protein | Res. | Native ligand | Role |
|---|---|---|---|---|
| **7JXX** | TTBK1 | 1.56 Å | VP7 | **Primary TTBK1** |
| 4BTK | TTBK1 | 2.00 Å | DTQ (240 nM) | TTBK1 cross-check |
| 7Q8V | TTBK1 | 2.13 Å | 9IV (330–530 nM) | Calibration pair |
| 7Q8Y | TTBK2 | 1.60 Å | 9IV (490 nM) | Calibration pair + anti-target |
| 2V5Z | MAO-B | 1.60 Å | SAG safinamide | Primary MAO-B (unchanged) |
| 2Z5X | MAO-A | 2.20 Å | HRM harmine | Anti-target |

**Two specific things A must check on Day 3:**

1. **Does 2Z5X actually contain FAD?** `grep " FAD " 03_receptors/2Z5X/raw.pdb | head` — the PDB entry lists FAD in chain B. If it is present and complete, the FAD transplant described in your Methods was unnecessary, and removing that custom step deletes a caveat from your Limitations.
2. **Do the 7JXX pocket residues match your 4NFM numbering?** Superpose in PyMOL and build a lookup table. Your manuscript names ILE40, ILE48, GLN89, GLN110, GLY111, ASN113, LEU175 in 4NFM numbering; these must be reported correctly in 7JXX numbering.

**GATE 0:** both environments reproduce 1.57 Å; all six structures downloaded and inspected; conventions agreed. → proceed.

---

## Part 3 — Phase 1: Rebuild the foundation (Week 1)

Everything from here is built fresh. You are not patching the old pipeline; you are rebuilding its cheap layers so you own and can describe every step.

### Track A — Week 1

**A1. Prepare all six receptors, in one batch, one command per receptor.**

```bash
bash scripts/prep_receptor.sh 7JXX A VP7 ""
bash scripts/prep_receptor.sh 4BTK A DTQ ""
bash scripts/prep_receptor.sh 7Q8V A 9IV ""
bash scripts/prep_receptor.sh 7Q8Y A 9IV ""
bash scripts/prep_receptor.sh 2V5Z A SAG "FAD"
bash scripts/prep_receptor.sh 2Z5X A HRM "FAD"
git add 03_receptors && git commit -m "receptors prepared, single protocol"
```

Same script, same flags, only the arguments differ. That is the rule from 1.5, enforced mechanically.

**A2. Redocking validation on both TTBK1 receptors.**

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

- **7JXX redocks < 2.0 Å** → you now have a validated on-target TTBK1 receptor. Limitations #2 and #3 are deleted from the manuscript. Proceed.
- **7JXX fails but 4BTK passes** → use 4BTK as primary. Proceed, note it.
- **Both fail** → do not proceed. Go to Part 7 troubleshooting. Work the list in order before concluding anything.

**A3. Rebuild the ligand set from scratch.**

```bash
python scripts/prep_ligands.py --csv 01_smiles/candidates_56.csv -o 02_ligands --ph 7.4
```

This regenerates 3D coordinates, protonation states and conformers for all 56 candidates. It takes minutes and it means every downstream number comes from a protocol you can describe exactly.

**A4. Re-run the filtering cascade.** It is a deterministic script over SMILES. Re-run it so the numbers in your Results are reproducible from your own current code, not from a pipeline you no longer control. Confirm you still get 369 Lipinski / 81 BBB / 334 GI / 253 alert-free / 56 passing all. **If any number differs, find out why before proceeding** — a discrepancy here means the filter definitions drifted, and you need to know which version produced the manuscript figures.

**A5. Consensus re-dock all 56 candidates on the new TTBK1** (and re-run on 2V5Z with the rebuilt ligands so both targets use identical ligand preparation).

```bash
for R in 7JXX 2V5Z; do
  for S in 11 22 33; do bash scripts/dock.sh $R candidates_56 $S; done
done
python scripts/collect_results.py 04_docking -o 08_analysis/consensus_new.csv
```

**Expect the ranking to change.** A holo pocket at 1.56 Å is a different shape from an apo pocket at 2.12 Å. If your lead compound drops in rank, that is the result — it means the old ranking was an apo-structure artifact, which is exactly what you set out to test. Report it openly.

### Track B — Week 1

**B1 — Build the reference compound set.** *(2 days)*

Create `01_smiles/references.csv` with columns `id, smiles, target, measured_value, measured_unit, source`. Get every SMILES from PubChem.

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

Then:

```bash
python scripts/prep_ligands.py --csv 01_smiles/references.csv -o 02_ligands --ph 7.4
```

**Checkpoint before docking anything:** open every generated 3D structure in PyMOL and look at it. Aromatic rings flat? No atoms overlapping? Charges sensible (secondary amines protonated at pH 7.4)? Ten minutes of looking catches errors that would otherwise survive to Week 5.

**B2 — Reference benchmark docking.** *(1 day of compute, unattended)*

```bash
for S in 11 22 33; do
  bash scripts/dock.sh 2V5Z references $S
  bash scripts/dock.sh 2Z5X references $S
done
python scripts/collect_results.py 04_docking -o 05_validation/benchmark_mao.csv
```

**Deliverable:** a CSV, plus a scatter plot of consensus Vina score against pIC50 for every compound with a measured value, Pearson r in the corner. This is the first time anyone will be able to tell whether −11.4 kcal/mol is good.

**B3 — Verification.** Safinamide docked into 2V5Z should reproduce the ≈1.57 Å redocking pose and score in the same range as the manuscript's MAO-B numbers. If it doesn't, something in the rebuilt ligand prep changed — flag it to A immediately, because A's candidate re-dock uses the same prep.

**SYNC POINT 1 (end of Week 1):** A reports the GATE 1 outcome and the new candidate ranking. B reports the benchmark plot. Decide together whether the lead compound is still the lead.

---

## Part 4 — Phase 2: Validation and calibration (Week 2)

### Track A — Week 2

**A6. Build four MD systems.**

| System | Receptor | Replicates | Est. GPU time |
|---|---|---|---|
| 1 | TTBK1 7JXX | 3 × 20 ns | ~10 h |
| 2 | MAO-B 2V5Z | 3 × 20 ns | ~18 h |
| 3 | TTBK2 7Q8Y | 2 × 20 ns | ~7 h |
| 4 | MAO-A 2Z5X | 2 × 20 ns | ~12 h |

Total ≈ 47 GPU-hours ≈ 7–8 overnight runs across Weeks 2–3. **Never run two systems at once on 6 GB.**

Reuse the lead's GAFF2/acpype parameters — same molecule, same charges, transferable to all four systems. Regenerate only if the recovered files are corrupt.

For MAO-B and MAO-A, **FAD belongs to the receptor**, not the ligand. It is covalently attached through an 8α-S-cysteinyl bond (Cys397 in MAO-B, Cys406 in MAO-A). State this explicitly in Methods.

**A7. Replicates differ only in the initial velocity seed:**

```python
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)
simulation.context.setVelocitiesToTemperature(310*kelvin, seed)   # 1, 2, 3
```

Same starting structure, same everything else. Three 20 ns replicates beat one 60 ns run, because only replicates give you a real error bar.

**A8. Start the queue.** Systems 1 and 2 first (they replace existing results). Systems 3 and 4 next (they are new science).

### Track B — Week 2

**B4 — The 9IV matched-pair calibration.** *(3 days — the most valuable thing B does)*

**Why it matters, in one sentence:** the paper claims a compound binds TTBK2 better than TTBK1 by 1.7 kcal/mol, but nobody has ever checked whether this protocol can measure such a difference correctly at all.

The test case: 9IV is crystallised in **both** proteins, with nearly identical measured potency (TTBK1 330–530 nM, TTBK2 490 nM). The true ΔΔG is close to zero — roughly 0.0–0.25 kcal/mol.

```bash
python scripts/prep_ligands.py --native 03_receptors/7Q8V/native_9IV.sdf -o 02_ligands
for R in 7Q8V 7Q8Y; do
  for S in 11 22 33; do bash scripts/dock.sh $R native_9IV $S; done
  python scripts/rmsd_check.py --ref 03_receptors/$R/native_9IV.sdf \
      --poses "04_docking/${R}_native_9IV_seed*/out.sdf"
done
python scripts/collect_results.py 04_docking -o 05_validation/calibration_9IV.csv
```

Then compute `margin = consensus_best(7Q8V) − consensus_best(7Q8Y)`.

**Interpreting the result:**

- **Margin near zero** → the protocol is unbiased. Every selectivity number in §3.7 stands as measured.
- **Margin clearly non-zero** → that is your protocol's **systematic bias**. Report it, and subtract it as a correction from all §3.7 margins.

Either outcome is publishable. The only way to get this wrong is to prepare the two receptors differently — which `prep_receptor.sh` prevents.

**B5 — Re-examine the old 7Q8Y redocking failure.** *(1 day, potentially removes a whole caveat)*

The manuscript attributes a 5.29 Å redock to a Vina scoring-function limitation. For a 1.60 Å structure of an ATP-competitive hinge binder, that is unlikely. Test, in this order, changing one thing at a time:

1. Was the original RMSD symmetry-corrected? A naive atom-order RMSD on a molecule with symmetric aromatic rings can inflate a correct pose by several Å. `rmsd_check.py` uses `spyrmsd`, which corrects for this.
2. **Protonation and tautomer of the pyrrolopyrimidine.** This is the most likely culprit. Open Babel routinely assigns the wrong tautomer for aminopyrimidines, and if the hinge N–H is on the wrong nitrogen, the key hydrogen bonds cannot form and the pose is garbage. Compare `obabel -p 7.4` against `dimorphite-dl`.
3. Box size — enlarge by 4 Å per dimension.
4. Structural waters bridging ligand to hinge — try retaining them.

If any of these fixes it, **that is a finding**: the anti-target validation failure was a preparation problem, now corrected, and the anti-target receptors are validated after all. That deletes Limitation #3 entirely.

**B6 — TTBK1 reference docking.**

```bash
for S in 11 22 33; do bash scripts/dock.sh 7JXX references $S; done
```

Fold into the benchmark plot from B2 so it covers all four receptors.

**SYNC POINT 2 (end of Week 2):** B reports the calibration margin and the 7Q8Y verdict. A folds the correction into the selectivity table. **This is where the paper's central claim becomes defensible.**

---

## Part 5 — Phase 3–5: Analysis and writing (Weeks 3–6)

### Week 3

**A** — MD runs continue overnight. As each completes, run MM-GBSA (igb=5, 0.150 M salt, 100 frames, **per replicate not pooled**). Re-run Vinardo on the new top-15 poses from the validated TTBK1.

**B** — ADMET re-screen: all 56 candidates plus kaempferol and quercetin through **ADMETlab 3.0** and **pkCSM** for BBB, P-glycoprotein substrate likelihood, hERG and hepatotoxicity. Then the important extra: run the **666 candidates from the first failed campaign** through the same BBB predictor. If ADMETlab also finds none permeant, your headline "absolute BBB ceiling" claim is confirmed by two independent models instead of one.

**B** — Reference chasing. Retrieve full bibliographic details for **PMID 17473466, 19748554, 23357036**; read each abstract and confirm it genuinely reports isolation of the flavonoid it's cited for. Add Xue 2013 (*ChemMedChem* 8:1846), Nozal 2022 (*J Med Chem* 65:1585), Bashore 2023 (*Sci Rep* 13:6118), Ahamad 2024 (*Pharmaceuticals* 17:952), Jo 2014 (*Nat Med* 20:886), and the two MAO-B astrocyte papers. Build the library in Zotero. Confirm the "Sharma 2020" → **Chaurasiya et al. 2020, *Molecules* 25:5358** correction.

### Week 4

**A** — Analysis. RMSD/RMSF per replicate; ProLIF interaction fingerprints on the new poses; the corrected selectivity table; and **honest error bars**:

- **Between-replicate SD/SEM** across the 2–3 independent runs — this is the real uncertainty, and it will be several kcal/mol, not 0.23
- Block-averaged error within each trajectory (5 blocks × 20 frames) as secondary

Your current `−36.85 ± 0.23` is a SEM over 100 correlated frames of a single trajectory. It is not a real uncertainty and a referee will say so. The 10.5 kcal/mol MAO-B preference will very likely survive an honest error bar — but it must be an honest one.

**A** — Compute the two decisive numbers: `ΔΔG(TTBK1 − TTBK2)` and `ΔΔG(MAO-B − MAO-A)` by MM-GBSA. This is what your own §4.5 identified as the decisive next step, and it converts the selectivity conclusion from docking-only to cross-method.

**B** — Optional: gnina CNN rescoring on the top-15 poses. Two disagreeing scoring functions cannot adjudicate anything; a third built on a different principle (machine-learned, not empirical) can break the tie and let A advance 2–3 candidates instead of 1.

```bash
./gnina -r receptor.pdbqt -l poses.sdf --autobox_ligand ref.sdf \
        --cnn_scoring rescore --cpu 4 -o rescored.sdf
```

**B** — Regenerate every figure from the new data, consistent style, 300 dpi minimum.

### Week 5 — Writing

**A writes** Results and Discussion. **B writes** Methods sub-sections for the work he did (ligand prep, calibration, benchmark, ADMET), which A then edits for consistency.

Manuscript checklist:

- [ ] **Title** — lead with protocol and constraints, not with a lead compound
- [ ] **Abstract** — BBB ceiling → redesign → dual engagement → *calibrated* selectivity ceiling reproducing published experimental behaviour
- [ ] **Intro ¶1–2** — reframe around Alzheimer's: neuronal TTBK1/tau + astrocytic MAO-B (Jo 2014 and the MAO-B astrocyte papers)
- [ ] **Intro ¶3** — fix the overstated claim that Biogen's is the only characterised TTBK1 series. Cite Xue, Nozal, Bashore, Ahamad.
- [ ] **§2.1** — compress species ranking to two sentences; apparatus → supplementary. Reframe *Evolvulus* as the source of the scaffold hypothesis, not of compounds. Be explicit that the lead is a synthetic trifluoromethylated derivative, not a natural product.
- [ ] **§2.2** — 7JXX replaces 4NFM; add the passing redock table for both on-targets
- [ ] **§2.6** — new subsection: the 9IV matched-pair calibration
- [ ] **§2.7** — replicates, seeds, block-averaging
- [ ] **§3.1–3.2** — unchanged; add second-BBB-predictor confirmation
- [ ] **§3.3** — updated Vina/Vinardo correlation on a validated receptor (+ gnina if done)
- [ ] **§3.4–3.6** — regenerated on new TTBK1; between-replicate error bars
- [ ] **New results subsection** — reference benchmark and score-vs-pIC50 plot
- [ ] **§3.7** — retitle to *"Anti-target counter-screening reproduces the experimentally documented TTBK1/TTBK2 selectivity ceiling"*. Open with the experimental precedent (Bashore 2023: 88% identity, 96% similarity, identical catalytic residues; Biogen/AZ/BMS compounds all roughly equipotent). Then present your margins as reproducing it. Then the calibration. Then the MM-GBSA anti-target numbers.
- [ ] **§4.2** — a field-wide structural limitation, not a defect unique to this series. **Keep the SCA11 safety point** — it is correct and well made.
- [ ] **§4.4** — Limitations #2 and #3 should be *deleted*, not softened. #1 and #4 stay. #8 resolved.
- [ ] **§4.5** — anti-target MM-GBSA moves from "next step" into Results; new next step is the selectivity-aware generative campaign
- [ ] **Data availability** — deposit on Zenodo, cite the DOI

### Week 6 — Buffer and submission

Internal read-through by both. Plagiarism check. Journal formatting. Cover letter that leads with the counter-screen and the two negative results. Submit.

**Target journals:** *Journal of Biomolecular Structure and Dynamics* · *Molecular Diversity* · *Computational Biology and Chemistry* · *Journal of Molecular Graphics and Modelling* · *Scientific Reports* · *RSC Advances*.

---

## Part 6 — Scripts

Put all of these in `scripts/`, commit them, and **never edit them mid-project without a commit and a logbook entry.** Meeko's CLI has changed across versions — check `mk_prepare_receptor.py --help` on your install and adjust the flags once, at the top, for both machines.

### `scripts/prep_receptor.sh`

```bash
#!/usr/bin/env bash
# Usage: prep_receptor.sh PDBID CHAIN LIGRESNAME "COFACTORS"
# Example: prep_receptor.sh 2V5Z A SAG "FAD"
set -euo pipefail
PDB=$1; CHAIN=$2; LIG=$3; COFACTORS=${4:-""}
D="03_receptors/$PDB"; mkdir -p "$D"
[ -f "$D/raw.pdb" ] || wget -qO "$D/raw.pdb" "https://files.rcsb.org/download/$PDB.pdb"

# 1. native ligand out, as its own file (reference pose for redocking validation)
grep "^HETATM" "$D/raw.pdb" | awk -v l="$LIG" -v c="$CHAIN" \
    'substr($0,18,3)==l && substr($0,22,1)==c' > "$D/native_$LIG.pdb"
obabel "$D/native_$LIG.pdb" -O "$D/native_$LIG.sdf" -h 2>/dev/null

# 2. receptor: chosen chain, protein + declared cofactors only. No waters, no cryoprotectant.
KEEP="$COFACTORS"
awk -v c="$CHAIN" -v keep="$KEEP" '
  /^ATOM/   && substr($0,22,1)==c {print}
  /^HETATM/ && substr($0,22,1)==c {
      r=substr($0,18,3); gsub(/ /,"",r)
      if (keep!="" && index(keep,r)>0) print
  }
  /^TER|^END/ {print}
' "$D/raw.pdb" > "$D/clean_noH.pdb"

# 3. protonate at pH 7.4  (swap for pdb2pqr30 if you prefer PROPKA pKa assignment)
obabel "$D/clean_noH.pdb" -O "$D/clean.pdb" -h -p 7.4 2>/dev/null

# 4. pdbqt
mk_prepare_receptor.py --read_pdb "$D/clean.pdb" -o "$D/receptor" -p \
  || prepare_receptor -r "$D/clean.pdb" -o "$D/receptor.pdbqt"

# 5. grid box from the native ligand: centroid + bounding box + 8 A padding
python - "$D" "$LIG" <<'PY'
import sys, json, numpy as np
d, lig = sys.argv[1], sys.argv[2]
c = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
              for l in open(f"{d}/native_{lig}.pdb") if l.startswith("HETATM")])
box = {"center": [round(float(x),2) for x in c.mean(0)],
       "size":   [round(float(max(s,18.0)),1) for s in (c.max(0)-c.min(0)+8.0)]}
json.dump(box, open(f"{d}/box.json","w"), indent=2); print(d, box)
PY
echo "[ok] $PDB prepared"
```

### `scripts/prep_ligands.py`

```python
#!/usr/bin/env python
"""Build 3D ligands from SMILES (or from an extracted native ligand).
   --csv  : CSV with columns id,smiles
   --native : an .sdf extracted from a crystal structure
Outputs <out>/sdf/<id>.sdf and <out>/pdbqt/<id>.pdbqt
"""
import argparse, os, subprocess, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

p = argparse.ArgumentParser()
p.add_argument("--csv"); p.add_argument("--native")
p.add_argument("-o", required=True); p.add_argument("--ph", type=float, default=7.4)
p.add_argument("--nconf", type=int, default=20)
a = p.parse_args()

os.makedirs(f"{a.o}/sdf", exist_ok=True); os.makedirs(f"{a.o}/pdbqt", exist_ok=True)

if a.csv:
    df = pd.read_csv(a.csv)
    items = list(zip(df["id"].astype(str), df["smiles"]))
else:
    m = Chem.MolFromMolFile(a.native)
    items = [(os.path.basename(a.native).replace(".sdf", ""), Chem.MolToSmiles(m))]

for lid, smi in items:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print(f"[FAIL parse] {lid} {smi}", file=sys.stderr); continue
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
    ids = AllChem.EmbedMultipleConfs(m, numConfs=a.nconf, params=ps)
    if not len(ids):
        print(f"[FAIL embed] {lid}", file=sys.stderr); continue
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = min(range(len(res)), key=lambda i: res[i][1])       # lowest MMFF94 energy
    sdf = f"{a.o}/sdf/{lid}.sdf"
    Chem.SDWriter(sdf).write(m, confId=ids[best])
    # protonate at target pH, then convert
    subprocess.run(["obabel", sdf, "-O", sdf, "-p", str(a.ph)],
                   check=True, capture_output=True)
    subprocess.run(["mk_prepare_ligand.py", "-i", sdf,
                    "-o", f"{a.o}/pdbqt/{lid}.pdbqt"], check=True)
    print(f"[ok] {lid}")
```

> **Protonation caveat, and it matters here:** `obabel -p` is fast but gets aminopyrimidine and other tautomers wrong fairly often — which is a leading candidate explanation for the 7Q8Y redocking failure. For the calibration ligand 9IV specifically, also generate states with `dimorphite-dl` and dock both. Whichever reproduces the crystal pose is the one to use, and say so in Methods.

### `scripts/dock.sh`

```bash
#!/usr/bin/env bash
# Usage: dock.sh RECEPTOR LIGANDSET SEED
set -euo pipefail
R=$1; SET=$2; SEED=$3
BOX="03_receptors/$R/box.json"
CX=$(python -c "import json;print(json.load(open('$BOX'))['center'][0])")
CY=$(python -c "import json;print(json.load(open('$BOX'))['center'][1])")
CZ=$(python -c "import json;print(json.load(open('$BOX'))['center'][2])")
SX=$(python -c "import json;print(json.load(open('$BOX'))['size'][0])")
SY=$(python -c "import json;print(json.load(open('$BOX'))['size'][1])")
SZ=$(python -c "import json;print(json.load(open('$BOX'))['size'][2])")
OUT="04_docking/${R}_${SET}_seed${SEED}"; mkdir -p "$OUT"

if [ -f "02_ligands/pdbqt/${SET}.pdbqt" ]; then LIGS="02_ligands/pdbqt/${SET}.pdbqt"
else LIGS=$(ls 02_ligands/pdbqt/*.pdbqt); fi

for L in $LIGS; do
  N=$(basename "$L" .pdbqt)
  vina --receptor "03_receptors/$R/receptor.pdbqt" --ligand "$L" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SX --size_y $SY --size_z $SZ \
       --exhaustiveness 32 --num_modes 9 --seed $SEED \
       --out "$OUT/${N}_out.pdbqt" > "$OUT/${N}.log" 2>&1
  obabel "$OUT/${N}_out.pdbqt" -O "$OUT/${N}_out.sdf" 2>/dev/null
done
echo "[ok] $R $SET seed$SEED"
```

### `scripts/rmsd_check.py`

```python
#!/usr/bin/env python
"""Symmetry-corrected RMSD of docked poses vs a crystal reference.
   Naive atom-order RMSD inflates correct poses on symmetric rings — always use this."""
import argparse, glob
from spyrmsd import io, rmsd

p = argparse.ArgumentParser()
p.add_argument("--ref", required=True); p.add_argument("--poses", required=True)
a = p.parse_args()

ref = io.loadmol(a.ref); ref.strip()
best = (None, 1e9)
for f in sorted(glob.glob(a.poses)):
    for i, m in enumerate(io.loadallmols(f)):
        m.strip()
        r = rmsd.symmrmsd(ref.coordinates, m.coordinates,
                          ref.atomicnums, m.atomicnums,
                          ref.adjacency_matrix, m.adjacency_matrix)
        print(f"{f}  pose {i+1}  RMSD {r:.2f}")
        if r < best[1]: best = (f"{f}:{i+1}", r)
print(f"\nBEST {best[0]}  {best[1]:.2f} A  ->  {'PASS' if best[1] < 2.0 else 'FAIL'}")
```

### `scripts/collect_results.py`

```python
#!/usr/bin/env python
"""Parse Vina logs into a consensus table: best, per-seed mean, inter-seed SD."""
import argparse, glob, os, re
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("root"); p.add_argument("-o", required=True)
a = p.parse_args()

rows = []
for log in glob.glob(f"{a.root}/*/*.log"):
    run = os.path.basename(os.path.dirname(log))
    m = re.match(r"(.+?)_(.+)_seed(\d+)$", run)
    if not m: continue
    rec, lset, seed = m.groups()
    scores = [float(x) for x in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)",
                                           open(log).read(), re.M)]
    if scores:
        rows.append({"receptor": rec, "ligandset": lset, "seed": int(seed),
                     "ligand": os.path.basename(log)[:-4], "best": min(scores)})

df = pd.DataFrame(rows)
out = (df.groupby(["receptor", "ligand"])["best"]
         .agg(best_overall="min", consensus="mean", seed_sd="std", n_seeds="count")
         .reset_index())
out.to_csv(a.o, index=False)
print(out.to_string(index=False))
```

---

## Part 7 — When things fail

**7JXX redock fails (> 2.0 Å).** Change one thing at a time, in this order: (a) confirm symmetry-corrected RMSD was used; (b) enlarge box by 4 Å per dimension; (c) check ligand protonation and tautomer — try `dimorphite-dl` instead of `obabel -p`; (d) retain any structural water bridging ligand to hinge; (e) switch primary to 4BTK. Apply this same list to the old 7Q8Y failure — it is very likely (c).

**The lead compound drops in rank on the validated receptor.** Good, and expected. Advance whatever is now top-ranked and robust across scoring functions. "The validated receptor changed the ranking" is a genuine methodological finding that strengthens your argument about receptor quality — write it up as such.

**MM-GBSA on TTBK2 shows the lead prefers TTBK1 after all.** Then the docking margins were a docking artifact, and cross-method disagreement is the finding. Report both. Do not suppress either.

**B's numbers disagree with A's.** Stop both tracks. Diff the two `LOGBOOK.md` entries, then diff the script versions (`git log scripts/`). Almost always one machine has a different Meeko or Open Babel version. Pin versions and re-run.

**You fall behind schedule.** Cut in this order: gnina → third MD replicate → 4BTK cross-check → ADMET re-screen. **Never cut:** the 7JXX receptor rebuild, the 9IV calibration, or the reference benchmark. Those three are what make the paper publishable.

---

## Part 8 — Master checklist

**Phase 0 — Setup**
- [ ] Artifacts recovered; manifest checked; git initialised
- [ ] Both environments installed
- [ ] **GATE 0:** both machines reproduce safinamide redock ≈1.57 Å
- [ ] Six structures downloaded and inspected; 2Z5X FAD checked; 7JXX residue map built

**Phase 1 — Rebuild (Week 1)**
- [ ] All six receptors prepared by one script
- [ ] **GATE 1:** 7JXX (or 4BTK) redocks < 2.0 Å
- [ ] All 56 candidate ligands rebuilt from SMILES
- [ ] Filtering cascade re-run; numbers reproduce
- [ ] 56 candidates re-docked on 7JXX + 2V5Z
- [ ] B: reference set built, docked, benchmark plot produced
- [ ] **SYNC 1**

**Phase 2 — Validation (Week 2)**
- [ ] Four MD systems built; replicate queue started
- [ ] B: 9IV calibration complete; margin known
- [ ] B: old 7Q8Y failure diagnosed
- [ ] **SYNC 2** — selectivity claim now calibrated

**Phase 3 — Compute (Week 3)**
- [ ] All MD replicates finished
- [ ] MM-GBSA per replicate, all four systems
- [ ] Vinardo re-run on validated receptor
- [ ] B: ADMET + second BBB predictor + 666-set confirmation
- [ ] B: all references chased and formatted

**Phase 4 — Analysis (Week 4)**
- [ ] RMSD/RMSF per replicate; ProLIF on new poses
- [ ] Between-replicate error bars replace single-trajectory SEM
- [ ] ΔΔG(TTBK1−TTBK2) and ΔΔG(MAO-B−MAO-A) computed
- [ ] All figures regenerated

**Phase 5 — Writing (Week 5)**
- [ ] Full manuscript checklist above completed
- [ ] Limitations #2 and #3 deleted, not softened

**Phase 6 — Submit (Week 6)**
- [ ] Internal read-through, plagiarism check, Zenodo DOI, cover letter, submitted

---

## The one-paragraph version

Rebuild the cheap layers from scratch so you own them — ligands, receptors, boxes, filters, docking — using one script per job so that any two receptors you compare are treated identically. Replace the empty TTBK1 structure with a drug-bound one and validate it by redocking. Calibrate your selectivity margin against a compound crystallised in both paralogs with known, nearly equal potency. Dock real drugs so your scores mean something. Extend MM-GBSA to the anti-targets. Then rewrite the selectivity result as a reproduction of a documented experimental ceiling rather than as a confession. Keep the generative run and the library — those were never the problem.
