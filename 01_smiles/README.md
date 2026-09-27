# Step 1 — Molecule Sets (SMILES)

## What we did

Defined the exact molecule sets the pipeline operates on: the candidate shortlist, the reference
inhibitor panel, and the RL-generated sets.

| File | Contents |
|---|---|
| [`candidates_56.csv`](candidates_56.csv) | **the 56-candidate shortlist** — SMILES, source campaign, full descriptor set, gate outcomes, RL score and its per-component breakdown |
| [`references.csv`](references.csv) | 16 reference inhibitors with measured activity and provenance |
| [`generation_top20.csv`](generation_top20.csv) | top 20 from the v3 RL campaign |
| [`generation_beats_cand003_all49.csv`](generation_beats_cand003_all49.csv) | all 49 generated molecules that out-scored the lead |

## Why we did it

**Candidates** need to be a frozen, identified set. Every downstream directory keys off the `id`
column (`cand_001` … `cand_056`), so docking runs, MD systems and analysis tables can be joined
without ambiguity. The descriptor and gate columns travel *with* the SMILES so a candidate's
eligibility never has to be re-derived.

**References** exist to anchor the scoring function against measured reality. Without compounds of
known activity there is no way to tell a meaningful docking score from a plausible-looking number —
they are what makes [step 5](../05_validation/)'s calibration possible.

**Generated sets** are kept separate from `candidates_56` because "beats the lead on the RL score"
and "passes the gates and docks well" are different claims. 49 molecules out-scored `cand_003` on
the RL objective; that is a statement about the reward function, not about binding.

## Reference

`references.csv` carries one row per compound with `measured_value`, `measured_unit`, `source`,
`mechanism`, `native_of_pdb` and a `use_in_correlation` flag.

**Every row is now traced to a source (2026-09-27).** No `source` field contains an unresolved
`CONFIRM` marker. Each value's provenance chain, assay format and species are recorded in the file
itself. Chasing them down changed three things materially and they are listed below because each
one is a number a reader might otherwise quote.

### Three corrections the sourcing exercise produced

**1. The 9IV calibration pair was mixing assays.** TTBK1 read 430 nM, which was the midpoint of
BindingDB's 330–530 nM range — a range that **aggregates two assays from different papers**. The
TTBK2 arm read 490 nM from a single paper. A cross-paper midpoint on one arm against a
single-paper value on the other breaks the matched pair, which is the same arm-asymmetry error
this project keeps catching elsewhere.

Both values exist in one paper and one assay format: Nozal et al. 2022 report **TTBK1 330 nM** and
**TTBK2 490 nM** for compound 42 / VNG2.73 (PDB ligand 9IV), inhibition of recombinant human
enzyme with the RICDLHDDEEDEAMSITA substrate, curated as `ChEMBL5200069` whose SMILES matches the
row exactly. Using the arm-symmetric pair moves the experimental ΔΔG from −0.077 to
**−0.234 kcal/mol**, and therefore the protocol's systematic bias from 0.861 to
**1.018 kcal/mol** — see [step 5](../05_validation/).

**2. DTQ's 240 nM is a *K*d, not an IC50.** Xue et al. 2013 measured binding by surface plasmon
resonance. RCSB 4BTK reports both a **Kd of 240 nM** and an **IC50 of 4610 nM** for the same
ligand; the row had the Kd labelled as an IC50, understating the enzymatic potency by **19-fold**.
It is excluded from the correlation, so no reported number moved, but 240 nM must not be quoted as
an IC50.

**3. Selegiline and rasagiline are still not primary-sourced, and cannot easily be.** Their 7.0 and
4.4 nM values trace to Cavalli et al. 2008, which PubMed types as a **Review** — so the values are
compiled there, not measured there. Both are irreversible covalent inhibitors whose IC50 is
preincubation-time dependent and not an equilibrium constant, which is why ChEMBL lists them across
2.76–52 nM and 4.0–46 nM respectively. Both are already excluded from the correlation. If either
number is quoted in the manuscript, chase the review to its own source first.

### Species audit of the correlation set

The docking is against **human** structures. Two of the seven correlation points are not:

| Compound | Target | Value | Assay species |
|---|---|---|---|
| safinamide | MAO-B | 7.67 nM | recombinant human |
| lazabemide | MAO-B | 0.063 µM | human |
| quercetin | MAO-A | 1.52 µM | recombinant human |
| 9IV | TTBK1 / TTBK2 | 330 / 490 nM | recombinant human |
| **kaempferol** | MAO-A | 700 nM | **rat brain** |
| **isatin** | MAO-B | ~3 µM | **rat brain** |

**The two rat points are the two flavonoids** — the compounds chemically closest to this
project's lead, and therefore the anchors a reader would weigh most. Either drop them and state
that the correlation rests on four human points, or keep them and label the correlation
mixed-species. Do not present it as a human-target correlation without saying so. The isatin value
is additionally reported only as "IC50 approximately 3 µM".

### Two standing cautions

- **Safinamide is not SAG.** The true safinamide (a secondary amine) is a different molecule from
  the PDB ligand `SAG` in 2V5Z. Treating them as the same compound conflates the reference with the
  crystallographic ligand.
- **A value is only meaningful with its assay attached.** Safinamide illustrates it: 7.67 nM here,
  but Binda et al. 2007 — the paper behind the 2V5Z structure this project docks into — reports
  *K*i 0.1–0.5 µM for the same compound, 13–65× weaker. Neither is wrong; they are different
  measurements.

## Inputs and parameters

Descriptors (MW, WLogP, TPSA, HBD, HBA, RotB, Fsp3, Rings, AromRings) are RDKit-computed. Gate
columns (`Lipinski_pass`, `GIA_pass`, `BBB_pass`, `alert_free`, `structural_alerts`) come from the
cascade in [step 8](../08_analysis/). `NLL`, `Score` and the `*_raw` columns are REINVENT4 outputs:
the reward components the generator actually optimised (`TTBK1_dock`, `MAOB_dock`, `SA_score`,
`BBB_signed_dist`).

Prepared 3D structures are produced from these SMILES in [step 2](../02_ligands/).

## Analysis performed

- Filtering cascade over the pool → [step 8](../08_analysis/)
- Docking of `candidates_56` and `references` against all receptors → [step 4](../04_docking/)
- Scoring-function calibration against `references.csv` measured values → [step 5](../05_validation/)

## Final result

A frozen 56-candidate shortlist, every member passing Lipinski, GI absorption, BBB and
structural-alert gates, with the lead being:

```
cand_003   Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F
```

A convergence worth noting: the top candidates across independent campaigns carry **CF3 or CHF2**
substituents. That pattern emerged from the generator rather than being designed in.

## Relevant files

| Path | Role |
|---|---|
| [`candidates_56.csv`](candidates_56.csv) | **the shortlist** — the canonical id → SMILES mapping used everywhere downstream |
| [`references.csv`](references.csv) | reference panel with measured activity; **read the caution columns** |
| `generation_top20.csv`, `generation_beats_cand003_all49.csv` | v3 RL output sets |
| [`prior_phase/`](prior_phase/) | phase 1 ligand set, protonation states, and run summary |
