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
`mechanism`, `native_of_pdb` and a `use_in_correlation` flag. MAO IC50 values for kaempferol,
quercetin, lazabemide and isatin were sourced via PubMed (see [`../LOGBOOK.md`](../LOGBOOK.md)).

Two cautions are recorded in the file itself and must be respected:

- **Safinamide is not SAG.** The true safinamide (a secondary amine) is a different molecule from
  the PDB ligand `SAG` in 2V5Z. Treating them as the same compound conflates the reference with the
  crystallographic ligand.
- Several `source` fields are marked *"value as cited in project manuscript Introduction — CONFIRM
  primary source before submission."* Those are **not yet primary-sourced** and are flagged rather
  than quietly trusted.

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
