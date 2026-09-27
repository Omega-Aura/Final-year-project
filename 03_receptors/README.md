# Step 3 — Receptor Structures and Preparation

## What we did

Selected, characterised and prepared six crystal structures as docking targets, then produced
several deliberate variants of each to serve as controls.

| Directory | Protein | Role |
|---|---|---|
| [`7JXX/`](7JXX/) | TTBK1 | primary target |
| [`7Q8V/`](7Q8V/) | TTBK1 | second TTBK1 structure; holds 9IV for calibration |
| [`7Q8Y/`](7Q8Y/) | TTBK2 | anti-target / selectivity comparison; holds 9IV |
| [`2V5Z/`](2V5Z/) | MAO-B | primary target |
| [`2Z5X/`](2Z5X/) | MAO-A | anti-target / selectivity comparison |
| [`4BTK/`](4BTK/) | — | additional redocking validation target |

Variant directories are controls, not duplicates:

| Suffix | Preparation | Purpose |
|---|---|---|
| *(none)* | crystallographic water shell retained | production runs |
| `dry` | all waters removed | makes a two-protein comparison **water-symmetric** |
| `noFAD` | FAD cofactor removed | tests whether candidate scores depend on the cofactor |
| `brg` | bridging water retained | isolates the effect of one specific water |

[`_cofactors/FAD_ideal.sdf`](_cofactors/) holds the RCSB ideal-chemistry FAD used as the single
source for cofactor parameterisation across both MAO systems.

## Why we did it

**Because structure choice silently determines the answer.** Every methodological failure this
project recorded traces back to a receptor-preparation asymmetry rather than to the science:

- Two kinase structures were compared while retaining **five waters in one site and none in the
  other**. That asymmetry alone produced an apparent selectivity margin. The `dry` variants exist
  so the comparison can be run symmetrically.
- The prepared MAO-A receptor was initially **missing its FAD cofactor** entirely (documented in
  [`receptor_characterization.md`](receptor_characterization.md) §5b). FAD forms one wall of the
  substrate cavity; docking into a site without it is docking into the wrong site.

Characterising each structure *before* docking — resolution, covalent linkages, incomplete side
chains, chain gaps near the site, residue numbering — is what makes these traps findable at all.

## Reference

[`receptor_characterization.md`](receptor_characterization.md) is the primary document for this
step and covers, per structure:

1. Protein identity and the **residue-numbering caveat** — numbering differs between structures,
   so pocket residues cannot be compared by number alone
2. Chain length, and incomplete side chains flagged by `REMARK 470`
3. Resolution and refinement quality
4. Bond inventory — covalent linkages **verified from coordinates, not from header records**,
   coordination/ionic bonds, and the non-covalent interaction types available in each pocket
5. Missing content — absent cofactors, absent inhibitors, chain gaps near the binding site
6. A consolidated list of issues found

Prior-phase target rationale and the literature behind target selection are in
[`prior_phase/`](prior_phase/): `phase0_targets_report.md` (disease biology and target
validation) and `phase0_target_citations.csv` (the citation list). `phase0_structures.csv` records
the structure survey; `grid_boxes.json` the earlier box definitions.

## Inputs and parameters

Prepared by [`../scripts/prep_receptor.sh`](../scripts/prep_receptor.sh). Per receptor directory:

| File | Role |
|---|---|
| `raw.pdb` | the deposited structure, unmodified |
| `clean_noH.pdb` | cleaned heavy-atom structure |
| `receptor.pdbqt` | the docking-ready receptor |
| `box.json` | grid box `center` and `size` (Å) |
| `native_<LIG>.pdb` / `.sdf` | the co-crystallised ligand as deposited |
| `native_<LIG>_ref.sdf` | bond-order-repaired reference for RMSD |

Box sizes are 18 Å cubes centred on the native ligand. Native ligand bond orders are repaired by
[`../scripts/fix_native_bondorders.py`](../scripts/fix_native_bondorders.py) before use as an RMSD
reference, since PDB-derived bond orders are unreliable and a wrong bond order changes the
symmetry graph the RMSD is computed over.

A note carried from the MD step: the **Cys-FAD covalent bond length differs between structures** —
2.31 Å in 2V5Z against 1.65 Å in 2Z5X. These are two independently refined copies of the same
cofactor, so cross-structure geometry comparisons need thresholds calibrated for that, not the
tight thresholds appropriate to a same-source transfer.

## Analysis performed

- Structural characterisation of all six targets → [`receptor_characterization.md`](receptor_characterization.md)
- Native-ligand redocking of every receptor → [step 5](../05_validation/)
- Paired wet/dry docking to quantify the water artifact → [step 4](../04_docking/) and
  [step 8](../08_analysis/)

## Final result

All six structures are docking-validated: native-ligand redocking recovers the crystallographic
pose in every case, from 0.27 Å (2Z5X) to 1.45 Å (7Q8Y). See [step 5](../05_validation/) for the
full table.

Two preparation defects were found and corrected before any production docking: the missing MAO-A
FAD cofactor, and the asymmetric water shell between the kinase pair. Both are documented with
dates in [`../LOGBOOK.md`](../LOGBOOK.md).

## Relevant files

| Path | Role |
|---|---|
| [`receptor_characterization.md`](receptor_characterization.md) | **primary document** — per-structure characterisation and the issues found |
| `<PDB>/receptor.pdbqt` + `<PDB>/box.json` | what docking actually consumes |
| `<PDB>/native_*_ref.sdf` | RMSD references for validation |
| [`_cofactors/FAD_ideal.sdf`](_cofactors/) | single source for FAD parameterisation |
| [`../scripts/prep_receptor.sh`](../scripts/prep_receptor.sh) | preparation driver |
| [`../scripts/prep_cofactor.py`](../scripts/prep_cofactor.py) | FAD extraction and parameterisation |
| [`prior_phase/`](prior_phase/) | target validation report, citations, structure survey, and the earlier receptor set (`4NFM`, `6U0K`, `2V60`) |
