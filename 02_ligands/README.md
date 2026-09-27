# Step 2 — Ligand Preparation

## What we did

Converted the SMILES sets from [step 1](../01_smiles/) into docking-ready 3D structures: conformer
generation, force-field minimisation, protonation at physiological pH, and conversion to PDBQT.

```
02_ligands/sdf/     77 minimised 3D structures
02_ligands/pdbqt/   77 docking-ready receptor-format ligands
```

## Why we did it

Docking needs a single, physically reasonable 3D conformer per ligand with correct protonation.
Three choices matter:

**Multiple conformers, then pick the lowest energy.** A single embedded conformer can land in a
strained local minimum, which biases the docking search before it starts. 20 conformers are
generated and the lowest-MMFF94-energy one is kept.

**Protonation at pH 7.4, not the neutral form.** Flavonols carry acidic hydroxyls; the neutral
SMILES form is not necessarily what exists at physiological pH, and charge state changes both the
electrostatics and the hydrogen-bond pattern in the site.

**Explicit conversion rather than letting the docking tool infer.** PDBQT atom typing and
rotatable-bond assignment are done by Meeko, so failures surface at preparation time instead of
appearing as an implausible score later.

## Reference

- **Conformers** — RDKit `EmbedMultipleConfs` (ETKDG)
- **Minimisation** — MMFF94 (`MMFFOptimizeMoleculeConfs`, 2000 iterations), lowest-energy conformer retained
- **Protonation** — Open Babel `-p 7.4`
- **PDBQT conversion** — Meeko `mk_prepare_ligand`

## Inputs and parameters

Driver: [`../scripts/prep_ligands.py`](../scripts/prep_ligands.py)

```
--csv     input SMILES table (01_smiles/*.csv)
--ph      7.4    protonation pH
--nconf   20     conformers generated per molecule
-o        02_ligands
```

Native crystallographic ligands take a different path (`--native`), since their 3D coordinates come
from the structure rather than from embedding — see
[`../scripts/fix_native_bondorders.py`](../scripts/fix_native_bondorders.py) and
[step 3](../03_receptors/).

## Analysis performed

Preparation is verified at the point of use rather than assumed: the docking driver checks for the
SDF record terminator (the last thing written, so a truncated file is detectable) and treats Vina's
`BEST None 1000000000.00` as a **preparation failure**, not a score. See
[step 4](../04_docking/README.md).

## Final result

77 ligands prepared in both SDF and PDBQT form — the 56 candidates plus the reference panel and
native ligands. These are the exact inputs every docking run consumes.

## Relevant files

| Path | Role |
|---|---|
| [`pdbqt/`](pdbqt/) | **docking inputs** (77 files) |
| [`sdf/`](sdf/) | minimised 3D structures, also used for RMSD and interaction profiling |
| [`../scripts/prep_ligands.py`](../scripts/prep_ligands.py) | the preparation driver |
| [`prior_phase/`](prior_phase/) | phase 0–8 prepared ligands (`pdbqt_*.pdbqt`) |

**Note:** `pdbqt/` holds more ligands than any single docking run should use. The docking driver
takes an explicit ligand list rather than globbing this directory — see
[step 4](../04_docking/README.md).
