# Week 5 Progress Report — Dual-Target (TTBK1 / MAO-B) Flavonoid Drug Discovery Project

## 1. Project background, in brief

The project's scientific goal is to computationally identify small molecules — built on
a flavonoid (plant-derived) chemical scaffold — that can bind **two different disease
targets at once**: **TTBK1** (a kinase implicated in Alzheimer's-related tau pathology)
and **MAO-B** (an enzyme implicated in Parkinson's-related neurodegeneration). This
"dual-target" idea has not been proposed before in the published literature, which makes
it a genuine, hypothesis-generating research question rather than a replication of known
work.

Before this week, the following had already been done, producing a full draft manuscript:

- Selected a source plant species (*Evolvulus alsinoides*) using a literature-gap
  scoring method, and built a small library of its known flavonoid compounds.
- Used generative AI (a reinforcement-learning model, REINVENT4) to design new
  drug-like molecules built around the flavonoid scaffold, optimized simultaneously for
  binding both targets, synthetic accessibility, and brain penetration.
- Filtered the generated molecules down to a shortlist of 56 candidates that pass
  standard drug-likeness, brain-penetration, and chemical-safety criteria.
- Docked the candidates against both targets computationally, ran molecular dynamics
  simulations and binding free-energy calculations on the top candidate, and drafted a
  full manuscript reporting the results.

Because this is a two-person continuation project, Week 1's job was to put that existing
work on firmer methodological footing before building further on top of it — validating
the computational protocol itself, and rebuilding the parts of the pipeline that are fast
and deterministic to run again, so every number in the eventual paper can be traced back
to a script we control and can re-run on demand.

## 2. Environment and infrastructure setup

- Set up a GPU-accelerated computing environment (CUDA-enabled) for molecular docking
  and molecular dynamics, verified working on an NVIDIA RTX 4050.
- Set up a second, Linux-based environment (via WSL2) specifically for the AmberTools
  software suite, which is used for ligand parameterization and binding free-energy
  calculations and does not have a native Windows version.
- Established a shared folder structure and a written logbook convention so that every
  computational step, and the reasoning behind it, is traceable and reproducible by
  either team member.

## 3. Structural biology setup

- Downloaded all six relevant protein crystal structures from the Protein Data Bank:
  three TTBK1 structures, one TTBK2 structure (the closest relative of TTBK1, used as an
  anti-target/selectivity check), and the MAO-B and MAO-A structures (MAO-A being the
  anti-target/selectivity check for MAO-B).
- Inspected each structure's resolution, bound ligand, and completeness, and confirmed
  that the MAO-A structure already contains its natural cofactor (FAD) fully and
  correctly, simplifying the modelling setup.
- Cross-checked that the new, higher-resolution TTBK1 structure uses the same amino-acid
  numbering as the one referenced in the existing manuscript text, so results can be
  reported consistently. Confirmed all binding-pocket residue positions match, with one
  residue-name correction identified and confirmed against the raw structural data.

## 4. Receptor preparation (single standardized protocol)

All six protein structures were prepared for docking using **one identical, scripted
procedure** applied to every structure the same way — this matters because the project's
central scientific claim is a *comparison* between targets, and that comparison is only
valid if every structure was cleaned, protonated, and formatted identically rather than
by hand.

## 5. Docking-protocol validation

Before trusting any new binding-score numbers, the docking protocol itself was validated
using a standard method called **redocking**: take a molecule whose real, experimentally
determined binding position is already known (from the crystal structure), remove it
computationally, and ask the software to predict where it binds. If the prediction lands
within about 2 Angstroms of the real position, the protocol is considered trustworthy for
that target.

- **TTBK1**: the new structure predicted the known binder's position to within **0.71
  Angstroms** — a high-confidence pass, and a substantial improvement over the older
  structure used in the original manuscript draft.
- A second TTBK1 structure was validated in parallel as a cross-check, also passing
  (0.75 Angstroms).
- **MAO-B**: validated independently by two different people/pipelines, both converging
  on the same result (0.56 Angstroms), confirming the setup is reproducible regardless of
  who runs it.

## 6. Candidate ligand rebuild

The full list of 56 shortlisted candidate molecules was rebuilt from scratch from their
chemical structures (SMILES strings), regenerating clean 3D geometries for every
compound, so this step of the pipeline is under our own control and independently
re-runnable.

## 7. Filtering-pipeline reconstruction and verification

The computational filters originally used to narrow the candidate pool — drug-likeness
rules (Lipinski's Rule of Five), a predicted blood-brain-barrier permeability model, a
predicted gastrointestinal-absorption model, and a chemical safety-alert screen (flagging
reactive or toxicologically risky substructures) — were reconstructed from their
underlying published scientific definitions and re-implemented as a standalone script.

This reconstructed filter was tested against the original filtering results and
reproduces them **exactly**, molecule by molecule, confirming that the filtering logic
behind the 56-candidate shortlist is scientifically sound, well-defined, and now fully
reproducible under our own code rather than depending on a black box.

## 8. Consensus re-docking of all 56 candidates

All 56 candidates were docked against both primary targets (the new TTBK1 structure and
MAO-B), with each docking run repeated three times using different random starting
conditions per compound — 336 individual docking calculations in total — to produce
robust, averaged binding-score estimates rather than single, potentially noisy numbers.
This also gives an internal measure of how consistent each result is.

## 9. Reference/benchmark compound set

A curated set of well-characterized reference compounds was built — approved drugs with
precisely known experimental potency (e.g. safinamide, selegiline, rasagiline), the
proteins' own native co-crystallized ligands, and several natural and endogenous
inhibitor compounds — to serve as a calibration ruler for the computational docking
scores.

- Verified the chemical correctness of every compound in this set (structure, molecular
  formula, molecular weight, and stereochemistry all independently confirmed).
- Ran benchmark docking of the full reference set against both MAO-B and MAO-A.
- Cross-validated the redocking accuracy for the MAO-B reference ligand using a second,
  independently-prepared copy of the molecule, obtaining the identical result as the
  original validation.

## 10. Literature verification of reference potency values

For several reference compounds whose potency values were not yet confirmed, primary
scientific literature was searched (via PubMed) to source verifiable, citable
experimental values, each traceable to a specific published paper (with PMID and DOI
recorded). This let us build a proper score-versus-potency comparison for the
computational docking method, with an honest, literature-grounded assessment of how well
the raw docking score currently tracks real experimental potency — an important
calibration result that will directly inform how the docking scores are interpreted and
reported in the manuscript going forward.

## 11. Where things stand

By the end of Week 1:

- Both computational protocols (docking on TTBK1, docking on MAO-B) are validated against
  experimental crystal structures with sub-1-Angstrom accuracy.
- The candidate filtering pipeline is fully reconstructed, independently verified, and
  reproducible.
- All 56 candidates have fresh, validated, multi-replicate binding scores on both primary
  targets.
- A properly calibrated reference/benchmark dataset is in place, grounded in verified
  primary literature.
- The team is positioned to move into Week 2: molecular dynamics simulations and
  selectivity calibration between each target and its closest biological relative.
