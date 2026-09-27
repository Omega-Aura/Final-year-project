---
title: "Generative design and structure-based evaluation of 7-deoxyflavonol dual-target candidates against TTBK1 and MAO-B, with explicit anti-target counter-screening"
subtitle: "A hypothesis-generating computational study"
date: 2026-08-19
---

# Abstract

**Background.** Tau-tubulin kinase 1 (TTBK1) and monoamine oxidase B (MAO-B) are both
CNS drug-discovery targets in neurodegeneration, but they belong to distinct disease
mechanisms and have never been proposed as a combined dual-target pair in the published
literature. Flavonoids are documented MAO-B inhibitors and promiscuous kinase binders,
making them a plausible starting scaffold for a multi-target-directed ligand (MTDL)
hypothesis.

**Methods.** We assembled a flavonoid library from *Evolvulus alsinoides*, a medicinal
species selected by a pre-registered two-axis rule that multiplies structurally-verified
flavonoid coverage by literature-gap score. Reinforcement-learning molecular generation
(REINVENT4 LibInvent) decorated a 7-deoxyflavonol core under a four-component reward
(TTBK1 docking, MAO-B docking, synthetic accessibility, BOILED-Egg blood-brain-barrier
signed distance). Candidates were filtered on Lipinski, BBB and PAINS/BRENK criteria,
re-docked at high exhaustiveness with multi-seed consensus, cross-validated with an
independent scoring function (Vinardo), counter-screened against the paralog/isoform
anti-targets TTBK2 and MAO-A, and the single most scoring-function-robust candidate was
advanced to 20 ns explicit-solvent molecular dynamics and MM-GBSA binding free energy
calculation on both on-targets.

**Results.** The original flavonol scaffold hit an absolute BBB ceiling: minimum
topological polar surface area across 666 unique generated candidates was 90.9 A^2,
outside the BOILED-Egg BBB region at any R-group decoration, yielding a null shortlist
(0/666). Redesign to a 7-deoxyflavonol core with a BBB-aware reward recovered 56
candidates passing all filters. Cross-scoring-function agreement over the top-15 was weak
(TTBK1 r = 0.56; MAO-B r = -0.28), so ranking was restricted to the one candidate robust
under both functions: a 4-trifluoromethyl/methyl-decorated 7-deoxyflavonol. Over 20 ns of
MD this candidate remained stably bound to both targets (last-quartile ligand RMSD
1.52 A at TTBK1, 0.68 A at MAO-B). MM-GBSA gave favorable binding free energies at
both targets, with MAO-B preferred: -36.8 +/- 0.2 kcal/mol
versus -26.3 +/- 0.5 kcal/mol for TTBK1 (mean +/- SEM,
n = 100 frames), driven mainly by van der Waals contact (-44.9 vs. -32.9
kcal/mol). Anti-target docking, however, showed essentially no isoform discrimination:
0/15 shortlisted candidates favored TTBK1 over TTBK2, and the lead candidate's MAO-B-over-
MAO-A margin was +0.09 kcal/mol.

**Conclusions.** Two independent computational methods concur that this 7-deoxyflavonol
engages both TTBK1 and MAO-B, but the same pipeline shows the series does not discriminate
either target from its closest paralog/isoform. We therefore report a validated dual
on-target engager with an unresolved selectivity liability, not a selective dual-target
lead. This is a hypothesis-generating computational study; no experimental validation was
performed and none of the reported affinities are measured values.

**Keywords:** TTBK1, MAO-B, flavonol, multi-target-directed ligand, generative
reinforcement learning, MM-GBSA, anti-target counter-screening, blood-brain barrier

---

# 1. Introduction

Tau-tubulin kinase 1 (TTBK1, UniProt Q5TCY1) phosphorylates tau at disease-relevant
epitopes including Ser422 and is genetically and biochemically implicated in Alzheimer's
disease and other tauopathies, with a smaller literature linking TTBK1 and its paralog
TTBK2 to TDP-43 proteinopathy relevant to ALS/FTLD-TDP. TTBK1 chemical-probe development
remains early-stage: the only well-characterized inhibitor series with public measured
potencies identified in this work is the Biogen brain-penetrant azaindazole/pyrrolopyridine
series (Halkina et al., *J Med Chem* 2021, 64:6358; PMID 33944571), whose lead reduces tau
pSer422 in vivo.

Monoamine oxidase B (MAO-B, UniProt P27338) is by contrast a clinically validated
Parkinson's disease target. It oxidatively deaminates dopamine and other monoamines,
producing H2O2 and contributing to oxidative stress and dopaminergic neurodegeneration, and
is upregulated in reactive astrocytes in neurodegenerative tissue. Three approved inhibitors
(safinamide, IC50 = 7.67 nM; selegiline, 7.0 nM; rasagiline, 4.4 nM) provide well-measured
reference pharmacology.

**The dual-target hypothesis examined here is novel and unvalidated.** Systematic PubMed
querying returned zero records combining flavonoids with TTBK1, and zero records proposing
a TTBK1 + MAO-B combination in any context. Flavonoids as MAO-B inhibitors, by contrast, is
an established sub-field. TTBK1 and MAO-B are not co-implicated in any single validated
disease mechanism in the retrieved literature — TTBK1 sits in tau/Alzheimer's-spectrum
biology, MAO-B in dopaminergic/Parkinson's-spectrum biology. The dual-target framing of this
study therefore rests on three weaker premises, stated plainly: both are licensed CNS
neurodegeneration targets; flavonoids are documented to be promiscuous across kinase and MAO
pharmacology; and multi-target-directed ligand design combining kinase inhibition with MAO-B
inhibition or antioxidant activity has general precedent in neurodegenerative drug design.
No claim of an established dual-target consensus is made or implied.

Two design constraints shaped the study. First, both targets are CNS targets, so
blood-brain-barrier permeability was treated as a hard requirement rather than a
post-hoc annotation. Second, because both targets have a close paralog or isoform whose
active site is highly similar — TTBK2 for TTBK1, MAO-A for MAO-B — selectivity cannot be
assumed and was tested explicitly by anti-target counter-screening rather than asserted.

# 2. Methods

## 2.1 Species selection and natural-product library

Six Indian medicinal species with traditional neurological indications were screened using
a pre-registered multiplicative ranking rule:

```
flavonoid_coverage = n_high_confidence_flavonoids / max(n_high_confidence_flavonoids)
literature_gap     = 1 - (PubMed CADD-term hits / PubMed species-alone hits)
composite          = flavonoid_coverage x literature_gap
```

The multiplicative form is deliberate: a species with a large literature gap but no
structurally-verified flavonoids scores zero, since a gap without dockable chemistry is not
an opportunity. *Bacopa monnieri* was included as a positive control for a crowded
literature (566 species records, 49 CADD records). *Evolvulus alsinoides* ranked first
(composite 0.933), carrying six structurally-confirmed flavonoid/flavonol compounds
traceable to primary isolation papers (PMIDs 17473466, 19748554, 23357036), a modest total
footprint (60 records), only 4 CADD-related records, and zero records combining the species
with either TTBK1 or MAO-B. Its dominant scaffold class is flavone/flavonol
(3-hydroxyflavone core), which set the scaffold constraint for generative design.
The final library comprised 24 entries; all 24 SMILES parsed in RDKit, with per-compound
confidence tiers (11 HIGH, 6 MEDIUM, 4 LOW, 3 DERIVED) recorded from source provenance.

## 2.2 Receptor preparation and redocking validation

On-targets: MAO-B PDB **2V5Z** (1.60 A, safinamide-bound, native FAD) and TTBK1 PDB
**4NFM** (2.12 A kinase domain). Anti-targets: MAO-A PDB **2Z5X** (2.20 A, harmine-bound)
and TTBK2 PDB **7Q8Y** (1.60 A) and **6U0K** (1.74 A). Ligands were protonated at pH 7.4.

Two structure-level corrections were made against the initial study design and are reported
rather than silently applied. First, the originally specified MAO-B structure 2V60 does not
contain safinamide — its co-crystallized ligand is a coumarin-4-carbaldehyde analog — so
2V5Z was used instead, which is both the genuine safinamide complex and higher resolution.
Second, 4NFM is apo (only cryoprotectant glycerol is bound), so native-ligand redocking
validation is impossible on that structure; the grid box was instead defined from the
resolved catalytic Lys38 and DFG176 landmarks and transferred from 7Q8Y by kinase-domain
superposition.

Redocking validation results are asymmetric and constrain interpretation throughout:

| Receptor | Native ligand | Best RMSD (A) | Pass (< 2.0 A) |
|---|---|---|---|
| MAO-B 2V5Z | safinamide | 1.57 | **yes** |
| TTBK2 7Q8Y | 9IV | 5.29 | **no** |
| TTBK1 4NFM | (apo) | not testable | n/a |
| MAO-A 2Z5X | (FAD transplanted from 2V5Z) | not performed | n/a |
| TTBK2 6U0K | (not performed) | not performed | n/a |

Only MAO-B carries a passing redocking validation. The TTBK2 7Q8Y failure was investigated
(no steric clash identified; re-run at exhaustiveness 32 with 20 poses) and attributed to a
genuine Vina scoring-function limitation for that ligand/pocket rather than a setup error.

## 2.3 Generative design under a BBB-aware reward

Molecular generation used REINVENT4 LibInvent to decorate a fixed scaffold core.
An initial campaign on the native flavonol core (three mandatory phenolic/carbonyl oxygens
at 3-OH, 5-OH, 7-OH) used a three-component reward (TTBK1 docking 0.35 / MAO-B docking 0.35
/ synthetic accessibility 0.15, plus BBB 0.15 introduced in the second campaign), Vina at
exhaustiveness 4 for tractability across ~1,000 candidates per run.

After that campaign failed the BBB gate outright (Section 3.1), the scaffold was redesigned
to a **mono-deoxy 7-deoxyflavonol** — retaining the 5-OH/4-carbonyl intramolecular hydrogen
bond and removing the 7-OH — and the reward was rebuilt with a fourth, direct BOILED-Egg BBB
signed-distance term (weights: TTBK1 0.35, MAO-B 0.35, SA 0.15, BBB 0.15; geometric mean).
An initial scaffold construction that inadvertently retained 7-OH and removed 5-OH was
detected and corrected, and the 32-step run was re-executed from scratch on the corrected
regiochemistry.

## 2.4 Filtering cascade

Applied in fixed order: RDKit sanitization and canonical deduplication with novelty check
against the natural-product seed library; Lipinski Ro5; BOILED-Egg GI-absorption and BBB
ellipse gates (Daina & Zoete 2016, computed from RDKit WLOGP/TPSA); PAINS (A/B/C) and BRENK
reactive-group alerts via the RDKit FilterCatalog. The BBB gate was applied as a literal
mandatory filter, not an advisory annotation.

## 2.5 Consensus re-docking and independent scoring-function cross-check

All shortlisted candidates were re-docked against both on-targets with AutoDock Vina at
exhaustiveness 32, 9 poses, and 3 independent random seeds (11/22/33), each seed using an
independent RDKit conformer embedding rather than a re-seeded search on the same starting
geometry. Reported per ligand: best score across all seeds and poses, consensus score (mean
of the per-seed best pose), and inter-seed standard deviation.

Because both the reward function and the re-dock used the same Vina scoring function, an
independent cross-check was run: the top-15 dual-target candidates were re-docked from
scratch and rescored with **Vinardo** (Quiroga & Villarreal 2016), a distinct functional form
and parameterization.

## 2.6 Anti-target counter-screening

The same consensus protocol was applied to all three anti-target receptors (MAO-A 2Z5X,
TTBK2 6U0K, TTBK2 7Q8Y). Selectivity margin was defined as
`on-target consensus best - anti-target consensus best`, so a negative margin indicates
favorable on-target selectivity.

## 2.7 Molecular dynamics and MM-GBSA

The single most scoring-function-robust candidate (Section 3.3) was advanced to MD against
both on-targets. Systems were built with PDBFixer, ligands parameterized with acpype/GAFF2,
topologies assembled with GROMACS `pdb2gmx`, solvated in explicit water with counterions,
and energy-minimized and propagated with **OpenMM** reading the GROMACS topology and
coordinates directly (the local GROMACS build was OpenCL-only and did not detect the
available CUDA device). Production runs were 20 ns per system on a single NVIDIA RTX 4050
(6 GB), executed sequentially (TTBK1 3.19 h; MAO-B 5.99 h). MAO-B's FAD cofactor was
retained throughout.

Trajectory analysis (MDAnalysis) computed protein and ligand RMSD against the **first
production frame**, per-residue RMSF, and hydrogen-bond occupancy using explicit
element-based donor/acceptor selections.

Binding free energies used AmberTools **MMPBSA.py** in single-trajectory MM-GBSA mode
(igb = 5, salt concentration 0.150 M, 100 evenly-spaced frames). GROMACS topologies were
converted to Amber `prmtop` with ParmEd, with three conversions required for correctness:
periodic box information was cleared for implicit solvent; dihedral types carried over with
GROMACS' zero-periodicity/zero-force-constant convention were reassigned a nonzero
periodicity (energetically inert, since the force constant is zero, but required by the
Amber parser); and Born radii and screening parameters — absent from any GROMACS topology —
were assigned with ParmEd `changeRadii mbondi2` to match igb = 5. For MAO-B, FAD was
assigned to the **receptor**, not scored as ligand, since it is a covalently-linked
biological component of the active site; only the small molecule was scored as ligand.

## 2.8 Reproducibility and limitations of the computational protocol

All reported affinities are computed, not measured. MM-GBSA in single-trajectory mode
neglects conformational entropy and receptor reorganization, and systematically overestimates
the magnitude of binding free energies; its values are used here for **relative** comparison
between the two on-targets and are not interpretable as absolute affinities or convertible
to Kd. Single 20 ns replicates per system, on one candidate, are a screening-level and not a
converged sampling protocol.

# 3. Results

## 3.1 The native flavonol scaffold has an absolute BBB ceiling

The first generative campaign produced 1,280 scored molecules, reducing to 666 unique,
valid, novel candidates. The filtering cascade gave 587 passing Lipinski, 198 passing
GI-absorption, and 217 alert-free — but **zero** passing the BOILED-Egg BBB gate.

The cause is structural, not statistical. Minimum TPSA across all 666 candidates was
**90.9 A^2**, a floor set by the flavonol core itself: the mandatory 3-OH, 4-C=O and 5-OH
oxygens plus the pyranone ring oxygen. The BOILED-Egg BBB region requires roughly
TPSA < 79 A^2 at the relevant WLogP, so no B-ring decoration can bring a native flavonol
into the permeant region. Both natural seed compounds (kaempferol, quercetin) and the
literature TTBK1/2 reference inhibitor also fall outside the ellipse. Under the literal
BBB gate the formal output of this campaign is a null shortlist.

## 3.2 Scaffold deoxygenation converts a null result into 56 candidates

Removing the 7-OH to give a 7-deoxyflavonol core, while retaining the 5-OH/4-carbonyl
intramolecular hydrogen bond, and adding a direct BBB signed-distance reward term produced
1,024 candidates over 32 steps (405 unique, valid, novel). Mean composite reward rose from
0.08 at step 1 to a 0.45-0.53 plateau by steps 20-32.

The BBB pass rate over the run is informative about the underlying trade-off rather than a
clean success curve: steps 1-2, still near the unmodified prior, start at 71.9% and 62.5%;
steps 3-8 collapse to a 0-34.4% range (mean 16.9%) as reinforcement learning optimizes
docking affinity at the expense of permeability; from step ~9 the BBB-aware term pulls the
rate back, averaging 45% over steps 13-32 but still ranging 12.5-65.6% between 32-candidate
batches, with no further downward trend. Overall run pass rate was 389/1024 (38%).

Final cascade: 369 Lipinski-passing, 81 BBB-passing, 334 GI-passing, 253 alert-free, and
**56 candidates passing all filters simultaneously** — the qualitative reversal the redesign
targeted. The 56 span MW 310-420 Da, TPSA 70.7-77.2 A^2 (inside the BBB ellipse), WLogP
3.6-4.8, and 8 unique Murcko scaffolds.

![Composite reward progression across the reinforcement-learning run.]({{artifact:art_575886ad-749c-4c0e-986b-3b8840a9af6a}})

## 3.3 Docking rank order is not robust to the choice of scoring function

Multi-seed consensus re-docking at exhaustiveness 32 reproduced the reward-time ranking
well within the same scoring function (TTBK1 Pearson r = 0.96; MAO-B r = 0.74) with high
inter-seed reproducibility (mean best-pose SD 0.029 kcal/mol at TTBK1, 0.094 at MAO-B).
This establishes search convergence, not scoring-function reliability.

The independent Vinardo cross-check on the top-15 shows weak agreement:

| Target | Vina-consensus vs. Vinardo (Pearson r) |
|---|---|
| TTBK1 | 0.56 |
| MAO-B | **-0.28** |

This is a negative result and is reported as such. All 15 candidates score favorably in
absolute terms under both functions (Vina consensus -7 to -12 kcal/mol; Vinardo -5 to -9.5
kcal/mol, a known scale offset), so binary "does this scaffold engage both sites plausibly"
calls survive. What does not survive is using the fine-grained numeric score to select a
single best compound: the Vina-consensus top pick and the Vinardo top pick are different
positional regioisomers of the same benzyl/methyl-decorated core.

![Vina-consensus versus Vinardo cross-check for the top-15 candidates.]({{artifact:art_1c8a5676-6a52-4bfe-9b11-145f993bc776}})

Ranking was therefore restricted to candidates robust under **both** functions
independently. Only one qualifies — `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`, a 4-trifluoromethyl/methyl-decorated
7-deoxyflavonol, ranked #2 under Vina-consensus alone and #2 under Vinardo alone. The
combined-average #1 candidate is Vinardo #1 but only ~8th of 15 under Vina-consensus; its
weak MAO-B Vina-consensus score (-10.32 kcal/mol, worst of the top-15 on that metric) is
masked by averaging rather than corroborated. This candidate ("the lead") was advanced to MD.

| Candidate | TTBK1 Vina | TTBK1 Vinardo | MAO-B Vina | MAO-B Vinardo | Robust under both? |
|---|---|---|---|---|---|
| `Cc1ccc(-c2oc3cccc(O)c3c(=O)c2O)cc1Cc1ccccc1` | -8.34 | -6.91 | -10.32 | -9.51 | no (Vina ~#8) |
| **`Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`** | **-7.95** | **-6.37** | **-11.42** | **-8.73** | **yes (#2 / #2)** |
| `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1Cc1ccccc1` | -8.15 | -6.94 | -11.49 | -6.46 | no |

## 3.4 The lead engages both on-target pockets with distinct interaction patterns

ProLIF profiling of the consensus best poses gives 10 scored contacts at TTBK1 and 14 at
MAO-B. At MAO-B the lead engages the aromatic cage adjacent to the FAD cofactor
(TYR326/PHE343/TYR435) plus the hydrophobic substrate cavity (PHE168, LEU171, CYS172,
ILE198/199, GLN206, ILE316). At TTBK1 it occupies the ATP pocket, contacting
ILE40/ILE48/GLN110/GLY111/LEU175 with a single hydrogen bond to ASN113 as acceptor.

![Lead candidate pose in the TTBK1 ATP pocket (4NFM).]({{artifact:art_423d699f-ba77-4102-97ac-e9a9fafa8d75}})

![Lead candidate pose in the MAO-B substrate cavity (2V5Z), FAD-adjacent aromatic cage.]({{artifact:art_d01c04f7-08c0-4d7a-bd39-1a9714d31f2b}})

## 3.5 Both complexes are stable over 20 ns of molecular dynamics

Both systems equilibrated and remained stable. Last-quartile means (relative to the first
production frame):

| System | Protein RMSD (A) | Ligand RMSD (A) |
|---|---|---|
| TTBK1 (4NFM) | 1.93 | 1.52 |
| MAO-B (2V5Z) | 1.25 | 0.68 |

The lead is more mobile within the TTBK1 ATP pocket (1.52 A, reflecting early
repositioning from the docked pose) than in the MAO-B cavity (0.68 A), consistent with
the shallower, more solvent-exposed kinase hinge region versus the enclosed MAO-B substrate
cage. Neither ligand dissociated or left its pocket.

Hydrogen-bond occupancies are low at both targets, indicating predominantly hydrophobic
and van der Waals-driven binding rather than a persistent hydrogen-bond network. The
highest-occupancy contacts are GLN89 at TTBK1 (18.5% for the ligand O4 -> GLN89 backbone
carbonyl; 7.3% and 6.2% for GLN89 backbone amide to ligand O3/O4) and TYR433 at MAO-B
(6.6% for TYR433 hydroxyl -> ligand O3).

![Protein and ligand RMSD and per-residue RMSF over 20 ns for both complexes.]({{artifact:art_935b21b8-cb4c-46ab-867a-c965fa477d83}})

## 3.6 MM-GBSA confirms engagement of both targets and prefers MAO-B

| Component (kcal/mol) | TTBK1 (4NFM) | MAO-B (2V5Z) |
|---|---|---|
| van der Waals | -32.89 | -44.87 |
| Electrostatic | -7.99 | -3.48 |
| GB solvation | +19.37 | +17.28 |
| Nonpolar surface | -4.81 | -5.78 |
| **Total dG (mean +/- SEM)** | **-26.32 +/- 0.53** | **-36.85 +/- 0.23** |
| Total dG SD across frames | 5.33 | 2.30 |

Both are favorable. MAO-B is preferred by -10.5 kcal/mol, the
same direction as the docking comparison (-3.5
kcal/mol), giving cross-method agreement on the on-target preference even though the
absolute scales are not comparable. The difference is driven by van der Waals contact
(-44.9 vs. -32.9 kcal/mol), consistent with the enclosing MAO-B cavity versus
the shallower kinase hinge. The frame-to-frame SD is also more than twice as large at TTBK1
(5.3 vs. 2.3 kcal/mol), matching the higher ligand RMSD
there — the lead samples a broader range of binding geometries in the kinase pocket.

![MM-GBSA total binding free energy and per-component decomposition for both targets.]({{artifact:art_c3d42f7a-fca5-42f0-822f-f741c8df8271}})

## 3.7 The series does not discriminate either target from its paralog/isoform

| Comparison | Margin (kcal/mol) | Favorable? | Shortlist-wide |
|---|---|---|---|
| TTBK1 vs. TTBK2 (6U0K) | +1.71 | no | 0/15 favorable; mean +1.92 |
| TTBK1 vs. TTBK2 (7Q8Y) | +2.23 | no | 0/15 favorable; mean +2.27 |
| MAO-B vs. MAO-A | +0.09 | no (no separation) | 2/15 favorable; mean +0.45 |

Every candidate in the shortlist is predicted to bind TTBK2 more strongly than TTBK1, by
1.5-2.7 kcal/mol. This is expected rather than surprising: TTBK1 and TTBK2 share a
near-identical ATP-binding kinase fold — which is precisely why the 7Q8Y grid box could be
transferred onto 4NFM by superposition — and a flavonoid occupying the conserved hinge
region has little to discriminate on. The MAO-B/MAO-A result is likewise consistent with
the well-documented difficulty of achieving MAO-B-selective inhibition with planar
polyphenolic scaffolds, given the two isoforms' highly homologous active sites.

**These margins carry an important caveat: none of the three anti-target receptors carries
a passing redocking validation** (7Q8Y failed at 5.29 A; 6U0K and 2Z5X were never
validated, and 2Z5X required a custom FAD transplant from the 2V5Z template). The margins
should be read as directional signal, not validated affinity predictions. That said, the
signal points uniformly in the unfavorable direction across 15 candidates and three
receptors, which is harder to attribute to receptor-preparation error than a marginal
result would be.

![On-target versus anti-target docking, and cross-method on-target confirmation.]({{artifact:art_d0e48a02-1f63-464d-9ec5-1e718c481e4a}})

# 4. Discussion

## 4.1 What this study establishes

The lead 7-deoxyflavonol engages both TTBK1 and MAO-B, and this conclusion rests on two
methodologically independent lines of evidence rather than one: consensus docking with two
distinct scoring functions, and 20 ns MD with MM-GBSA free energies. Both agree on the
direction of on-target preference (MAO-B > TTBK1). Cross-method agreement of this kind is
the minimum standard for a docking-derived claim, and most published virtual-screening
reports of natural-product MTDLs do not meet it.

Two negative results are as informative as the positive one, and both are reported at
equal weight. The scaffold-level BBB ceiling (Section 3.1) is a hard structural constraint
on the entire flavonol class for CNS targets: it cannot be optimized away by decoration,
only by removing a core oxygen. Any CNS-directed flavonoid design campaign that reports BBB
permeability without confronting the TPSA floor of the native scaffold should be read
skeptically. Separately, the Vinardo cross-check (Section 3.3) shows that the fine-grained
docking rank order within an already narrow high-scoring slice of chemical space is not
robust to scoring-function choice — negatively correlated at MAO-B. Reporting a "top hit"
from single-scoring-function virtual screening without such a check risks selecting on
scoring-function idiosyncrasy.

## 4.2 The selectivity liability is the dominant risk

The pipeline's own counter-screen contradicts a selective-dual-target interpretation. The
honest reading is that this series is a **dual on-target engager with unresolved isoform
selectivity**, not a selective dual-target lead. This distinction matters
pharmacologically: MAO-A inhibition without MAO-B selectivity raises the tyramine
pressor ("cheese effect") liability that isoform-selective MAO-B inhibitors were developed
to avoid, and pan-TTBK activity would engage TTBK2, whose loss-of-function causes
spinocerebellar ataxia type 11, making unintended TTBK2 inhibition a specific safety
concern rather than a generic off-target note.

## 4.3 The dual-target rationale itself remains a hypothesis

No literature was found proposing or validating TTBK1 + MAO-B as a coherent dual-target
strategy. The two targets belong to different disease-mechanism families. A compound that
engages both would need a therapeutic context in which simultaneous tau-kinase and
MAO-B inhibition is desirable — plausible in principle for mixed-pathology
neurodegeneration, but not demonstrated. Readers should treat the dual-target framing as
the study's hypothesis, not its finding.

## 4.4 Limitations

1. **No experimental validation.** Every affinity reported is computed. No compound was
   synthesized, and no biochemical or cellular assay was performed.
2. **TTBK1's receptor is apo and unvalidatable.** 4NFM has no co-crystal ligand, so the
   docking protocol could not be validated on the actual on-target; the box was transferred
   from a TTBK2 structure that itself failed redocking validation.
3. **No anti-target receptor passed validation.** All selectivity margins inherit this.
4. **MM-GBSA is relative, not absolute.** Single-trajectory mode omits conformational
   entropy and receptor reorganization and overestimates binding magnitudes; the values
   support the TTBK1-vs-MAO-B comparison, not affinity prediction.
5. **Single 20 ns replicate, single candidate.** No replicate MD, no MD on anti-targets, no
   MD on the other shortlisted candidates — a consequence of the 6 GB single-GPU budget.
6. **BBB assessment is a 2D property model.** BOILED-Egg uses WLOGP/TPSA and does not model
   efflux transporters (notably P-glycoprotein), to which flavonoids are known substrates.
7. **Literature quantification is PubMed-only** (as queried 2026-08-17) and used name-string
   counts; OpenAlex was unavailable and IMPPAT was not queried, so compound provenance rests
   on primary-abstract cross-referencing against PubChem rather than a curated phytochemical
   database.
8. **One brief-supplied citation was mis-attributed** and is corrected here: the reference
   given as "Sharma 2020, *Molecules*" resolves by journal/volume/pages/topic to
   Chaurasiya et al. 2020, *Molecules* 25:5358.

## 4.5 Next steps

The decisive next computational step is extending MM-GBSA to TTBK2 and MAO-A complexes,
which would test whether the unfavorable docking-derived selectivity margins survive a
higher level of theory — the one place where the current evidence is weakest and where the
anti-target receptors' lack of validation matters most. Replicate MD on the lead, and MD on
the two benzyl/methyl regioisomers, would establish whether the on-target result is
candidate-specific or a scaffold property. Experimentally, the informative first assay is
not a single-target potency measurement but a side-by-side isoform panel: MAO-A/MAO-B and
TTBK1/TTBK2 in parallel, since selectivity and not potency is the open question.

# 5. Conclusions

A generative, BBB-constrained design campaign on a deoxygenated flavonol core produced a
4-trifluoromethyl/methyl-substituted 7-deoxyflavonol that engages both TTBK1 and MAO-B by
two independent computational methods, with MAO-B the more robustly engaged target
(-36.8 vs. -26.3 kcal/mol by MM-GBSA). The same pipeline shows
that neither target is discriminated from its closest paralog or isoform. Two methodological
negative results — the flavonol scaffold's absolute BBB ceiling and the failure of docking
rank order to survive a scoring-function change — are reported as primary findings. This
work is hypothesis-generating: it identifies a chemically tractable dual-engagement scaffold
and, equally, identifies isoform selectivity as the liability that any development of this
series must solve first.

# Data and code availability

All intermediate and final data products are available as project artifacts, including the
flavonoid library with provenance and confidence tiers, the full generative-run output, the
filtering-cascade tables, consensus and Vinardo docking results for on- and anti-targets,
prepared receptor and ligand files, MD trajectory analyses, and raw MMPBSA.py output for
both systems.

# References

1. Halkina T, et al. Discovery of potent and brain-penetrant tau tubulin kinase 1 (TTBK1)
   inhibitors that lower tau phosphorylation in vivo. *J Med Chem*. 2021;64(9):6358-6380.
   PMID 33944571.
2. Binda C, et al. Structures of human monoamine oxidase B complexes with selective
   noncovalent inhibitors: safinamide and coumarin analogs. *J Med Chem*. 2007;50(23):
   5848-5852. (PDB 2V5Z, 2V60.)
3. Son SY, et al. Structure of human monoamine oxidase A at 2.2 A resolution: the control of
   opening the entry for substrates/inhibitors. *Proc Natl Acad Sci USA*. 2008;105:5739-5744.
   (PDB 2Z5X.)
4. Daina A, Zoete V. A BOILED-Egg to predict gastrointestinal absorption and brain
   penetration of small molecules. *ChemMedChem*. 2016;11(11):1117-1121.
5. Quiroga R, Villarreal MA. Vinardo: a scoring function based on AutoDock Vina improves
   scoring, docking, and virtual screening. *PLoS One*. 2016;11(5):e0155183.
6. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new
   scoring function, efficient optimization, and multithreading. *J Comput Chem*.
   2010;31(2):455-461.
7. Loeffler HH, He J, Tibo A, et al. REINVENT 4: modern AI-driven generative molecule design.
   *J Cheminform*. 2024;16:20.
8. Miller BR III, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE. MMPBSA.py: an
   efficient program for end-state free energy calculations. *J Chem Theory Comput*.
   2012;8(9):3314-3321.
9. Eastman P, et al. OpenMM 7: rapid development of high performance algorithms for
   molecular dynamics. *PLoS Comput Biol*. 2017;13(7):e1005659.
10. Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O. MDAnalysis: a toolkit for the
    analysis of molecular dynamics simulations. *J Comput Chem*. 2011;32(10):2319-2327.
11. Bouysset C, Fiorucci S. ProLIF: a library to encode molecular interactions as
    fingerprints. *J Cheminform*. 2021;13:72.
12. Roe DR, Cheatham TE III. PTRAJ and CPPTRAJ: software for processing and analysis of
    molecular dynamics trajectory data. *J Chem Theory Comput*. 2013;9(7):3084-3095.
13. Chaurasiya ND, et al. Selective inhibition of human monoamine oxidase B by
    O-methylated flavonoids. *Molecules*. 2020;25(22):5358. (Cited in the project brief as
    "Sharma 2020"; attribution corrected here.)
14. Primary phytochemical isolation reports for *Evolvulus alsinoides* flavonoid
    constituents: PMID 17473466, PMID 19748554, PMID 23357036. (Full bibliographic details
    to be completed from the PubMed records at submission; only PMIDs were captured during
    library assembly.)
