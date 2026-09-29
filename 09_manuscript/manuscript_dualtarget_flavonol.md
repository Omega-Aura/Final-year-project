---
title: "A generative 7-deoxyflavonol against TTBK1 and MAO-B: isoform selectivity resolvable for one target pair and not the other, and why docking rank does not predict pose stability"
subtitle: "A hypothesis-generating computational study"
date: 2026-09-27
---

# Abstract

**Background.** Tau-tubulin kinase 1 (TTBK1) and monoamine oxidase B (MAO-B) are both CNS
drug-discovery targets in neurodegeneration, but they belong to distinct disease mechanisms and
have never been proposed as a combined dual-target pair in the published literature. Flavonoids
are documented MAO-B inhibitors and promiscuous kinase binders, making them a plausible starting
scaffold for a multi-target-directed ligand (MTDL) hypothesis. Both targets have a close
paralog or isoform — TTBK2 and MAO-A — whose inhibition carries a specific, named liability, so
selectivity rather than potency is the load-bearing question.

**Methods.** A flavonoid library from *Evolvulus alsinoides* seeded reinforcement-learning
molecular generation (REINVENT4 LibInvent) on a 7-deoxyflavonol core under a four-component
geometric-mean reward (TTBK1 docking, MAO-B docking, synthetic accessibility, BOILED-Egg
blood–brain-barrier signed distance). Candidates were filtered on Lipinski, BBB, GI-absorption
and PAINS/BRENK criteria, then re-docked with AutoDock Vina at exhaustiveness 32 over three
independent seeds against six receptors — TTBK1 (7JXX, 7Q8V), TTBK2 (7Q8Y), MAO-B (2V5Z),
MAO-A (2Z5X) and 4BTK — **all six validated by native-ligand redocking**, and with the
crystallographic water shell made symmetric between any two receptors being compared. The
scoring function was calibrated against a compound crystallised in both TTBK paralogs with
same-assay measured potencies. The lead was advanced to **seventeen** 10.1 ns explicit-solvent MD
simulations: **three independent docked poses in each of the four targets**, plus a velocity-seed
replicate of the best on-pose trajectory in **every one of the four arms**. MM-GBSA was computed **only on trajectories that held
their docked pose**.

**Results.** The native flavonol scaffold has an absolute BBB ceiling — minimum TPSA 90.9 Å²
across 666 generated candidates, a floor set by the core oxygens, giving a null shortlist.
Removing the 7-OH recovered 56 candidates passing every filter. All six receptors reproduced
their crystallographic pose (0.27–1.45 Å). **Docking rank did not predict which pose survived
dynamics: the top-ranked pose was the most stable one in one of four targets, and the
best-holding MAO-A pose was its worst-ranked.** Nine of seventeen runs held their pose and two
dissociated outright. For the MAO isoform pair the intended selectivity is supported in both
direction and size — all twelve on-pose pairings favour MAO-B, and velocity replicates put
ΔΔG at **2.55 kcal/mol against a measured run-to-run spread of 0.24–0.26**, with docking
agreeing independently (−3.28 kcal/mol mean over 48 quantifiable candidates, favourable 48/48).
For the TTBK paralog pair there is **no resolvable discrimination**: with both arms replicated the
ΔΔG is **0.53 kcal/mol** against measured run-to-run spreads of **2.08 and 3.97 kcal/mol**, pose
stability is identical (1 of 3 poses each), and the protocol carries a calibrated
**~1.0 ± 0.25 kcal/mol systematic bias toward TTBK2**, placing any margin below ~1.3 kcal/mol
inside its own bias. Separately, `MMPBSA.py`'s reported standard error is **not** a usable error
bar, and across five replicate pairs the true spread rises monotonically with pose instability —
2× the reported SEM for tightly held poses, 33× for a pair that drifted and dissociated.

**Conclusions.** This 7-deoxyflavonol engages both on-targets, and the two target pairs must be
reported separately because they reach opposite conclusions: MAO-B over MAO-A is supportable in
direction and magnitude, while TTBK1 over TTBK2 is not supportable at all and the pipeline
lacks the resolution to speak on it. Two methodological findings — that docking rank does not
predict pose stability, and that a single-trajectory SEM understates the true uncertainty by up
to 33-fold in a way that tracks pose instability — bound what any study of this design may
claim. This is a hypothesis-generating computational study; no affinity reported here is a
measured value.

**Keywords:** TTBK1, MAO-B, flavonol, multi-target-directed ligand, generative reinforcement
learning, MM-GBSA, pose stability, anti-target counter-screening, blood–brain barrier,
protocol calibration

---

# 1. Introduction

Tau-tubulin kinase 1 (TTBK1, UniProt Q5TCY1) phosphorylates tau at disease-relevant epitopes
including Ser422 and is genetically and biochemically implicated in Alzheimer's disease and
other tauopathies. TTBK1 chemical-probe development remains early-stage: the
best-characterised inhibitor series with public measured potencies is the brain-penetrant
azaindazole/pyrrolopyridine series of Halkina et al. (2021), whose lead lowers tau pSer422
in vivo.

Monoamine oxidase B (MAO-B, UniProt P27338) is by contrast a clinically validated Parkinson's
disease target. It oxidatively deaminates dopamine and other monoamines, producing H₂O₂ and
contributing to oxidative stress and dopaminergic neurodegeneration, and is upregulated in
reactive astrocytes in neurodegenerative tissue. Approved inhibitors provide well-measured
reference pharmacology, though — as Section 2.2 records — only with their assay attached.

**The dual-target hypothesis examined here is novel and unvalidated.** Systematic PubMed
querying returned zero records combining flavonoids with TTBK1, and zero proposing a
TTBK1 + MAO-B combination in any context. Flavonoids as MAO-B inhibitors is by contrast an
established sub-field. The two targets are not co-implicated in any single validated disease
mechanism in the retrieved literature. The dual-target framing therefore rests on three weaker
premises, stated plainly: both are licensed CNS neurodegeneration targets; flavonoids are
documented to be promiscuous across kinase and MAO pharmacology; and MTDL design combining
kinase inhibition with MAO-B inhibition has general precedent in neurodegenerative drug design.
No claim of an established dual-target consensus is made or implied.

Two constraints shaped the study, and a third shaped its methodology.

First, both targets are CNS targets, so blood–brain-barrier permeability was treated as a hard
requirement rather than a post-hoc annotation.

Second, **both targets have a close relative whose inhibition is a specific liability, not a
generic off-target note.** TTBK2 loss-of-function causes spinocerebellar ataxia type 11, so
unintended TTBK2 inhibition is a defined safety concern; MAO-A inhibition without MAO-B
selectivity raises the tyramine pressor ("cheese effect") liability that isoform-selective MAO-B
inhibitors were developed to avoid. Selectivity was therefore tested explicitly rather than
asserted, and — the point this paper turns on — tested at a level of theory sufficient to know
whether the test could resolve the effect at all.

Third, and arising during the work: a computational selectivity claim is only as good as the
symmetry of the comparison that produced it and the error bar placed on it. This study reports
two findings of that kind at equal weight with its chemistry, because in the course of the work
each of them silently inverted a conclusion before being caught (Sections 3.6, 3.9 and 4.3).

# 2. Methods

## 2.1 Species selection, natural-product library and generative design

Six Indian medicinal species with traditional neurological indications were ranked by a
pre-registered multiplicative rule, `flavonoid_coverage × literature_gap`. The multiplicative
form is deliberate: a species with a large literature gap but no structurally verified
flavonoids scores zero, since a gap without dockable chemistry is not an opportunity.
*Evolvulus alsinoides* ranked first (composite 0.933), carrying six structurally confirmed
flavonoid/flavonol compounds traceable to primary isolation papers, only four CADD-related
records, and zero records combining the species with either target. Its dominant scaffold class
(3-hydroxyflavone) set the scaffold constraint. The library comprised 24 entries, all parsing in
RDKit, with per-compound confidence tiers recorded from source provenance.

Generation used REINVENT4 LibInvent `staged_learning` (DAP, σ = 128, rate 1×10⁻⁴) with an
`IdenticalMurckoScaffold` diversity filter (bucket 25, minscore 0.4) to stop convergence onto a
single scaffold. Reward components were TTBK1 docking (0.35), MAO-B docking (0.35), synthetic
accessibility (0.15) and BOILED-Egg BBB signed distance (0.15), **aggregated as a geometric
mean** so that a molecule cannot win by excelling at one target and failing the other — an
arithmetic mean permits exactly that trade-off. Docking inside the reward used Vina at
exhaustiveness 4 for tractability; all reported scores come from the independent re-dock at
exhaustiveness 32 (Section 2.4).

An initial campaign on the native flavonol core preceded this one and failed the BBB gate
outright (Section 3.1). The scaffold was then redesigned to a **mono-deoxy 7-deoxyflavonol**,
retaining the 5-OH/4-carbonyl intramolecular hydrogen bond and removing the 7-OH. An early
scaffold construction that inadvertently retained 7-OH and removed 5-OH was detected and the
32-step run re-executed from scratch on the corrected regiochemistry.

## 2.2 Receptors, and the reference set with its assays

Six structures were prepared by a single script in one run, changing only the filename, so that
any two receptors whose scores are compared are treated identically:

| Role | Protein | PDB | Native ligand |
|---|---|---|---|
| **On-target** | TTBK1 | **7JXX** | VP7 |
| On-target, calibration arm | TTBK1 | 7Q8V | 9IV |
| **Anti-target** | TTBK2 | **7Q8Y** | 9IV |
| **On-target** | MAO-B | **2V5Z** | SAG (safinamide Schiff base) |
| **Anti-target** | MAO-A | **2Z5X** | HRM (harmine) |
| Cross-check | — | 4BTK | DTQ |

This replaces the apo TTBK1 structure used in an earlier pass of this work: 7JXX is
drug-bound, so the docking protocol can be validated on the actual on-target rather than having
a grid box transferred onto it by superposition.

**Water shell.** Crystallographic waters were handled by an explicit minimal-water protocol:
start dry and add waters back only until the native redock passes. The rule adopted after
Section 3.4 is that **any two receptors whose scores will be compared must carry the same water
treatment**, because an asymmetric shell alone produces an apparent selectivity margin.

**Reference compounds and their assays.** Every value in the reference set was traced to a
primary source. Three provisos travel with the set and are stated here rather than in a
footnote:

- The MAO-B reference safinamide is **not the same molecule as the PDB ligand SAG** in 2V5Z,
  which is its Schiff base; both are in the set, separately.
- Safinamide's IC50 of 7.67 nM (Stössel et al. 2013, recombinant human MAO-B, p-tyramine
  substrate) coexists with a *K*i of 0.1–0.5 µM from Binda et al. 2007 — the very paper behind
  the 2V5Z structure used here, 13–65× weaker. Neither is wrong; they are different
  measurements. **A potency value is only meaningful with its assay attached.**
- The score-versus-potency correlation set is **mixed-species**: two of its points (kaempferol,
  isatin) are rat-brain assays and are the two flavonoids, i.e. the compounds closest to this
  study's lead. It is reported as mixed-species, never as a human-target correlation.
  Selegiline and rasagiline are excluded: both are irreversible covalent inhibitors whose IC50
  is preincubation-time dependent and therefore not an equilibrium constant, and the available
  values are review-compiled rather than primary. The TTBK1 reference DTQ is likewise excluded,
  and its 240 nM figure is a surface-plasmon-resonance *K*d — the enzymatic IC50 for the same
  compound is 4610 nM, 19× weaker — so it must not be quoted as an IC50.

## 2.3 Protocol validation and scoring-function calibration

Two independent validations were run before any candidate score was interpreted.

**Geometry.** Each receptor's own crystallographic ligand was redocked into its own site, with
heavy-atom RMSD computed by `spyrmsd` (symmetry-corrected; plain atom-order RMSD inflates values
for symmetric groups). Bond orders in ligands extracted from PDB files are unreliable and were
repaired before RMSD, since a wrong bond order changes the symmetry graph and therefore the
number.

**Ranking.** Geometry can be right while ranking is backwards, so the scoring function was
calibrated on the compound 9IV (VNG2.73), which is crystallised in **both** TTBK paralogs —
7Q8V (TTBK1) and 7Q8Y (TTBK2) — and whose potency against both was measured **in one assay in
one paper** (Nozal et al. 2022: IC50 330 nM TTBK1, 490 nM TTBK2). Both arms of a calibration
pair must come from the same assay or the comparison measures the difference between two
laboratories; an earlier version of this calibration mixed a cross-paper midpoint against a
single-paper value and is corrected here.

The assumed uncertainty on the experimental ΔΔG is a named parameter, not a measurement:
Nozal et al. publish no error on either IC50, so a conventional 30% within-assay relative
precision is propagated through the log. That assumed term is **16× larger than the docking
SEM**, which means the calibration's precision is set by the literature value and no amount of
further docking would tighten it.

## 2.4 Docking

AutoDock Vina, `--exhaustiveness 32 --num_modes 9`, three independent seeds (11/22/33), each
from an independent conformer embedding rather than a re-seeded search on one geometry. Reported
per ligand: best score over all seeds and poses, consensus (mean of per-seed best pose), and
inter-seed SD. Seed-to-seed spread is typically ≤ 0.1 kcal/mol, so consensus with a reported SD
makes it evident when a difference is smaller than the search noise.

Selectivity margin is `on-target consensus − anti-target consensus`, so a negative margin is
favourable. **Candidates with no viable pose in an anti-target are excluded from the margin
statistics and reported separately.** Vina returns a "best" mode even when it cannot place a
ligand — near-zero or frankly positive — and averaging such a non-measurement into a ΔΔG
manufactures selectivity out of a docking failure. Eight of the 56 candidates fail this way in
MAO-A; the margin is therefore quoted over the 48 quantifiable candidates, with the eight
reported as a separate qualitative result.

An independent scoring-function cross-check (Vinardo) was run in an earlier pass of this work,
on the earlier receptor set. It has **not** been repeated on the validated receptors, and
Section 3.3 is reported with that limitation explicit.

## 2.5 Molecular dynamics

Seventeen systems, each built independently and self-contained (own protein PDB, own ligand
parameters, own `tleap` input, topology and trajectory): `cand_003` in TTBK1 (7JXX), TTBK2
(7Q8Y), MAO-A (2Z5X) and MAO-B (2V5Z), from **three independent docked poses per target**, plus a
velocity-seed replicate of the best on-pose trajectory in **each of the four arms**. A replicate
re-runs byte-identical topology and coordinates and differs only in the random velocity seed,
which the integrator draws itself.

Protein force field AMBER ff14SB; ligand and cofactor GAFF2 with `acpype` charges; TIP3P water
with neutralising Na⁺/Cl⁻; systems built with `tleap`, protein prepared with `pdb4amber`
(HIS→HIE) and hydrogens added with OpenMM `Modeller`. Propagated in OpenMM on the CUDA
platform: PME, 1.0 nm cutoff, HBonds constrained, LangevinMiddle integrator at 300 K with
1 ps⁻¹ friction and a 2 fs timestep, Monte Carlo barostat, 2000 minimisation steps, 100 ps
equilibration, then **10.1 ns production** with a frame every 10 ps (1000 frames). All sixteen
runs completed; `tleap` reported zero errors and every system was neutral.

Three design choices each exist because of a defect found earlier in this work, and each is
load-bearing for the results:

1. **Three poses per protein, not one.** An earlier run tested TTBK2 from three poses but TTBK1
   from one pose twice, and produced an apparent "TTBK2 is less stable" result that was purely
   an artifact of comparing a *pose scan* against a *replicate scan*. Made symmetric, the
   difference vanished.
2. **A velocity replicate of the best on-pose trajectory in every arm.** The reported SEM is not
   the error bar (Section 3.9); the run-to-run spread is, and it is knowable only from a
   replicate. Replicating only three of the four arms would have left one comparison judged by
   the other arm's uncertainty — and when the fourth was run, it moved the central value as well
   as supplying the error bar (Section 3.8).
3. **Ligand parameters generated from the docked pose**, not from a relaxed free conformer, so
   the MD starting geometry is the geometry the docking produced. Alternative poses are
   transferred onto the parameter template by positional coordinate swap — never re-embedded —
   and the transfer refuses to proceed unless pose 1 round-trips onto the template (it does, at
   0.0002 Å) and the element sequences agree, so a silent atom-order mismatch cannot yield a
   scrambled ligand that still builds cleanly.

**The FAD cofactor** was parameterised as an 84-atom GAFF2 residue from RCSB ideal chemistry
with **zero guessed parameters**, and positionally restrained to its crystallographic
coordinates. It is treated as part of the receptor on both sides of the MM-GBSA subtraction, not
scored as ligand. The 8α-S-cysteinyl covalent bond is **not** modelled; this is an
approximation, labelled as one, and its consequences are in Section 4.4.

Analysis used `cpptraj`: ligand RMSD after a backbone fit, per-residue RMSF, and — for the MAO
systems — a refit on the ordered core, because both MAO constructs end in a solvent-exposed
C-terminal tail that has no membrane to sit in here and dominates any whole-protein fit.
**Stability was judged on the mean of the final 100 frames**, not the whole-run mean: a ligand
that leaves late still shows a low whole-run average.

## 2.6 MM-GBSA

AmberTools `MMPBSA.py` v14.0, single-trajectory Generalized Born, `igb=5`, `saltcon=0.150`,
LCPO surface area, topologies from `ante-MMPBSA.py` with `--radii=mbondi2`. Every fifth frame of
1000 gives ~200 frames 10 ps apart.

**MM-GBSA was computed only where it means something.** No reported binding energy comes from a
trajectory that left its docked pose: an energy averaged over such a trajectory describes a
structure that was never docked. Three off-pose runs do carry computed values, retained
deliberately as the evidence for Sections 3.6 and 3.9, and none of them is quoted as a binding
energy.

## 2.7 What this protocol cannot do

Stated here rather than only in the Discussion, because it governs how the Results must be read.
All reported affinities are computed. Single-trajectory MM-GBSA neglects conformational entropy
and receptor reorganisation and overestimates the magnitude of binding free energies; values are
interaction-energy estimates, not thermodynamic binding free energies, and are not convertible
to *K*d. **Absolute values are not comparable between proteins** — protein-specific desolvation
and surface terms do not cancel — so a MAO number may not be set against a TTBK number. ΔΔG
*within* a target pair is the only free-energy quantity this study interprets.

# 3. Results

## 3.1 The native flavonol scaffold has an absolute BBB ceiling

The first generative campaign produced 1,280 scored molecules, reducing to 666 unique, valid,
novel candidates. The cascade gave 587 passing Lipinski, 198 passing GI absorption and 217
alert-free — but **zero** passing the BBB gate.

The cause is structural, not statistical. Minimum TPSA across all 666 candidates was
**90.9 Å²**, a floor set by the flavonol core itself: the mandatory 3-OH, 4-C=O and 5-OH oxygens
plus the pyranone ring oxygen. The BBB region requires roughly TPSA < 79 Å² at the relevant
WLogP, so **no B-ring decoration can bring a native flavonol into the permeant region.** Both
natural seed compounds and the literature TTBK1/2 reference inhibitor also fall outside. Under a
literal BBB gate the formal output of this campaign is a null shortlist.

## 3.2 Scaffold deoxygenation converts a null result into 56 candidates

Removing the 7-OH while retaining the 5-OH/4-carbonyl hydrogen bond, and adding a direct BBB
signed-distance reward term, produced 1,024 candidates over 32 steps (405 unique, valid, novel).
Mean composite reward rose from 0.08 at step 1 to a 0.45–0.53 plateau by steps 20–32.

The BBB pass rate over the run is informative about the trade-off rather than a clean success
curve: steps 1–2, still near the unmodified prior, start at 71.9% and 62.5%; steps 3–8 collapse
to 0–34.4% (mean 16.9%) as reinforcement learning optimises affinity at the expense of
permeability; from step ~9 the BBB term pulls the rate back, averaging 45% over steps 13–32 but
still ranging 12.5–65.6% between batches, with no further downward trend.

Final cascade over the 405-molecule pool: 369 Lipinski-passing, 81 BBB-passing, 334 GI-passing,
253 alert-free, and **56 candidates passing all filters simultaneously** — the qualitative
reversal the redesign targeted. The 56 span MW 310–420 Da, TPSA 70.7–77.2 Å² (inside the BBB
region), WLogP 3.6–4.8, and 8 unique Murcko scaffolds.

## 3.3 The lead compound, and the one cross-check that is missing

The compound carried forward is `cand_003`:

```
Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F
```

a 4-trifluoromethyl/methyl-decorated 7-deoxyflavonol.

| MW | WLogP | TPSA | HBD | HBA | Lipinski | BBB | GI | Alerts |
|---|---|---|---|---|---|---|---|---|
| 336.27 | 4.20 | 70.67 Å² | 2 | 4 | 0 violations | pass | pass | 0 |

Consensus docking over three seeds on the validated receptors: **−8.39 kcal/mol at TTBK1
(7JXX), rank 8 of 56**, and **−11.41 kcal/mol at MAO-B (2V5Z), rank 2 of 56**.

**The limitation to state here rather than bury.** In an earlier pass of this work the top-15
candidates were rescored with an independent scoring function (Vinardo), which agreed poorly
with Vina — Pearson r = 0.56 at TTBK1 and **−0.28** at MAO-B — and that negative result was a
primary finding: within an already narrow high-scoring slice of chemical space, fine-grained
docking rank order does not survive a change of scoring function, even though all candidates
score favourably in absolute terms under both. That cross-check was run on the **earlier
receptor set** and has not been repeated on the validated receptors reported here. The
scoring-function caveat therefore stands as a general warning supported by this project's
earlier data, but it is **not** a cross-check of the present rankings, and `cand_003`'s rank of
8 of 56 at TTBK1 should be read with that gap in mind.

## 3.4 All six receptors reproduce their crystallographic pose

| Receptor | Protein | Native ligand | Best RMSD | Result |
|---|---|---|---|---|
| 2Z5X | MAO-A | HRM | **0.27 Å** | pass |
| 7JXX | TTBK1 | VP7 | **0.55 Å** | pass |
| 2V5Z | MAO-B | SAG | **0.66 Å** | pass |
| 7Q8V | TTBK1 | 9IV | **0.68 Å** | pass |
| 4BTK | — | DTQ | **0.79 Å** | pass |
| 7Q8Y | TTBK2 | 9IV | **1.45 Å** | pass |

All six pass, four under 0.7 Å. **Pose-level conclusions therefore rest on solid ground** — a
material improvement on the earlier pass of this work, in which only one receptor carried a
passing redocking validation and the TTBK1 structure was apo and could not be validated at all.

## 3.5 The scoring function carries a measured ~1.0 kcal/mol bias toward TTBK2

Geometry being right does not make ranking right. On the 9IV pair, crystallised in both paralogs
and measured in one assay:

| Quantity | Value |
|---|---|
| Docking consensus, 7Q8V (TTBK1) | −8.115 ± 0.015 |
| Docking consensus, 7Q8Y (TTBK2) | −8.899 ± 0.022 |
| **Docking margin** | **+0.784 kcal/mol favouring TTBK2** |
| **Experimental ΔΔG** | **−0.234 kcal/mol** (330 nM vs 490 nM, same assay) |
| **Systematic bias** | **~1.0 ± 0.25 kcal/mol** |
| **Interpretability floor** | **~1.3 kcal/mol** (bias + 1 SD) |

Experiment says 9IV slightly prefers **TTBK1**. The protocol reports the opposite, by about
1 kcal/mol, on a compound where the true difference is small and points the other way.

**Consequence, and it governs Section 3.8:** any TTBK2-over-TTBK1 margin smaller than
~1.3 kcal/mol is inside the protocol's demonstrated bias and cannot be read as selectivity. The
bias is quoted as ~1.0 and never to three figures: its inputs are two unreplicated IC50s, and a
precise-looking number with its real uncertainty unstated is the specific failure this study
documents in Section 3.9.

## 3.6 The TTBK selectivity margin is largely a property of the water shell

On the validated receptor pair the TTBK1-vs-TTBK2 margin is +0.12 kcal/mol. That comparison is
confounded: it pairs a receptor carrying **zero** crystallographic waters against one carrying
**five**. Since waters were independently shown to penalise these candidates by of order
1 kcal/mol — individual ligands move by up to 3 kcal/mol depending purely on whether waters were
retained, with a mean wet-minus-dry shift of only −0.08 kcal/mol but a range running to +3.04 —
that asymmetry is of the same size as the effect being measured.

Re-run **water-symmetrically**, the margin reverses and grows: TTBK2 is favoured by
~1.6 kcal/mol, with **0 of 56 candidates** favourable toward TTBK1. Subtracting the calibrated
bias of Section 3.5 leaves ~0.6 kcal/mol, comfortably inside the interpretability floor.

For the MAO pair the comparison is water-symmetric by construction — both receptors carry zero
waters — so its margin is not subject to this confound: **−3.28 kcal/mol mean over the 48
quantifiable candidates, favourable in 48 of 48**, range −6.04 to −1.80. The other eight
candidates, the largest in the set (MW 339–420), have **no viable MAO-A pose at all** while
still binding MAO-B at −8.7 to −11.3 kcal/mol. That is a qualitative selectivity result, and
arguably a stronger one than a ΔΔG, but it is not a ΔΔG and is not averaged into one.

## 3.7 Docking rank does not predict which pose survives dynamics

![Pose stability across all seventeen MD systems](../10_results/fig1_pose_stability.png)

Ligand RMSD from the docked pose, mean over the final 100 frames of 10.1 ns:

| System | Target | Structure | Start | Last-100 RMSD | Verdict |
|---|---|---|---|---|---|
| `system` | TTBK1 | 7JXX | pose 1, run 1 | 1.88 Å | stable |
| `system_TTBK1_r2` | TTBK1 | 7JXX | pose 1, run 2 | 1.38 Å | stable |
| `system_TTBK1_p2` | TTBK1 | 7JXX | pose 2 | 4.18 Å | drifts |
| `system_TTBK1_p3` | TTBK1 | 7JXX | pose 3 | 3.42 Å | drifts |
| `system_TTBK2` | TTBK2 | 7Q8Y | pose 1, run 1 | 3.72 Å | drifts |
| `system_TTBK2m` | TTBK2 | 7Q8Y | pose 1, run 2 | 5.69 Å | dissociates |
| `system_TTBK2_p2` | TTBK2 | 7Q8Y | pose 2, run 1 | 1.99 Å | stable |
| `system_TTBK2_p2_r2` | TTBK2 | 7Q8Y | pose 2, run 2 | 1.89 Å | stable |
| `system_TTBK2_p3` | TTBK2 | 7Q8Y | pose 3 | 5.79 Å | dissociates |
| `system_MAOA` | MAO-A | 2Z5X | pose 1 | 3.24 Å | drifts off pose |
| `system_MAOA_p2` | MAO-A | 2Z5X | pose 2 | 2.43 Å (1.79 core-fit) | holds, loose |
| `system_MAOA_p3` | MAO-A | 2Z5X | pose 3, run 1 | 1.43 Å | stable |
| `system_MAOA_p3_r2` | MAO-A | 2Z5X | pose 3, run 2 | 1.07 Å | stable |
| `system_MAOB` | MAO-B | 2V5Z | pose 1 | 2.28 Å | stable |
| `system_MAOB_p2` | MAO-B | 2V5Z | pose 2, run 1 | 1.10 Å | stable |
| `system_MAOB_p2_r2` | MAO-B | 2V5Z | pose 2, run 2 | 1.37 Å | stable |
| `system_MAOB_p3` | MAO-B | 2V5Z | pose 3 | 2.57 Å | holds, loose |

**Nine of seventeen runs hold their pose; two leave the site entirely.** With three poses in each
of four targets, **the top-ranked docked pose was the most stable one in one of four targets**
(TTBK1, the exception). The best-holding MAO-A pose is its *worst*-ranked of three — pose 3, Vina
−7.40, against pose 1's −8.24. This is the most strongly supported methodological result in the
study, and the reason every energy below is tied to a named pose.

Protein backbones were stable throughout: 1.34–1.90 Å across the eight TTBK runs, and
1.06–1.77 Å core-fit across the eight MAO runs. The two highest whole-protein figures — 3.64 Å
and 3.44 Å, both MAO-A — are the solvent-exposed C-terminal tail flailing without a membrane,
not unstable folds: per-residue RMSF puts the maximum at residues 506–513 of 513 in both, and
their cores sit at 1.50 and 1.28 Å. The more instructive of the two carries the study's **most
stable ligand (1.07 Å) alongside its highest backbone number** — a whole-protein RMSD says
nothing on its own about the binding site. The restrained FAD stayed within 0.46–0.73 Å of its
crystallographic position across all eight MAO systems, so the flavin wall of the cavity is
reproduced.

**Pose stability is also the discriminating filter no docking score supplied.** Two poses that
docking ranked at or near the top dissociated outright, and a pose docking ranked last held
tightest.

## 3.8 MAO-B over MAO-A is supported in direction and in size; TTBK1 over TTBK2 is not supported at all

![MAO binding free energies, MAO-A against MAO-B](../10_results/fig2_mao_binding_energy.png)

| System | Pose | ΔG (kcal/mol) | SD | Usable for selectivity? |
|---|---|---|---|---|
| TTBK1 | pose 1, run 1 | −31.36 | 3.82 | yes |
| TTBK1 | pose 1, run 2 | −33.44 | 2.63 | yes |
| **TTBK2** | **pose 2, run 1** | **−30.95** | 3.81 | yes — best on-pose TTBK2 |
| **TTBK2** | **pose 2, run 2** | **−34.92** | 3.14 | yes — its velocity replicate |
| MAO-A | pose 2 | −35.54 | 2.73 | direction only |
| **MAO-A** | **pose 3, run 1** | **−36.67** | 2.40 | yes — best on-pose MAO-A |
| **MAO-A** | **pose 3, run 2** | **−36.40** | 2.38 | yes — its velocity replicate |
| MAO-B | pose 1 | −39.48 | 2.20 | yes |
| **MAO-B** | **pose 2, run 1** | **−39.20** | 2.13 | yes — best on-pose MAO-B |
| **MAO-B** | **pose 2, run 2** | **−38.97** | 2.11 | yes — its velocity replicate |
| MAO-B | pose 3 | −39.04 | 2.86 | direction only |

**MAO-A versus MAO-B — a direction and a size.** On-pose, MAO-A lands at
−35.54/−36.67/−36.40 and MAO-B at −39.48/−39.20/−38.97/−39.04, so **all twelve on-pose pairings
favour MAO-B, by 2.30 to 3.94 kcal/mol** — the intended direction, MAO-A being the anti-target.
Docking independently favours MAO-B by 3.03–3.22 kcal/mol pose for pose on this ligand. Velocity
replicates of both best poses then supply the error bar: comparing replicate means, −36.53
against −39.09 is a **ΔΔG of 2.55 kcal/mol, 9.8× the larger of the two measured spreads
(0.24 and 0.26 kcal/mol)**. Both replicates also reproduced the pose stability (1.07 and
1.37 Å). For this pair, neither "no discrimination" nor "no *significant* discrimination"
survives.

One caveat travels with the magnitude: FAD is positionally restrained in both arms, which
suppresses receptor motion and damps run-to-run variation. **0.24–0.26 kcal/mol is the correct
error bar for this protocol, not evidence that MM-GBSA is intrinsically this precise.** The ΔΔG
is fair because both arms carry the same restraint — the symmetry is what makes it meaningful.

**TTBK1 versus TTBK2 — no resolvable discrimination, with both arms replicated.** Each kinase arm
has a velocity replicate of its best on-pose trajectory, so the comparison is symmetric in
uncertainty as well as in design. TTBK1's replicate mean is **−32.40** (spread **2.08**); TTBK2
pose 2's is **−32.93** (spread **3.97**). Comparing replicate means gives a ΔΔG of
**0.53 kcal/mol — 7.4× smaller than the larger of the two measured spreads.** Pose stability says
the same: 1 of 3 poses holds on each kinase under a symmetric test, and TTBK2 pose 2 held in both
replicates (1.99 and 1.89 Å). Taken with the ~1.0 ± 0.25 kcal/mol calibrated bias of Section 3.5
and the water-shell dependence of Section 3.6, four independent lines agree that **dynamics lacks
the resolution to speak on this pair at all.**

**Replicating the second arm did not merely add an error bar — it moved the central value.** With
a single TTBK2 trajectory the ΔΔG was 1.45 kcal/mol, judged against TTBK1's spread alone; that is
a two-arm difference assessed with one arm's uncertainty, the same class of asymmetry catalogued
in Section 4.3. Measured symmetrically the difference is a third of that size, and the honest
outcome is reported as a finding rather than buried: the original 1.45 kcal/mol quoted against the
reported SEM of 0.27 would have looked like a five-sigma selectivity result.

**The two target pairs therefore reach opposite conclusions and are not reported as one
finding.**

## 3.9 The reported standard error is not the error bar, and its failure scales with pose instability

![Reported SEM against the spread measured by re-running with a new velocity seed](../10_results/fig3_sem_vs_replicate.png)

`MMPBSA.py` computes its standard error as if ~200 frames sampled 10 ps apart were independent
draws. They are not. Across five replicate pairs — each two runs differing **only** in the
random velocity seed — the reported SEM stays in a narrow 0.15–0.27 kcal/mol band while the
measured spread ranges from 0.24 to 6.89 kcal/mol, a **29-fold range**.

The discrepancy is not random. **Ordered by how well the pair held its pose, the measured spread
rises monotonically across all five pairs:**

| Pair | Mean ligand RMSD | Measured spread | Ratio to reported SEM |
|---|---|---|---|
| MAO-B pose 2 | 1.24 Å | 0.24 | 2× |
| MAO-A pose 3 | 1.25 Å | 0.26 | 2× |
| TTBK1 pose 1 | 1.63 Å | 2.08 | 8× |
| TTBK2 pose 2 | 1.94 Å | **3.97** | 15× |
| TTBK2 pose 1 | 4.71 Å | 6.89 | 33× |

The fifth pair is an out-of-sample test of the relationship: TTBK2 pose 2 was replicated to give
that arm an error bar for Section 3.8, not to probe this trend, and its spread fell between its
neighbours in exactly the predicted order.

The mechanism is straightforward once seen: a trajectory leaving the site samples structures
that were never the complex, so its energy swings between seeds, while a tightly held pose
samples one basin and returns nearly the same number twice.

**But a stable pose is not a guarantee of a small spread.** TTBK2 pose 2 held in both replicates
— 1.99 and 1.89 Å, comfortably inside this study's stability band — and still returned a
3.97 kcal/mol spread against a reported SEM of 0.27. Pose stability bounds how bad the SEM can
be; it does not make the SEM usable.

**The rule this yields is sharper than "do not trust the SEM": the SEM is only as good as the
pose is stable, and you cannot tell which case you are in without running a replicate.** The MAO
magnitude in Section 3.8 is quotable because the replicate was run — not because its poses
looked stable. With five pairs across two protein families and three proteins this is a consistent
pattern, not a calibration curve, and it should not be used to predict a spread from an RMSD.

# 4. Discussion

## 4.1 What this study establishes, separately for each target pair

`cand_003` engages both on-targets, is drug-like and BBB-permeant on computed gates, and forms a
stable complex with MAO-B (2 of 3 poses hold; no MAO-B run left the site) and with TTBK1 in a
pose-specific way (pose 1 holds in both replicates; poses 2 and 3 drift).

Beyond that the two target pairs diverge, and **the central reporting decision of this paper is
to keep them apart:**

- **MAO-A/MAO-B: supportable in direction and magnitude.** Docking and MM-GBSA agree
  independently on the sign, and replicates put the margin at 2.55 kcal/mol against a measured
  spread of 0.24–0.26 — roughly a tenfold separation — with the protocol-bound caveat stated in
  Section 3.8.
- **TTBK1/TTBK2: not supportable.** The effect size is below the pipeline's own resolution, the
  margin's sign depends on the water shell, and the protocol is calibrated to be biased toward
  the anti-target by about the size of the effect. Any TTBK2 liability statement for this series
  rests on **docking scores alone**, and must say so.

An earlier pass of this work reported a TTBK2 liability as a finding. It does not survive, and
the correction is part of the result rather than an erratum.

## 4.2 Two methodological findings that bound studies of this design

**Docking rank does not predict pose stability.** Across three poses in each of four targets,
the top-ranked docked pose was the most stable one in one of four, and one target's best-holding
pose was its worst-ranked. A docking score ranks poses; it does not tell you whether a pose is a
physical minimum. Any pipeline that advances a single top-ranked pose into MD or free-energy
calculation is making an assumption this dataset contradicts three times out of four — and the
failure is silent, because an off-pose trajectory still produces a plausible-looking energy.

**A single-trajectory SEM understates the true uncertainty by 2- to 33-fold, in proportion to
pose instability.** This makes the SEM actively misleading in exactly the cases where a
selectivity claim is most tempting: the reported error is smallest-looking when the pose is
least stable relative to the truth. Quoting it on a ΔΔG is the easiest way to manufacture a
significant selectivity result from this class of pipeline. In this study, doing so would have
turned a 0.53 kcal/mol non-result into an apparent five-sigma finding. And a stable pose is not
a licence to skip the replicate: the TTBK2 pose-2 pair held at 1.99 and 1.89 Å and still spread
by 3.97 kcal/mol.

## 4.3 Comparison arms must be symmetric — four times over

Five separate false positives arose in this work, and all five had the same shape: an asymmetry
between the two things being compared.

1. A **0-versus-5 crystallographic water shell** between two receptors produced an apparent TTBK
   selectivity margin (Section 3.6).
2. A **pose scan in one arm against a replicate scan in the other** produced an apparent
   pose-stability difference between the kinases that vanished when both arms were scanned over
   three poses (Section 2.5).
3. An **off-pose trajectory against an on-pose one** inverted the MAO conclusion. MAO-A pose 1
   returns −44.22 kcal/mol — the most favourable MAO-A value in the dataset, by 7.55 kcal/mol,
   and the only one from a run that had drifted off its docked pose (3.24 Å). Taken at face value
   it says the compound prefers the anti-target: a selectivity *and* safety liability that does
   not exist. Every on-pose MAO-A trajectory lands near −36. **The value is named here once, so
   that a reader comparing against the supplementary table knows why it is excluded, and it is
   not a binding energy** — it describes a structure that was never docked.
4. A **cross-paper midpoint against a single-paper value** in the calibration pair itself. The
   two IC50s now come from one assay in one paper, which is what makes it a calibration rather
   than a comparison of two laboratories.
5. A **replicated arm against an unreplicated one.** Until 2026-09-28 the TTBK ΔΔG of
   1.45 kcal/mol was judged against TTBK1's replicate spread alone, because TTBK2 had a single
   on-pose trajectory. Replicating it put the ΔΔG at 0.53 with spreads of 2.08 and 3.97: the
   conclusion did not change, but the number it rested on did, by a factor of nearly three. The
   asymmetry here was in the *uncertainty*, not in the systems — which makes it the easiest of the
   five to miss, since both arms looked properly matched in design.

Checking arm symmetry before interpreting any comparison is, on this evidence, not a refinement
but a precondition. It is cheap, and each of these four cost days.

## 4.4 Limitations

1. **No experimental validation.** Every affinity reported is computed. No compound was
   synthesised; no biochemical or cellular assay was performed. The literature potencies are
   anchors for reference compounds, not measurements of `cand_003`.
2. **One ligand through the dynamics.** Every MD and MM-GBSA number is for `cand_003`. Nothing
   generalises to the series without more runs.
3. **10.1 ns per run.** Long enough to see a pose fail, too short for a binding/unbinding
   equilibrium or a slow conformational change.
4. **Single-trajectory MM-GBSA, no entropy term.** Interaction-energy estimates, not binding
   free energies; absolute values are not comparable to experiment, and **not comparable between
   proteins** (Section 2.7).
5. **Two runs per replicated arm, not more.** All four arms now have a velocity-replicate pair,
   which is enough to measure a spread — and to show the MAO ΔΔG clears its spread ~10-fold while
   the TTBK ΔΔG does not — but two runs do not give a distribution, and the spreads themselves
   are therefore point estimates.
6. **FAD is a restrained GAFF2 residue, not the covalent 8α-S-cysteinyl cofactor it is.** The
   flavin wall is reproduced and validated (0.46–0.73 Å across eight systems), but FAD cannot
   relax in response to the ligand, so induced fit involving the flavin is suppressed.
   Defensible for ligand MM-GBSA, where FAD is part of the receptor on both sides of the
   subtraction; not a substitute for covalent parameterisation. It is also why the MAO replicate
   spread is so tight.
7. **No independent scoring-function cross-check on the validated receptors** (Section 3.3).
   The Vinardo result stands as a general warning from this project's earlier data, not as
   validation of the rankings reported here.
8. **BBB assessment is a 2D property model.** It uses WLogP/TPSA and does not model efflux
   transporters, notably P-glycoprotein, to which flavonoids are known substrates.
9. **The score-versus-potency correlation is mixed-species**, and its two rat-brain points are
   the two flavonoids — the compounds closest to the lead (Section 2.2).
10. **The dual-target rationale is the study's hypothesis, not its finding.** No literature
    proposes or validates TTBK1 + MAO-B as a coherent strategy, and the two targets sit in
    different disease-mechanism families. A compound engaging both would need a therapeutic
    context in which simultaneous tau-kinase and MAO-B inhibition is desirable — plausible in
    principle for mixed-pathology neurodegeneration, not demonstrated.

## 4.5 What would change the answer

- **To settle the TTBK pair either way** needs much longer sampling or an alternative
  free-energy method, and the replicate of 2026-09-28 made that conclusion stronger rather than
  weaker. The pipeline's resolution floor is now ~4 kcal/mol — the larger of the two measured
  TTBK spreads — against an effect size of ~0.53, with the ~1.3 kcal/mol calibration floor
  sitting between them. More runs of the same protocol will not close a gap of that shape. The
  calibration floor is the harder of the two limits: it is set by a published IC50 pair, not by
  compute, so no amount of further docking or sampling would tighten it.
- **A third and fourth replicate per arm** would turn the spreads from point estimates into
  distributions. That is the only way to put a confidence interval on the error bars themselves,
  which is what a formal claim about the MAO ΔΔG would eventually need.
- **To generalise beyond `cand_003`**, run MD on the next two or three shortlisted candidates.
  Pose stability has been the discriminating filter at every step and is cheap relative to its
  value.
- **To close the scoring-function gap**, repeat the independent-function cross-check on the
  validated receptors.
- **Experimentally**, the informative first assay is not a single-target potency measurement but
  a side-by-side isoform panel — MAO-A/MAO-B and TTBK1/TTBK2 in parallel — since selectivity,
  not potency, is the open question. Covalent FAD parameterisation would be the corresponding
  computational upgrade.

# 5. Conclusions

A BBB-constrained generative campaign on a deoxygenated flavonol core produced a
4-trifluoromethyl/methyl-substituted 7-deoxyflavonol that engages both TTBK1 and MAO-B, on
receptors that all reproduce their crystallographic poses.

The selectivity question, which is the one that matters for this pair of targets, resolves
differently for each. For the MAO isoforms the intended direction holds and now carries a size:
MAO-B favoured by 2.55 kcal/mol against a measured replicate spread of 0.24–0.26, with docking
agreeing independently. For the TTBK paralogs there is no resolvable discrimination: with both arms
replicated the difference is 0.53 kcal/mol against measured spreads of 2.08 and 3.97, and saying
so required calibrating the protocol against a compound crystallised in both paralogs and
finding a ~1.0 kcal/mol bias toward the anti-target.

Two methodological results are reported as primary findings because they bound what any study of
this design may claim. Docking rank did not predict which pose survived dynamics in three of
four targets. And the single-trajectory standard error is not an error bar: across five replicate
pairs it understates the true run-to-run spread by 2- to 33-fold, rising monotonically with pose
instability, so it is least trustworthy exactly where it is most tempting to quote — and even a
pose that holds through both replicates can spread by 4 kcal/mol.

This work is hypothesis-generating. It identifies a chemically tractable dual-engagement
scaffold with a resolvable MAO-B preference, and identifies TTBK isoform selectivity as
unresolved at this level of theory — not as absent, but as beyond the reach of the method used
to look for it.

# Data and code availability

Every number in this manuscript is traceable to a file in the project workspace, and every
derived table is generated by a script rather than transcribed. The MD and MM-GBSA summary table
is regenerated from the primary `MMPBSA.py` and `cpptraj` output by
`scripts/collect_md_summary.py`, which fails if the table is stale; all three figures are drawn
from that same table by `10_results/make_figures.py`, so a figure cannot disagree with the text
and neither can disagree with what was computed. Docking consensus, selectivity margins,
filtering cascade and the 9IV calibration are produced by `scripts/collect_results.py`,
`scripts/analyze_selectivity.py`, `scripts/filter_cascade.py` and `scripts/calibration_9iv.py`
respectively.

Available: the flavonoid library with provenance and confidence tiers; the full generative-run
output; filtering-cascade tables; prepared receptors and ligands; all docking runs including the
wet/dry, no-cofactor and bridging-water controls; redocking validation logs; per-system MD
trajectories, RMSD/RMSF analyses and raw `MMPBSA.py` output for all seventeen systems; and a
chronological record of every decision and every defect found and corrected.

# References

1. Halkina T, et al. Discovery of potent and brain-penetrant tau tubulin kinase 1 (TTBK1)
   inhibitors that lower tau phosphorylation in vivo. *J Med Chem*. 2021;64(9):6358–6380.
   PMID 33944571.
2. Nozal V, et al. Tau Tubulin Kinase 1 (TTBK1), a new player in the fight against
   neurodegenerative diseases. *J Med Chem*. 2022;65(2):1585–1607. PMID 34978799.
   (Compound 42 / VNG2.73 = PDB ligand 9IV; TTBK1 IC50 330 nM, TTBK2 IC50 490 nM, same assay.
   Source of the calibration pair in Section 3.5.)
3. Binda C, et al. Structures of human monoamine oxidase B complexes with selective noncovalent
   inhibitors: safinamide and coumarin analogs. *J Med Chem*. 2007;50(23):5848–5852.
   PMID 17915852. (PDB 2V5Z.)
4. Stössel A, et al. Development of molecular probes for the imaging of monoamine oxidase.
   *J Med Chem*. 2013;56(11):4580–4596. PMID 23631427. (Safinamide IC50 7.67 nM, recombinant
   human MAO-B, p-tyramine substrate.)
5. Son SY, et al. Structure of human monoamine oxidase A at 2.2 Å resolution: the control of
   opening the entry for substrates/inhibitors. *Proc Natl Acad Sci USA*. 2008;105:5739–5744.
   (PDB 2Z5X.)
6. Xue F, et al. Discovery of macrocyclic inhibitors of tau tubulin kinase.
   *ChemMedChem*. 2013;8(11):1846–1854. PMID 24039150. (PDB 4BTK ligand DTQ; *K*d 240 nM by
   surface plasmon resonance, distinct from its 4610 nM enzymatic IC50.)
7. Daina A, Zoete V. A BOILED-Egg to predict gastrointestinal absorption and brain penetration
   of small molecules. *ChemMedChem*. 2016;11(11):1117–1121.
8. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new
   scoring function, efficient optimization, and multithreading. *J Comput Chem*.
   2010;31(2):455–461.
9. Quiroga R, Villarreal MA. Vinardo: a scoring function based on AutoDock Vina improves
   scoring, docking, and virtual screening. *PLoS One*. 2016;11(5):e0155183.
10. Loeffler HH, He J, Tibo A, et al. REINVENT 4: modern AI-driven generative molecule design.
    *J Cheminform*. 2024;16:20.
11. Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: improving
    the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory
    Comput*. 2015;11(8):3696–3713.
12. Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general
    Amber force field. *J Comput Chem*. 2004;25(9):1157–1174. (GAFF; GAFF2 as distributed with
    AmberTools.)
13. Sousa da Silva AW, Vranken WF. ACPYPE — AnteChamber PYthon Parser interfacE.
    *BMC Res Notes*. 2012;5:367.
14. Eastman P, et al. OpenMM 7: rapid development of high performance algorithms for molecular
    dynamics. *PLoS Comput Biol*. 2017;13(7):e1005659.
15. Miller BR III, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE. MMPBSA.py: an
    efficient program for end-state free energy calculations. *J Chem Theory Comput*.
    2012;8(9):3314–3321.
16. Roe DR, Cheatham TE III. PTRAJ and CPPTRAJ: software for processing and analysis of
    molecular dynamics trajectory data. *J Chem Theory Comput*. 2013;9(7):3084–3095.
17. Meli R, Biggin PC. spyrmsd: symmetry-corrected RMSD calculations in Python.
    *J Cheminform*. 2020;12:49.
18. Bouysset C, Fiorucci S. ProLIF: a library to encode molecular interactions as fingerprints.
    *J Cheminform*. 2021;13:72.
19. Chaurasiya ND, et al. Selective inhibition of human monoamine oxidase B by O-methylated
    flavonoids. *Molecules*. 2020;25(22):5358.
