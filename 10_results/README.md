# Step 10 — Results and Interpretation

## What we did

Collected every quantitative result the pipeline produced into one place and asked the only
question that matters for the project: **what can this study actually claim about `cand_003`, and
with what confidence?**

Nothing new was computed here. Every number is read from the primary files in steps 4–7 — mostly
from [`../07_mmgbsa/md_mmgbsa_summary.csv`](../07_mmgbsa/md_mmgbsa_summary.csv), which is itself
regenerated from the `.dat` files by
[`../scripts/collect_md_summary.py`](../scripts/collect_md_summary.py). The figures here are drawn
by [`make_figures.py`](make_figures.py) from that same table, so a figure cannot disagree with the
text, and neither can disagree with what was computed.

## Why we did it

The step READMEs each answer *what happened at that step*. None of them answers **what the study
concludes**, and that gap is where a project like this goes wrong: a reader assembles the headline
from whichever numbers are largest, rather than from the numbers that are valid. This project has
already produced four documented instances of exactly that failure, all with the same shape — a
number that looked like a result until someone checked whether the two things being compared were
comparable.

So this directory is organised by **claim**, not by method, and every claim carries its status.

---

## The project's question

A dual-target design: one molecule intended to inhibit **TTBK1** (a tau kinase implicated in
Alzheimer's disease) and **MAO-B**, while *avoiding* their close relatives.

| Protein | Role | Why the distinction matters |
|---|---|---|
| **TTBK1** (7JXX) | intended target | tau phosphorylation |
| TTBK2 (7Q8Y) | **anti-target** | loss-of-function mutations cause spinocerebellar ataxia type 11, so unintended TTBK2 inhibition is a specific safety liability |
| **MAO-B** (2V5Z) | intended target | the neurodegeneration-relevant isoform |
| MAO-A (2Z5X) | **anti-target** | MAO-A inhibition without MAO-B selectivity raises the tyramine pressor risk ("cheese effect") |

Both safety rationales are the manuscript's own (§4). They are what make *selectivity*, not
potency, the load-bearing claim of the study — and therefore what the MD and MM-GBSA work was
built to test.

### The lead compound

`cand_003` — `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`, a trifluoromethyl/methyl-substituted
flavonol from the REINVENT4 campaign.

| MW | WLogP | TPSA | HBD | HBA | Lipinski | BBB | GI absorption | Structural alerts |
|---|---|---|---|---|---|---|---|---|
| 336.27 | 4.20 | 70.67 Å² | 2 | 4 | 0 violations | pass | pass | 0 |

Docking consensus over 3 seeds on the validated receptors: **−8.39 kcal/mol on TTBK1 (7JXX),
rank 8 of 56** and **−11.41 on MAO-B (2V5Z), rank 2 of 56**. It is drug-like and BBB-permeant on
the computed gates, which for a CNS target is the precondition for anything else mattering.

---

## Result 1 — Most docked poses do not survive dynamics

![Pose stability across all fourteen MD systems](fig1_pose_stability.png)

Sixteen systems, 10.1 ns each, ligand RMSD from the docked pose averaged over the final 100
frames. Stability is judged on that window rather than the whole run, because a ligand that leaves
late still shows a low whole-run average.

**Eight of sixteen runs hold their pose; two leave the site entirely.** The immediate consequence
is that a docking score is not evidence of binding: it ranks poses, and the ranking does not
predict which pose is physically stable.

**Docking rank predicted the most stable pose in one of four targets.** Only TTBK1's top-ranked
pose won. MAO-A's best-holding pose is its *worst*-ranked of three (pose 3, Vina −7.40 against
pose 1's −8.24). With 3 poses × 4 targets this is the most strongly supported methodological
finding in the project, and it is why every energy below is tied to a named pose.

**The MAO complexes are the more stable pair.** No MAO run left the site — the worst is 3.24 Å,
against two TTBK runs past 5 Å. Counting strictly, MAO-A holds 1 of 3 and MAO-B 2 of 3.

Protein backbones were stable throughout (1.12–1.81 Å; 1.06–1.77 Å core-fit for MAO), so the
ligand motion is ligand motion, not a collapsing binding site. (Core-fit RMSD across the eight
MAO runs is 1.06–1.77 Å.) One apparent exception —
`system_MAOA_p2` at 3.44 Å, the highest in the project — is the solvent-exposed C-terminal tail of
the MAO-A construct flailing without a membrane to sit in, not an unstable fold; per-residue RMSF
puts the maximum at residues 510–513 of 513, and the core is 1.28 Å.

## Result 2 — For MAO, the intended selectivity is supported in direction

![MAO binding free energies, MAO-A against MAO-B](fig2_mao_binding_energy.png)

This is the study's central new result, and it reversed itself once the poses were done properly.

On-pose, MAO-A lands at −35.54 / −36.67 / −36.40 and MAO-B at −39.48 / −39.20 / −38.97 /
−39.04. **All twelve on-pose pairings favour MAO-B, by 2.30 to 3.94 kcal/mol** — the intended
direction, since MAO-A is the anti-target. Docking independently favours MAO-B by 3.03–3.22
kcal/mol pose for pose. Two methods that disagreed now agree.

**And the size of the margin now has a measured error bar.** Velocity replicates of both
best-pose systems — identical topology and coordinates, differing only in the random velocity
seed — give run-to-run spreads of **0.26** (MAO-A pose 3) and **0.24** (MAO-B pose 2) kcal/mol.
Comparing replicate means, MAO-A −36.53 against MAO-B −39.09 is a **ΔΔG of 2.55 kcal/mol,
9.8× the larger spread.** Both replicates also reproduced the pose stability (1.07 and 1.37 Å).

The caveat that must travel with that number: FAD is positionally restrained in both arms, which
suppresses receptor motion and so damps the run-to-run variation. 0.24–0.26 is the right error bar
*for this protocol*, not evidence that MM-GBSA is intrinsically this precise. The ΔΔG is fair
because **both arms carry the same restraint.**

**What made them disagree was one off-pose trajectory.** MAO-A pose 1 gives −44.22, the most
favourable MAO-A value by 7.55 kcal/mol and the only one from a run that had drifted off its docked
pose (3.24 Å). Taken at face value it says the compound prefers the anti-target — a selectivity
*and* safety liability that does not exist. Every on-pose MAO-A trajectory lands near −36.

This is the clearest case in the project of an off-pose trajectory producing not a noisy number but
an **inverted conclusion**, and it is the strongest argument for the rule the project already
follows: never compute a binding energy on a trajectory that has left its docked pose.

## Result 3 — For TTBK, there is no resolvable discrimination

The best on-pose comparison is TTBK1 −31.36 / −33.44 against TTBK2 −30.95: a ΔΔG of about
1.45 kcal/mol, *smaller* than the 2.08 kcal/mol spread between two TTBK1 runs that differed only
in their random velocity seed. Pose stability says the same thing — 1 of 3 poses hold on each
kinase, under a symmetric test.

Three independent lines say the original TTBK2 liability claim does not survive:

- The +0.12 kcal/mol docking margin came from a receptor pair with **five waters in one site and
  none in the other**; run water-symmetrically, the margin *reverses* to favour TTBK2 by ~1.6.
- The protocol carries a **1.018 kcal/mol systematic bias toward TTBK2** on the 9IV pair,
  against a same-assay experimental ΔΔG of −0.234 kcal/mol. Any margin below ~1.0 is inside
  the bias.
- MD and MM-GBSA find no resolvable difference.

**The honest outcome is that dynamics lacks the resolution to speak on TTBK1 vs TTBK2 at all.**
That is not a null result to be buried — a 1.45 kcal/mol margin quoted against the reported SEM
would have looked like a five-sigma selectivity finding.

## Result 4 — The reported SEM is not the error bar

![SEM against the spread measured by re-running with a new velocity seed](fig3_sem_vs_replicate.png)

`MMPBSA.py` computes its standard error as if ~200 frames sampled 10 ps apart were independent
draws. They are not. Across the project's four replicate pairs the reported SEM stays in a narrow
0.15–0.27 kcal/mol band while the measured spread between two runs differing *only* in velocity
seed ranges from 0.24 to 6.89 — a **29-fold** range.

**And the discrepancy is not random: it tracks pose stability, monotonically across all four
pairs.**

| Pair | Mean ligand RMSD | Measured spread | Ratio to SEM |
|---|---|---|---|
| MAO-B pose 2 | 1.23 Å | 0.24 | 2× |
| MAO-A pose 3 | 1.25 Å | 0.26 | 2× |
| TTBK1 pose 1 | 1.63 Å | 2.08 | 8× |
| TTBK2 pose 1 | 4.70 Å | 6.89 | 33× |

The mechanism is straightforward once seen: a trajectory leaving the site samples structures that
were never the complex, so its energy swings between seeds; a tightly held pose samples one basin
and returns nearly the same number twice.

**Quoting the SEM as the uncertainty on a ΔΔG is the single easiest way to manufacture a
significant selectivity result from this pipeline.** Without the second replicate, that is exactly
what would have happened for TTBK.

The rule this yields is sharper than "never trust the SEM": **the SEM is only as good as the pose
is stable, and you cannot tell which case you are in without running a replicate.** The MAO
magnitude is quotable because the replicate was run — not because its poses looked stable. With
n = 4 pairs across two protein families this is a consistent pattern, not a calibration curve.

---

## What the project can and cannot claim

This is the section to write the thesis and the manuscript Results from.

| # | Claim | Status | Rests on |
|---|---|---|---|
| 1 | `cand_003` is drug-like and BBB-permeant on computed gates | **Supportable** | 0 Lipinski violations, BBB + GI pass, 0 structural alerts |
| 2 | Docking rank does not predict which pose survives dynamics | **Supportable, strongly** | 1 of 4 targets, 3 poses each, 14 runs |
| 3 | `cand_003` forms a stable complex with MAO-B | **Supportable** | 2 of 3 poses hold (1.10, 2.28 Å); no MAO-B run left the site |
| 4 | `cand_003` forms a stable complex with TTBK1 | **Supportable, pose-specific** | pose 1 holds in both replicates (1.88, 1.38 Å); poses 2–3 drift |
| 5 | The compound favours MAO-B over MAO-A — the intended direction | **Supportable** | all 12 on-pose pairings agree; docking agrees independently |
| 6 | …by about 2.55 kcal/mol | **Supportable, protocol-bound** | 9.8× the measured 0.24–0.26 replicate spread; both arms share the FAD restraint that damps it |
| 7 | The compound is selective for TTBK1 over TTBK2 | **Not supportable** | 1.45 < 2.08 replicate spread; no pose-stability difference |
| 8 | TTBK2 is an off-target liability for this series | **Docking only** | reverses under a symmetric water shell; inside the 1.018 protocol bias |
| 9 | The reported SEM is not a usable error bar, and its failure scales with pose instability | **Supportable, strongly** | 4 replicate pairs, ratio 2× to 33×, monotonic in ligand RMSD |
| 10 | Absolute ΔG values are comparable between targets | **Not supportable** | protein-specific desolvation/surface terms do not cancel |

### Two things not to write

**Do not quote the −44.22.** It is the most attractive number in the dataset and it describes a
structure that was never docked. It is kept on record only as the cautionary case.

**Do not put TTBK and MAO binding energies on the same axis.** They carry protein-specific terms
that do not cancel; MAO-A's −44.22 against TTBK1's −31.36 compares a flavoenzyme to a kinase and
means nothing. Figure 2 is MAO-only for this reason, and figure 1 shares an axis legitimately
because ligand RMSD is a geometric measure of one ligand against its own starting pose.

---

## Limitations

Stated plainly, because each one bounds a claim above.

1. **One ligand.** Every MD and MM-GBSA number is for `cand_003`. Nothing here generalises to the
   series without more runs.
2. **10.1 ns per run.** Long enough to see a pose fail, too short for a binding/unbinding
   equilibrium or a slow conformational change.
3. **Single-trajectory MM-GBSA, no entropy term.** These are interaction-energy estimates, not
   binding free energies in the thermodynamic sense, and their absolute values are not comparable
   to experiment.
4. **Two replicates per MAO arm, not more.** Enough to measure a spread and show the ΔΔG clears
   it ~10-fold; not enough for a distribution. The TTBK arms still have no on-pose replicate pair
   beyond TTBK1 pose 1.
5. **FAD is a restrained GAFF2 residue, not the covalent 8α-S-cysteinyl cofactor it really is.**
   Beyond the induced-fit point below, this is also why the MAO replicate spread is so tight — the
   restraint removes receptor motion that would otherwise vary between seeds.
   The flavin wall of the cavity is reproduced and validated (0.46–0.73 Å from crystal across six
   systems), but FAD cannot relax in response to the ligand, so induced fit involving the flavin is
   suppressed. Defensible for ligand MM-GBSA, where FAD is part of the receptor on both sides of
   the subtraction; not a substitute for covalent parameterisation.
6. **No experimental validation.** Everything here is computational. The IC50 values in
   [`../01_smiles/`](../01_smiles/) are literature anchors for reference compounds, not measurements
   of `cand_003`.

## What would change the answer

- **Done 2026-09-27:** velocity replicates of the two best on-pose MAO systems
  (`system_MAOA_p3_r2`, `system_MAOB_p2_r2`) via
  [`../scripts/run_mao_replicate_queue.sh`](../scripts/run_mao_replicate_queue.sh). They measured
  the MAO arms' own spread — 0.26 and 0.24 kcal/mol — and moved claim 6 from unsupportable to
  supportable.
- **The obvious next replicate is TTBK2 pose 2**, the one on-pose TTBK2 trajectory. Its ΔG
  (−30.95) currently carries no error bar of its own, and it is the number the TTBK1-vs-TTBK2
  comparison rests on.
- **To support claim 7 either way**, the TTBK arm needs either much longer sampling or an
  alternative free-energy method. The current pipeline's resolution floor (~2 kcal/mol) is above
  the effect size (~1.45), so more of the same will not settle it.
- **To generalise beyond `cand_003`**, MD the next two or three shortlisted candidates. Pose
  stability has been the discriminating filter at every step and is cheap relative to its value.

## Relevant files

| File | Role |
|---|---|
| [`make_figures.py`](make_figures.py) | draws all three figures from the collected table; rerun after any new MD |
| [`fig1_pose_stability.png`](fig1_pose_stability.png) | result 1 — which poses survive |
| [`fig2_mao_binding_energy.png`](fig2_mao_binding_energy.png) | result 2 — the MAO comparison |
| [`fig3_sem_vs_replicate.png`](fig3_sem_vs_replicate.png) | result 4 — the error-bar finding |
| [`../07_mmgbsa/md_mmgbsa_summary.csv`](../07_mmgbsa/md_mmgbsa_summary.csv) | **the source table** — every RMSD and energy quoted here |
| [`../scripts/collect_md_summary.py`](../scripts/collect_md_summary.py) | regenerates that table from primary `.dat` files; `--check` fails if stale |
| [`../08_analysis/selectivity_margins.csv`](../08_analysis/selectivity_margins.csv) | per-candidate docking margins against each anti-target |
| [`../08_analysis/filter_cascade_candidates_56.csv`](../08_analysis/filter_cascade_candidates_56.csv) | the ADMET/BBB properties in the lead table above |
| [`../LOGBOOK.md`](../LOGBOOK.md) | the chronological record, and the authoritative account of *why* — including every defect found |
| [`../09_manuscript/README.md`](../09_manuscript/README.md) | which manuscript claims need revising, section by section |

Where this file and the logbook disagree, **the logbook is the primary source.**
