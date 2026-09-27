# Step 9 — Manuscript and Reporting

## What we did

Drafted the study write-up, a figure/table index, and a target-journal assessment.

| File | Contents |
|---|---|
| [`manuscript_dualtarget_flavonol.md`](manuscript_dualtarget_flavonol.md) | the manuscript source — Abstract, Introduction, Methods (2.1–2.8), Results (3.1–3.7), Discussion |
| `manuscript_dualtarget_flavonol.html` / `.docx` | rendered outputs |
| [`manuscript_figures_tables_and_journals.md`](manuscript_figures_tables_and_journals.md) | figure/table index, supplementary data list, journal recommendations, title options |
| [`Week1_Progress_Report.md`](Week1_Progress_Report.md) | Week 1 progress report |

## Why we did it

The manuscript is the deliverable, but it is also the place where this project's methodological
findings either get stated honestly or get quietly smoothed over. The Results sections are written
as claims, and each one has to be traceable to a specific directory in this workspace.

## Reference

Methods sections map onto the workflow steps directly:

| Manuscript section | Workflow step |
|---|---|
| 2.1 Species selection and natural-product library | [`00_library/`](../00_library/) |
| 2.2 Receptor preparation and redocking validation | [`03_receptors/`](../03_receptors/), [`05_validation/`](../05_validation/) |
| 2.3 Generative design under a BBB-aware reward | [`generation/`](../generation/), [`00_library/reinvent4_output/`](../00_library/) |
| 2.4 Filtering cascade | [`08_analysis/`](../08_analysis/), [`filtering/`](../filtering/) |
| 2.5 Consensus re-docking and scoring cross-check | [`04_docking/`](../04_docking/), [`08_analysis/`](../08_analysis/) |
| 2.6 Anti-target counter-screening | [`08_analysis/`](../08_analysis/) |
| 2.7 Molecular dynamics and MM-GBSA | [`06_md/`](../06_md/), [`07_mmgbsa/`](../07_mmgbsa/) |

Target journals, in the document's own order of preference: *Journal of Cheminformatics*
(recommended primary), *Molecules*, *Frontiers in Chemistry* / *Molecular Biosciences*, *IJMSci*
as fallback.

## ⚠ Claims requiring revision before submission

**The manuscript predates the Week 1–2 MD work and three of its Results claims are now out of
date.** This is the most important thing in this directory.

### §3.5 — "Both complexes are stable over 20 ns of molecular dynamics"

**Must be rewritten.** That claim rests on the prior-phase MD, run on a different (apo) structure.
The current, structurally validated MD ([`06_md/`](../06_md/)) shows stability is **pose-dependent**:

- TTBK1 holds the ligand in **1 of 3** docked poses
- TTBK2 holds it in **1 of 3** docked poses
- Runs from the other TTBK poses drift (3.4–4.2 Å) or dissociate outright (5.7–5.8 Å)
- MAO-A holds it in **1 of 3** strictly (pose 3, 1.43 Å), 2 of 3 counting pose 2 at 2.43 Å
- MAO-B holds it in **2 of 3** strictly (poses 1 and 2, 2.28 and 1.10 Å); no MAO run dissociated

"Both complexes are stable" is true only of specific poses, and which pose that is differs between
proteins. The honest version of this claim is narrower and more interesting. Note that it is
closest to true for the MAO pair — all six MAO runs stayed in the site, the worst at 3.24 Å,
whereas two of the eight TTBK runs passed 5 Å — and least true for TTBK2, where two of three
poses dissociated outright.

### §3.6 — "MM-GBSA confirms engagement of both targets and prefers MAO-B"

**Directionally supportable as of 2026-09-27; the magnitude is not.** This section was previously
flagged because MM-GBSA appeared to point the *opposite* way to docking (MAO-A −44.22 against
MAO-B −39.48). The MAO pose scan showed that inversion was an artifact: the −44.22 came from the
one MAO-A trajectory that had drifted 3.24 Å off its docked pose, and it is the only off-pose
MAO-A number. On-pose, MAO-A lands at −35.54/−36.67 and MAO-B at −39.48/−39.20/−39.04, so
**all six on-pose pairings favour MAO-B by 2.37–3.94 kcal/mol**, agreeing in sign with docking,
which favours MAO-B by 3.03–3.22 kcal/mol pose for pose on this ligand's Vina affinities. The
library-wide docking margin quoted elsewhere (−3.28, favourable 48/48) is a mean over all 56
candidates and a different quantity; do not present it as the `cand_003` per-pose figure.

So "prefers MAO-B" may now be stated as a direction that docking and MM-GBSA agree on, with two
constraints that must travel with it:

- **No magnitude.** Neither MAO arm has a velocity replicate, and the margin is comparable to the
  2.08 kcal/mol replicate spread measured on TTBK1. Quote no number and no significance.
- **Still a cross-protein comparison.** Absolute MM-GBSA values carry protein-specific desolvation
  and surface terms that do not cancel, which is why the direction is quotable and the size is not.

See [`07_mmgbsa/README.md`](../07_mmgbsa/README.md). **Any −44.22 figure must not be cited as a
MAO-A binding energy** — it describes a structure that was never docked.

### §3.7 — "The series does not discriminate either target from its paralog/isoform"

**This claim now has stronger support than when it was written**, and should be strengthened rather
than softened. Four independent lines converge on it:

| Evidence | Result |
|---|---|
| Docking, water-symmetric | TTBK2 favoured ~1.6 kcal/mol; 0/56 favourable to TTBK1 |
| Protocol calibration (9IV) | 0.861 kcal/mol systematic bias toward TTBK2 where experiment says none |
| MD pose stability | no difference — 1/3 poses each under a symmetric test |
| MM-GBSA | no resolvable difference (1.45 < 2.08 replicate spread) |

Any TTBK2 "liability" statement must be attributed to **docking scores alone**, with the explicit
note that dynamics lacks the resolution to confirm or refute it.

**This claim must now be split by pair, because the two halves no longer agree.** The four lines
above all concern the TTBK1/TTBK2 paralog pair, where "does not discriminate" still holds. For the
MAO-A/MAO-B isoform pair it no longer reads the same way: after the pose scan, docking and MM-GBSA
agree in sign on a MAO-B preference across all six on-pose pairings (§3.6). The margin is still
not *resolvable* — it is comparable to the TTBK1 replicate spread and the MAO arms have no
replicates — so "no significant discrimination" survives for MAO, but "no discrimination" and any
wording implying the isoform comparison came out symmetric do not. Write the two pairs separately.

### Two methods points that belong in §2.8

1. **The `MMPBSA.py` SEM is not a valid error bar** — it treats 200 frames 10 ps apart as
   independent, giving 0.15–0.27 kcal/mol where the true replicate spread is 2.08–6.90. Report the
   replicate spread.
2. **Comparison arms must be symmetric.** Three separate false positives in this project all came
   from asymmetry between the two things compared (a 0-vs-5 water shell; a pose scan against a
   replicate scan; an on-pose run against an off-pose run).

### One reference to verify

Several `source` fields in [`../01_smiles/references.csv`](../01_smiles/references.csv) are marked
*"value as cited in project manuscript Introduction — CONFIRM primary source before submission."*
Also note that **safinamide is not the same molecule as the PDB ligand `SAG`** in 2V5Z.

## Final result

A complete draft covering the full pipeline, with a figure/table index and journal shortlist. Three
Results claims need revision against the Week 1–2 data before submission, as detailed above. The
central negative result (§3.7) is the most defensible claim in the paper.

## Relevant files

| Path | Role |
|---|---|
| [`manuscript_dualtarget_flavonol.md`](manuscript_dualtarget_flavonol.md) | **edit this**; HTML/DOCX are rendered from it |
| [`manuscript_figures_tables_and_journals.md`](manuscript_figures_tables_and_journals.md) | figure/table index and journal targets |
| [`Week1_Progress_Report.md`](Week1_Progress_Report.md) | Week 1 progress report |
| [`../LOGBOOK.md`](../LOGBOOK.md) | the authoritative record of what was done and why — **primary source for Methods** |
| [`../07_mmgbsa/md_mmgbsa_summary.csv`](../07_mmgbsa/md_mmgbsa_summary.csv) | all MD/MM-GBSA numbers in one table, for Results |
