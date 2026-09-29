# Step 9 — Manuscript and Reporting

## What we did

Wrote the study up from the Week 1–2 dataset, with a companion figure/table index and a
target-journal assessment.

| File | Contents |
|---|---|
| [`manuscript_dualtarget_flavonol.md`](manuscript_dualtarget_flavonol.md) | the manuscript source — Abstract, Introduction, Methods (2.1–2.7), Results (3.1–3.9), Discussion, Conclusions, References |
| [`manuscript_figures_tables_and_journals.md`](manuscript_figures_tables_and_journals.md) | figure/table index, supplementary data list, journal assessment, title options, framing requirements |
| [`prior_phase/`](prior_phase/) | the superseded August draft and its figure index, plus the Week 1 progress report |

**Edit the `.md`.** There are no rendered outputs at present: the August `.html` and `.docx` were
deleted on 2026-09-27 rather than carried forward, because they were generated from the
superseded draft and a stale render that looks authoritative is worse than no render. They remain
in git history (commit `e442b33`) if ever needed. Re-render from the current `.md` when the text
is final.

## Why we did it

The manuscript is the deliverable, but it is also where this project's methodological findings
either get stated honestly or quietly smoothed over. The Results sections are written as claims,
and each is traceable to a specific directory in this workspace.

## Reference

Methods sections map onto the workflow steps directly:

| Manuscript section | Workflow step |
|---|---|
| 2.1 Species selection, library and generative design | [`00_library/`](../00_library/), [`generation/`](../generation/) |
| 2.2 Receptors, and the reference set with its assays | [`03_receptors/`](../03_receptors/), [`01_smiles/`](../01_smiles/) |
| 2.3 Protocol validation and scoring-function calibration | [`05_validation/`](../05_validation/) |
| 2.4 Docking | [`04_docking/`](../04_docking/), [`08_analysis/`](../08_analysis/) |
| 2.5 Molecular dynamics | [`06_md/`](../06_md/) |
| 2.6 MM-GBSA | [`07_mmgbsa/`](../07_mmgbsa/) |
| 2.7 What this protocol cannot do | [`10_results/`](../10_results/) |
| 3.1–3.2 BBB ceiling, filtering cascade | [`08_analysis/`](../08_analysis/), [`filtering/`](../filtering/) |
| 3.3–3.6 Lead, validation, calibration, water shell | [`05_validation/`](../05_validation/), [`04_docking/`](../04_docking/) |
| 3.7–3.9 Pose stability, ΔΔG, the SEM finding | [`06_md/`](../06_md/), [`07_mmgbsa/`](../07_mmgbsa/), [`10_results/`](../10_results/) |

## The rewrite of 2026-09-27, and what changed

The August draft predated the Week 1–2 rebuild and **three of its Results claims did not
survive**. It was archived to [`prior_phase/`](prior_phase/) rather than edited, because the
receptor set, the MD protocol and the conclusions all changed together and a patched draft would
have mixed two incompatible passes.

| August draft | Current manuscript |
|---|---|
| Receptors `4NFM` (apo TTBK1), `6U0K`, `2V60`; one passing redocking validation | `7JXX`/`7Q8V`/`7Q8Y`/`2V5Z`/`2Z5X`/`4BTK`, **all six pass** (0.27–1.45 Å) |
| §3.5 "Both complexes are stable over 20 ns" | §3.7 stability is **pose-dependent**: 8 of 16 runs hold, 2 dissociate; docking rank predicted the most stable pose in **1 of 4** targets |
| §3.6 MM-GBSA −26.32 ± 0.53 / −36.85 ± 0.23, quoted as SEM | §3.8 on-pose only, with **replicate spread** as the error bar; the SEM is shown to be invalid in §3.9 |
| §3.7 one combined "no isoform discrimination" claim | §3.8 **split by pair** — MAO supportable in direction *and* size (2.55 vs 0.24–0.26); TTBK not supportable at all |
| No protocol calibration | §3.5 a measured **~1.0 ± 0.25 kcal/mol bias toward TTBK2**, giving a ~1.3 kcal/mol interpretability floor |
| Selectivity margins from unvalidated anti-targets, asymmetric water shells | §3.6 water-symmetric throughout, with the asymmetry itself reported as a false-positive source |
| Limitation "single 20 ns replicate, no replicate MD" | **Deleted, not softened** — the replicates were run |
| 15-candidate shortlist | 56 candidates; MAO margin over the **48 quantifiable** ones, 8 sterically excluded |

Two findings are new and are now reported as primary results: **docking rank does not predict
pose stability** (§3.7) and **the single-trajectory SEM understates the true spread by 2–33×, in
proportion to pose instability** (§3.9).

## Standing rules for any edit

These are not style preferences; each one exists because violating it produced a wrong result in
this project. The full account is in [`../LOGBOOK.md`](../LOGBOOK.md).

1. **Never quote the −44.22** as a MAO-A binding energy. It is the only off-pose MAO-A value and
   it inverted the MAO conclusion for two days. It appears once, in §4.3, as the cautionary case.
2. **Never put TTBK and MAO binding energies on the same axis or in the same comparison.**
   Protein-specific desolvation and surface terms do not cancel.
3. **Report the replicate spread, never the `MMPBSA.py` SEM**, as the uncertainty on a ΔΔG.
4. **Quote the calibration bias as ~1.0, never 1.018.**
5. **Keep the two target pairs separate.** They reach opposite conclusions.
6. **Attribute any TTBK2 liability claim to docking scores alone.**
7. **Label the score-versus-potency correlation as mixed-species** (two rat-brain points, and
   they are the two flavonoids).
8. **Check arm symmetry before interpreting any comparison.** Four false positives in this
   project all came from an asymmetry — water shell, pose-scan versus replicate-scan, on-pose
   versus off-pose, and cross-paper versus single-paper IC50.

## Final result

A complete manuscript written from the current dataset, with the two target pairs reported
separately and the two methodological findings foregrounded. The most defensible claims in the
paper are the methodological ones; the most defensible chemical claim is the MAO-B preference,
in both direction and magnitude.

**Known gap:** no independent scoring-function cross-check exists on the validated receptors. The
Vinardo result quoted in §3.3 was run on the earlier receptor set and is presented as a general
warning, not as validation of the current rankings. Closing that gap is listed in
[`../WORKFLOW.md`](../WORKFLOW.md) Phase 3.

## Relevant files

| Path | Role |
|---|---|
| [`manuscript_dualtarget_flavonol.md`](manuscript_dualtarget_flavonol.md) | **edit this** |
| [`manuscript_figures_tables_and_journals.md`](manuscript_figures_tables_and_journals.md) | figures, tables, supplementary index, journals |
| [`prior_phase/`](prior_phase/) | the superseded August draft, its figure index, and the Week 1 progress report |
| [`../10_results/README.md`](../10_results/README.md) | the claims table with each claim's status — **write Results from this** |
| [`../07_mmgbsa/md_mmgbsa_summary.csv`](../07_mmgbsa/md_mmgbsa_summary.csv) | every MD/MM-GBSA number in one generated table |
| [`../LOGBOOK.md`](../LOGBOOK.md) | the authoritative record of what was done and why — **primary source for Methods** |
