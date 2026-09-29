# Step 8 — Consensus Ranking, Filtering and Selectivity

## What we did

Turned raw docking output into decisions: consensus scores across seeds, a drug-likeness and
BBB filtering cascade, selectivity margins between each target and its anti-target, and a
dedicated control quantifying how much of the selectivity signal is a crystallographic-water
artifact.

## Why we did it

Three questions have to be answered before a candidate can be called a lead, and none of them is
answerable from a single docking score:

1. **Is the score real or is it noise?** Vina's search is stochastic, so scores are consensus
   across three seeds with a reported SD.
2. **Could the molecule ever be a CNS drug?** Both target families are intracranial. A potent
   compound that cannot cross the blood–brain barrier is not a candidate, so the cascade gates on
   Lipinski, GI absorption, BBB permeability and structural alerts *before* potency is considered.
3. **Is the selectivity margin a property of the chemistry or of the file?** This is the question
   the water test exists to answer, and the answer turned out to be "largely the file".

## Reference

- **Consensus** — mean Vina score across seeds 11/22/33, with per-ligand SD
- **Independent scoring cross-check** — Vinardo, a different scoring function on the same poses
- **Drug-likeness** — Lipinski rule-of-five
- **BBB / absorption** — BOILED-Egg (WLogP vs TPSA), implemented in
  [`../filtering/boiled_egg_coords.py`](../filtering/) and
  [`../scripts/bbb_score.py`](../scripts/bbb_score.py)
- **Structural alerts** — PAINS and reactive-group alerts
- **Descriptors** — RDKit (MW, WLogP, TPSA, HBD, HBA)

Prior-phase reports in [`prior_phase/`](prior_phase/) carry the full method text for the earlier
pass: `phase4_report.md` (filtering), `phase5_shortlist_report.md` (novelty and shortlist
selection), `phase6_consensus_report.md` and `phase6_final_report.md` (consensus re-docking plus
independent scoring validation), `phase6b_selectivity_report.md` (anti-target counter-screening)
and `phase8_selectivity_triangulation_report.md`.

## Inputs and parameters

| Script | Produces |
|---|---|
| [`../scripts/collect_results.py`](../scripts/collect_results.py) | `consensus_*.csv` from Vina logs |
| [`../scripts/filter_cascade.py`](../scripts/filter_cascade.py) | `filter_cascade_*.csv` |
| [`../scripts/analyze_selectivity.py`](../scripts/analyze_selectivity.py) | `selectivity_margins.csv` |
| [`../scripts/analyze_water_test.py`](../scripts/analyze_water_test.py) | `water_test_*.csv` |
| [`../scripts/bbb_score.py`](../scripts/bbb_score.py) | BBB / BOILED-Egg classification |

Sign convention in `selectivity_margins.csv`: a **positive** `margin_TTBK1_vs_TTBK2` favours
TTBK1; a **negative** `margin_MAO-B_vs_MAO-A` favours MAO-B.

## Analysis performed

**Filtering cascade.** 405 pool molecules screened; the surviving `candidates_56` set passes all
gates — 56/56 Lipinski, 56/56 GI absorption, 56/56 BBB, 56/56 alert-free. The cascade is a hard
gate, not a ranking: potency is only considered afterwards.

**Consensus scoring.** `consensus_week2.csv` (165 receptor–ligand pairs) is the production table;
`consensus_new.csv` is the earlier pass. Seed SDs are typically ≤ 0.1 kcal/mol.

**Selectivity margins.** The TTBK pair is over all 56 candidates. The MAO pair is over the **48
quantifiable** ones, because 8 candidates have no viable MAO-A pose at all — best mode worse than
−5.0 kcal/mol, a steric failure rather than weak binding — and averaging a non-measurement into a
margin manufactures selectivity out of a docking failure:

| Comparison | n | Mean | Range |
|---|---|---|---|
| TTBK1 vs TTBK2 (validated pair) | 56 | **+0.12** | −0.70 to +1.04 |
| MAO-B vs MAO-A | 48 | **−3.28** (favours MAO-B) | −6.04 to −1.80 |

The 8 excluded candidates — `cand_013`, `020`, `023`, `025`, `029`, `032`, `054`, `056`, the
largest in the set — are *sterically excluded* from MAO-A while still binding MAO-B at −8.7 to
−11.3. That is a qualitative selectivity result, and arguably a stronger one, but it is not a
ΔG and must not be reported as one.

> **Do not average `selectivity_margins.csv` directly.** It carries all 56 rows and does not flag
> the 8 failures, so a plain mean of its MAO column returns **−3.89** — the figure this project
> withdrew on 2026-09-26 for exactly that reason (see `LOGBOOK.md`).
> [`../scripts/analyze_selectivity.py`](../scripts/analyze_selectivity.py) applies the threshold
> and prints the corrected statistics.

**Water test.** Paired wet/dry docking of all 56 candidates and the reference set. The mean
wet-minus-dry shift is only −0.08 kcal/mol, but the range runs to **+3.04** — individual ligands
move by up to 3 kcal/mol depending purely on whether crystallographic waters were retained.

## Final result

**The TTBK1/TTBK2 selectivity claim does not survive scrutiny.** Three independent lines converge:

- The +0.12 kcal/mol mean margin in `selectivity_margins.csv` comes from a receptor pair that
  retained **five waters in one site and none in the other**. That asymmetry, not the chemistry,
  produces the margin.
- Re-run **water-symmetrically**, TTBK2 is favoured by ~1.6 kcal/mol and **0 of 56 candidates**
  are favourable toward TTBK1 — the opposite direction.
- The protocol carries a **~1.0 ± 0.25 kcal/mol systematic bias toward TTBK2** on the 9IV pair, where
  the same-assay experimental ΔΔG is −0.234 kcal/mol, a slight preference the other way
  ([step 5](../05_validation/)). Any margin below **~1.3** — the bias plus its own uncertainty —
  is inside it. The ~1.6 water-symmetric margin clears that floor only barely, and subtracting
  the bias leaves ~0.6.

MD and MM-GBSA were then run to test the claim independently, and found **no resolvable
difference** in either pose stability or binding free energy ([step 6](../06_md/),
[step 7](../07_mmgbsa/)).

**The MAO-B preference is the more robust signal**, at −3.28 kcal/mol mean across the 48
quantifiable candidates, well outside seed noise and favourable in 48 of 48 — with the other 8
sterically excluded from MAO-A altogether. **As of 2026-09-27
MM-GBSA agrees with it.** The earlier sign disagreement (+4.74 favouring MAO-A) was traced to the
one MAO-A trajectory that had drifted off its docked pose; on-pose, **all twelve** MAO-A/MAO-B
pairings favour MAO-B by 2.30–3.94 kcal/mol. The *direction* is corroborated by docking and
MM-GBSA independently, and **since the velocity replicates landed on 2026-09-27 so is the
magnitude**: ΔΔG 2.55 kcal/mol against a measured run-to-run spread of 0.24–0.26, roughly a
tenfold separation. That spread is protocol-bound — FAD is restrained in both arms — so quote
2.55 ± ~0.26 as this protocol's figure. See [step 7](../07_mmgbsa/README.md).

**Lead selected:** `cand_003`, `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F` — #8/56 on 7JXX
(−8.39) and #2/56 on **2V5Z** (−11.41). Ranks are over the 56 candidates, excluding the native
reference ligand in the same file. The lead also held its rank under Vinardo — but that
cross-check was run in the prior phase, on `4NFM`/`2V5Z`
([`prior_phase/phase6_top15_with_vinardo.csv`](prior_phase/phase6_top15_with_vinardo.csv)), **not
on the validated receptors**. Re-running Vinardo on `7JXX` is still open (`WORKFLOW.md`, Phase 3).

*Corrected 2026-09-27:* the −11.41 was previously attributed to 7Q8Y. It is 2V5Z. There is no
`7Q8Y`/`cand_003` row in `consensus_new.csv`. The error credited the lead's best score to the
TTBK2 anti-target instead of MAO-B, the intended target.

## Relevant files

| File | Role |
|---|---|
| [`consensus_week2.csv`](consensus_week2.csv) | **production consensus scores**, all receptor–ligand pairs |
| `consensus_new.csv` | earlier consensus pass, kept for comparison |
| [`selectivity_margins.csv`](selectivity_margins.csv) | per-candidate margins, both target pairs |
| [`filter_cascade_candidates_56.csv`](filter_cascade_candidates_56.csv) | the shortlist with all descriptors and gate outcomes |
| [`filter_cascade_A4_check.csv`](filter_cascade_A4_check.csv) | the full 405-molecule pool screen |
| [`water_test_candidates.csv`](water_test_candidates.csv) | **the water artifact**, per candidate: wet, dry, delta |
| [`water_test_references.csv`](water_test_references.csv) | same for reference inhibitors, with pIC50 |
| `generation_top20_filtered.csv`, `generation_all49_filtered.csv` | RL-generated molecules through the cascade |
| [`prior_phase/`](prior_phase/) | phase 4–8 reports, consensus JSONs, Vinardo cross-check, selectivity triangulation |
