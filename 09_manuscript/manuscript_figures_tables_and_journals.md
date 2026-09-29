# Figures, tables, supplementary data and target journals

Companion to [`manuscript_dualtarget_flavonol.md`](manuscript_dualtarget_flavonol.md).
Regenerated 2026-09-27 against the Week 1–2 dataset; the August version is in
[`prior_phase/`](prior_phase/) and refers to figures and receptors that are no longer used.

---

## Figures

Three figures exist as files and are drawn from the collected summary table by
[`../10_results/make_figures.py`](../10_results/make_figures.py), so a figure cannot disagree
with the text. Re-run that script after any new MD; it regenerates all three.

| # | File | Section | Content |
|---|---|---|---|
| 1 | [`../10_results/fig1_pose_stability.png`](../10_results/fig1_pose_stability.png) | 3.7 | Last-100-frame ligand RMSD for all seventeen systems, banded by verdict. **The paper's lead figure** — it carries the docking-rank-versus-stability result. |
| 2 | [`../10_results/fig2_mao_binding_energy.png`](../10_results/fig2_mao_binding_energy.png) | 3.8 | MAO-A against MAO-B on-pose ΔG with replicate spread. **MAO-only by design** — see the axis rule below. |
| 3 | [`../10_results/fig3_sem_vs_replicate.png`](../10_results/fig3_sem_vs_replicate.png) | 3.9 | Reported `MMPBSA.py` SEM against the spread measured from a second velocity seed, for all five replicate pairs. |
| 4 | [`../05_validation/benchmark_mao_scatter.png`](../05_validation/benchmark_mao_scatter.png) | 2.2, 3.4 | MAO reference-ligand docking score against measured potency. **Label as mixed-species** — two of its points are rat-brain assays. |

**Two figures still to draw**, both optional and neither load-bearing:

- **Filtering cascade / BBB ceiling** (§3.1–3.2): the TPSA-versus-WLogP plane with the BBB
  region marked, showing the native-flavonol candidates sitting wholly outside it and the
  7-deoxy set inside. This is the most publication-friendly way to show §3.1, which is currently
  text-only. Data: `../08_analysis/filter_cascade_candidates_56.csv` and the 405-molecule pool
  screen.
- **Lead poses in both sites** (§3.3): the docked pose in the TTBK1 ATP pocket and the MAO-B
  substrate cavity. The August draft had these rendered from `4NFM`; they must be **re-rendered
  on `7JXX`** if used at all, and the caption must name the pose number, since §3.7 establishes
  that pose identity matters.

**Do not reuse any figure from the August draft.** All of them were rendered on the earlier
receptor set (`4NFM`, `6U0K`, `2V60`) or from the 20 ns prior-phase MD, and several encode
numbers that have since been withdrawn.

### The one hard axis rule

**Never put TTBK and MAO binding energies on the same axis.** They carry protein-specific
desolvation and surface terms that do not cancel; MAO-A's −44.22 against TTBK1's −31.36 compares
a flavoenzyme to a kinase and means nothing. Figure 2 is MAO-only for exactly this reason.
Figure 1 shares an axis legitimately, because ligand RMSD from the docked pose is the same
quantity in every system.

**And never plot or quote the −44.22** as a MAO-A binding energy. It is the most attractive
number in the dataset, it is the only off-pose MAO-A value, and it inverted the MAO conclusion
for two days. It appears in the manuscript once, in §4.3, as the cautionary case.

---

## Tables

| # | Section | Content | Source |
|---|---|---|---|
| 1 | 2.2 | The six receptors, their proteins and native ligands | `../03_receptors/README.md` |
| 2 | 3.3 | Lead compound properties (MW, WLogP, TPSA, HBD/HBA, Lipinski, BBB, GI, alerts) | `../08_analysis/filter_cascade_candidates_56.csv` |
| 3 | 3.4 | Redocking validation, all six receptors, best RMSD and verdict | `../05_validation/redock/` |
| 4 | 3.5 | The 9IV calibration: docking margin, experimental ΔΔG, derived bias, interpretability floor | `../05_validation/calibration_9IV_margin.json` |
| 5 | 3.7 | All seventeen MD systems: target, structure, starting pose, last-100 RMSD, verdict | `../07_mmgbsa/md_mmgbsa_summary.csv` |
| 6 | 3.8 | On-pose MM-GBSA ΔG per system with SD, flagged for selectivity usability | same |
| 7 | 3.9 | Five replicate pairs: mean ligand RMSD, measured spread, ratio to reported SEM | same |

Tables 5–7 all read from one generated CSV. That is deliberate: it is the single file the whole
Results section is written from, and `scripts/collect_md_summary.py --check` fails if it is stale
relative to the primary `.dat` files.

---

## Supplementary data

| Item | Location |
|---|---|
| Flavonoid library with provenance and confidence tiers | `../00_library/flavonoid_library.csv` |
| Reference set with traced primary sources, assay format and species | `../01_smiles/references.csv` |
| Full generative-run output and reward trace | `../00_library/reinvent4_output/`, `../generation/` |
| Filtering cascade, 405-molecule pool and the final 56 | `../08_analysis/filter_cascade_*.csv` |
| Docking consensus, every receptor–ligand pair | `../08_analysis/consensus_week2.csv` |
| Selectivity margins per candidate, both pairs | `../08_analysis/selectivity_margins.csv` |
| Water-shell test, candidates and references | `../08_analysis/water_test_*.csv` |
| Prepared receptors (incl. dry / noFAD / bridging-water controls) | `../03_receptors/` |
| Prepared ligands | `../02_ligands/` |
| All docking runs | `../04_docking/` |
| Per-system MD trajectories, RMSD/RMSF, raw `MMPBSA.py` output | `../06_md/systems/` |
| MD + MM-GBSA summary table | `../07_mmgbsa/md_mmgbsa_summary.csv` |
| Chronological record of every decision and defect | `../LOGBOOK.md` |
| Protocol as specified, with divergences flagged | `../WORKFLOW.md` |

**On quoting `selectivity_margins.csv` directly:** it carries all 56 rows and does not flag the
eight candidates with no viable MAO-A pose, so a plain mean of its MAO column returns −3.89
rather than the correct −3.28 over 48 quantifiable candidates. Use
`scripts/analyze_selectivity.py`, which applies the exclusion.

---

## Target journals

The contribution has shifted since the August assessment. It is now **as much a methods paper as
a compound paper**: the two strongest results are that docking rank does not predict pose
stability (3 of 4 targets) and that a single-trajectory SEM understates the true spread by 2–33×
in proportion to pose instability. That changes which venues fit.

| Venue | Fit | Notes |
|---|---|---|
| ***J Cheminformatics*** | **Recommended primary** | Publishes negative and methodological results; the two methods findings are the paper's strongest claims and this is the natural home for them. Open access, no wet-lab expectation. |
| ***J Chem Inf Model*** | **Strong alternative** | Better reach for the pose-stability/SEM findings specifically. Expects methodological rigour, which the calibration and replicate design supply. |
| ***Molecules*** / ***IJMS*** | Viable fallback | Fast and receptive to computational MTDL studies, but tends to reward positive framing — there is a real risk the negative findings get read as weaknesses rather than results. |
| ***Frontiers in Chemistry*** / ***Mol Biosciences*** | Viable | Similar profile to the above. |
| **Any medicinal-chemistry journal expecting synthesis and assay** | **Do not submit** | A purely computational lead with one target pair's selectivity unresolved will not clear review. |
| **Any venue where the dual-target rationale must be presented as established** | **Do not submit** | It is a hypothesis with zero supporting literature, and the manuscript says so. |

### Title options

1. *(current)* **A generative 7-deoxyflavonol against TTBK1 and MAO-B: isoform selectivity
   resolvable for one target pair and not the other, and why docking rank does not predict pose
   stability**
2. **Docking rank does not predict pose stability: lessons from a dual-target flavonol against
   TTBK1 and MAO-B** — leads with the methods finding; best for *JCIM*.
3. **What a docking-plus-MM-GBSA pipeline can and cannot resolve about isoform selectivity:
   a calibrated case study on TTBK1/TTBK2 and MAO-A/MAO-B** — most honest framing of the
   contribution; weakest as a discovery narrative.

### Framing requirements, whichever venue

These follow from §4 and are not stylistic preferences:

1. **The abstract and title must foreground the negative and methodological findings.** Deferring
   them to the Discussion misrepresents the contribution — and, on this dataset, the negatives
   are better supported than the positives.
2. **The two target pairs must never be summarised as one selectivity result.** They reach
   opposite conclusions.
3. **Any TTBK2 liability statement must be attributed to docking scores alone**, with the
   explicit note that dynamics lacks the resolution to confirm or refute it.
4. **Report the replicate spread, never the `MMPBSA.py` SEM**, as the uncertainty on any ΔΔG.
5. **Quote the calibration bias as ~1.0, never 1.018.** Its inputs are two unreplicated IC50s
   and do not support three significant figures.
6. **State that the score-versus-potency correlation is mixed-species**, and that its two
   rat-brain points are the two flavonoids.
7. **Resolved limitations were deleted, not softened.** The August draft's "no replicate MD" and
   "no anti-target validation" limitations are gone because the work was done, not because the
   wording was weakened. Do not reinstate them.
