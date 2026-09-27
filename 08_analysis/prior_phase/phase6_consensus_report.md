# Phase 6 — High-Exhaustiveness Consensus Re-Docking

## Method
All 56 Phase-5 shortlisted candidates (pass Lipinski + literal BOILED-Egg BBB gate + PAINS/reactive-alert
free) were re-docked against both on-target receptors at much higher rigor than the RL-time scoring pass:
- AutoDock Vina exhaustiveness = 32 (vs. 4 during RL training)
- 9 output poses per docking run
- 3 independent random seeds (11, 22, 33) per ligand x receptor pair, each with an independent RDKit
  conformer embedding (not just a re-seeded Vina search on the same starting geometry)
- Reported per-ligand: best_score (minimum energy across all seeds/poses), consensus_score (mean of the
  best pose per seed -- a robustness-weighted estimate), and top1_std (std dev across seeds -- pose/seed
  reproducibility)

Receptors and grid boxes (unchanged from Phase 1 preparation / RL training):
- TTBK1: 4NFM (apo; box transferred from 7Q8Y via kinase-domain superposition), center (9.25, 23.18,
  26.94), box 20x20x20 A
- MAO-B: 2V5Z (native FAD cofactor), center (51.89, 156.45, 28.56), box 28x28x28 A (widened from the
  original 20x20x20 A grid-box spec before RL launch)

## Validation against RL-time scoring
The RL reward function used exhaustiveness=4 (single seed) for computational tractability across
~1,000 candidates per run. Comparing the fast RL-time docking scores to the consensus re-dock:
- TTBK1: Pearson r = 0.96 (strong agreement)
- MAO-B: Pearson r = 0.74 (moderate-to-strong agreement, more scatter than TTBK1)

This confirms the RL reward signal was a reasonably faithful proxy for binding affinity ranking,
though MAO-B shows more re-ranking under higher-rigor docking than TTBK1 -- consistent with MAO-B's
larger, more flexible binding cavity.

Seed reproducibility (mean std-dev of best-pose energy across the 3 seeds): TTBK1 0.029
kcal/mol, MAO-B 0.094 kcal/mol -- both indicate stable, reproducible
docking poses (not stochastic search noise driving the rankings).

## Final ranked shortlist
Candidates ranked by combined z-score across both targets' consensus best-pose energies (more negative
combined z-score = stronger dual-target binder relative to the shortlist distribution). Top-15 saved to
`phase6_top15_final_shortlist.csv`; full 56-candidate consensus-scored set in `phase6_consensus_shortlist.csv`.

Top-ranked candidate: `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1Cc1ccccc1` (7-deoxyflavonol core, benzyl + methyl
B-ring decoration) -- TTBK1 consensus best -8.15 kcal/mol, MAO-B consensus best -11.49 kcal/mol, MW
358.4, TPSA 70.7 A^2 (inside BBB-permeable window), WLogP 4.77.

## Next steps
Phase 6b (counter-screening against MAO-A / TTBK2 anti-targets for selectivity) and Phases 7-8
(MD system building, production MD, MM-PBSA binding free energies) remain the next steps. MD phases
are blocked pending an attached GPU compute target.
