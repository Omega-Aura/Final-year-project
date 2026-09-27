# Phase 5 — Novelty Check & Shortlist Selection (Corrected v2, BBB-Aware RL)

## Context
Phase 4 (flavonol scaffold, docking-only reward) produced **zero** candidates passing the literal
BOILED-Egg BBB gate — the flavonol's mandatory 3-OH/5-OH/7-OH substitution pattern fixes TPSA at a
minimum of 90.9 Å², outside the BBB-permeable region regardless of R-group decoration.

In response, the scaffold was redesigned to a mono-deoxy flavonol (**7-deoxyflavonol**, retaining
the 5-OH/4-carbonyl intramolecular H-bond and removing the 7-OH) and the RL reward function was
rebuilt with a fourth, direct BBB signed-distance component (TTBK1 0.35 / MAO-B 0.35 / SA 0.15 /
BBB 0.15, geometric mean). A first scaffold-construction attempt inadvertently built the wrong
regiochemistry (retained 7-OH, dropped 5-OH) — this was caught, corrected, and the full 32-step RL
run (4 chunks x 8 steps) was executed from scratch on the corrected scaffold.

## RL run summary (corrected v2)
- 1,024 total generated candidates across 32 steps, 412 unique SMILES
- Mean composite reward score rose from 0.08 (step 1) to a plateau of 0.45-0.53 (steps ~20-32)
- Overall BBB-pass rate across the run: 389/1024 (38%). Per-step BBB-pass fraction is noisy
  throughout, not a tight plateau. Steps 1-2 (still close to the unmodified prior) start high
  (71.9%, 62.5%), then drop sharply as RL begins optimizing docking affinity (steps 3-8 range
  0-34.4%, mean 16.9%) — this is the expected affinity/permeability trade-off the BBB reward term
  is meant to correct. From step ~9 onward the BBB-aware reward pulls the pass rate back up; across
  steps 13-32 it averages 45% but still ranges from 12.5% to 65.6% step-to-step (32-candidate
  batches per step), with no further downward trend

## Phase 5 filtering (novelty + drug-likeness + BBB + PAINS)
Run on the full corrected v2 combined library (1,024 rows) using the same pipeline validated in Phase 4:

| Filter stage | Count |
|---|---|
| Input rows | 1,024 |
| Invalid SMILES | 7 |
| Duplicate of seed (natural product) library | 0 |
| Internal duplicates | 612 |
| Unique, valid, novel | 405 |
| Lipinski pass | 369 |
| **BBB pass (literal BOILED-Egg gate)** | **81** |
| GI absorption pass | 334 |
| PAINS/reactive-alert free | 253 |
| **Pass ALL filters (Lipinski + BBB + alert-free)** | **56** |

This is the qualitative reversal the scaffold redesign targeted: **56 non-empty shortlisted
candidates**, versus zero for the original flavonol scaffold.

## Shortlist characteristics (56 candidates passing all filters)
- MW range: 310-420 Da (all within Lipinski bounds)
- TPSA: 70.7-77.2 A^2 (inside the BOILED-Egg BBB ellipse)
- cLogP (Wildman-Crippen): 3.6-4.8
- TTBK1 docking (raw): mean 7.49, range 6.90-8.32
- MAO-B docking (raw): mean 10.52, range 6.83-11.76
- 8 unique Murcko scaffolds (R-group decoration on the shared 7-deoxyflavonol core drives most
  of the structural diversity, consistent with the LibInvent decoration strategy)

## Top-20 shortlist
Ranked by composite RL score (balances TTBK1 docking, MAO-B docking, synthetic accessibility, and
BBB permeability). Leading candidates carry benzyl, methyl/trifluoromethyl, or cyclic-amine
substituents on the pendant B-ring, retaining strong dual-target docking scores while sitting
inside the BBB-permeable TPSA/LogP window.

Full 56-candidate shortlist: `rl_v2_shortlist_56.csv`
Top-20 subset: `rl_v2_top20_shortlist.csv`

## Next steps
Proceed to Phase 6 (high-exhaustiveness consensus re-docking) and Phase 6b (counter-screen against
anti-targets for selectivity) on this shortlist. Phases 7-8 (MD system building, production MD,
MM-PBSA) remain blocked pending an attached GPU compute target.
