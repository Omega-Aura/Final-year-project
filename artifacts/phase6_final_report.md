# Phase 6 (final) — Consensus Re-Docking + Independent Scoring-Function Validation

## Part A: Multi-seed Vina consensus re-docking (all 56 shortlist candidates)
See `phase6_consensus_report.md` for full method. Summary: AutoDock Vina, exhaustiveness=32, 9 poses,
3 seeds/ligand/target, both TTBK1 (4NFM) and MAO-B (2V5Z). Correlation with the RL-time fast-docking
reward (exhaustiveness=4, same Vina scoring function): TTBK1 r=0.96, MAO-B r=0.74.
This validates only that increasing search effort within the SAME scoring function preserves ranking --
it does not test whether the Vina scoring function itself is a reliable affinity predictor.

## Part B: Independent scoring-function cross-check (Vinardo, top-15 dual-target candidates)
To break this circularity, the top-15 dual-target-ranked candidates were re-docked from scratch (fresh
RDKit conformer embedding per seed, not reused RL-time or Part-A geometries) and scored with the
**Vinardo** empirical scoring function (Quiroga & Villarreal 2016) -- a distinct functional form and
parameterization from the standard Vina scoring function used throughout RL training and Part A.

**Result: cross-scoring-function agreement is weak.**
- TTBK1: Vina-consensus vs. Vinardo, r = 0.56 (moderate)
- MAO-B: Vina-consensus vs. Vinardo, r = -0.28 (no agreement -- essentially uncorrelated,
  slightly negative)

This is a genuine negative finding, not a validation. It indicates that among this top-15 set -- which
is already a narrow, high-scoring slice of chemical space by construction -- the fine-grained rank
ordering by absolute docking score is not robust to the choice of scoring function, particularly for
MAO-B. The two scoring functions likely weight the flavonoid scaffold's H-bond donors/acceptors and the
halogenated/lipophilic B-ring substituents differently, and MAO-B's larger, more hydrophobic and
flexible active-site cavity appears more sensitive to this than TTBK1's kinase ATP pocket.

## Interpretation and what this does NOT invalidate
- Binary hit calls (does this scaffold dock plausibly into both sites at all) are still supported --
  all 15 candidates score favorably by both functions in absolute terms (Vina consensus -7 to -12
  kcal/mol; Vinardo -5 to -9.5 kcal/mol -- Vinardo scores are systematically less negative, a known
  scale offset between the two functions, not a sign of non-binding).
- The Lipinski/PAINS/BOILED-Egg BBB filtering from Phase 5 is unaffected -- those are property-based,
  not docking-score-based.
- What IS undermined: using the fine-grained numeric docking score to pick a single "best" compound
  with confidence. The top-ranked-by-Vina-consensus compound (benzyl+methyl regioisomer,
  `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1Cc1ccccc1`) is NOT the top-ranked-by-Vinardo compound
  (`Cc1ccc(-c2oc3cccc(O)c3c(=O)c2O)cc1Cc1ccccc1`, its positional regioisomer) -- both are the same
  benzyl/methyl-decorated flavonol differing only in substitution pattern on the B-ring, one edging out
  under Vina-only ranking and the other under the combined Vina+Vinardo ranking below.

## Final combined ranking
Rank by the average of z-scored Vina-consensus and z-scored Vinardo scores across both targets
(`final_dual_target_score`, more negative = better under both scoring functions simultaneously) --
this is the most defensible ranking given the cross-function disagreement, since it does not let either
single scoring function's idiosyncrasies dominate. Saved to `phase6_top15_with_vinardo.csv`.

Top-3 by combined score:
1. `Cc1ccc(-c2oc3cccc(O)c3c(=O)c2O)cc1Cc1ccccc1` -- TTBK1 Vina -8.34 / Vinardo -6.91 kcal/mol;
   MAO-B Vina -10.32 / Vinardo -9.51 kcal/mol
2. `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F` -- TTBK1 Vina -7.95 / Vinardo -6.37 kcal/mol;
   MAO-B Vina -11.42 / Vinardo -8.73 kcal/mol
3. `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1Cc1ccccc1` -- TTBK1 Vina -8.15 / Vinardo -6.94 kcal/mol;
   MAO-B Vina -11.49 / Vinardo -6.46 kcal/mol

**Correction (verified by direct recomputation from `phase6_top15_with_vinardo.csv`):** only the combined
#2 candidate (`Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`) is actually robust to scoring-function choice --
it ranks #2 under Vina-consensus-only AND #2 under Vinardo-only. The combined #1 candidate
(`Cc1ccc(-c2oc3cccc(O)c3c(=O)c2O)cc1Cc1ccccc1`) is Vinardo-only #1 but only ~8th of 15 under
Vina-consensus-only -- its weak MAO-B Vina-consensus score (-10.315 kcal/mol, the worst of the top-15 on
that metric) is masked by averaging with its strong Vinardo score, not corroborated by it. The averaged
score is therefore not evidence of cross-function agreement for combined #1; it is evidence of that
agreement only for combined #2, which should be treated as the most defensible single pick from docking
alone. Confirming any ranking further would require orthogonal evidence -- MM-PBSA/MM-GBSA free energies
from MD (blocked, no GPU compute attached) or, ultimately, biochemical assay data.

## Recommendation
Carry forward `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F` (combined rank #2, and the only candidate
ranked #2 or better under BOTH scoring functions independently) as the primary candidate for the
pipeline's next phase (MD/MM-PBSA). Treat the combined-score #1 pick as secondary -- its high average
rank is an artifact of averaging a strong Vinardo score with a weak Vina-consensus MAO-B score, not
evidence of consistent dual-function support. The recurring 7-deoxyflavonol core with 4-methyl/benzyl-type
B-ring decoration remains the structural motif of interest across the shortlist regardless of exact
compound ranking.
