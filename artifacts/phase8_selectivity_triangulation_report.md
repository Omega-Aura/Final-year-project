# Phase 8 — Dual-Target Selectivity Triangulation

## Candidate
Primary lead: the 4-CF3/methyl-substituted flavonol carried through molecular dynamics
(canonical SMILES `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`), the single candidate advanced to MD given the
6 GB single-GPU compute budget (see Phase 6 consensus docking / Phase 6b anti-target
counter-screen for the full 15-candidate shortlist context).

## Method
Three independent lines of evidence are triangulated for this candidate:
1. **On-target consensus docking** (Vina, exhaustiveness=32, 9 poses, 3 seeds; Phase 6)
2. **Anti-target consensus docking** against TTBK2 (two receptor options, 6U0K/7Q8Y) and
   MAO-A (2Z5X) (Phase 6b)
3. **20 ns explicit-solvent MD + single-trajectory MM-GBSA** binding free energy on both
   on-targets (this phase)

## On-target engagement: confirmed by two independent methods

| Target | Docking (Vina consensus best, kcal/mol) | MM-GBSA ΔG (kcal/mol, mean ± SEM, n=100 frames) |
|---|---|---|
| TTBK1 (4NFM) | -7.95 | -26.32 ± 0.53 |
| MAO-B (2V5Z) | -11.42 | -36.85 ± 0.23 |

Both docking and post-MD MM-GBSA agree on the qualitative ranking: **MAO-B engagement is
predicted stronger than TTBK1 engagement** for this candidate, by both methods (docking
delta -3.47 kcal/mol; MM-GBSA delta
-10.53 kcal/mol). The MM-GBSA magnitudes are not directly
comparable to docking scores (different scoring functions, absolute values are not on the
same scale; MM-GBSA also implicitly favors the larger/more polar cofactor-adjacent MAO-B
pocket), but the *direction* of the on-target preference is consistent, which is the
triangulation this step is meant to establish rather than asserting dual-target activity
from docking alone.

The MM-GBSA per-component decomposition (`mmgbsa_comparison.png`) shows MAO-B's more
favorable ΔG is driven primarily by a larger van der Waals contribution
(-44.9 vs.
-32.9 kcal/mol for TTBK1), consistent with
the larger, more enclosing MAO-B/FAD-adjacent hydrophobic pocket versus the shallower
kinase hinge region in TTBK1.

## Anti-target selectivity (docking only — MD not run on anti-targets given compute budget)

| Comparison | Selectivity margin (kcal/mol, on-target − anti-target consensus best) | Interpretation |
|---|---|---|
| TTBK1 vs. TTBK2 (6U0K) | +1.71 | positive = TTBK2 binds more strongly (unfavorable) |
| TTBK1 vs. TTBK2 (7Q8Y) | +2.23 | positive = TTBK2 binds more strongly (unfavorable) |
| MAO-B vs. MAO-A | +0.09 | near-zero = no separation |

Per the Phase 6b counter-screen, **none of the three anti-target receptor structures
passed redocking validation**, so these margins carry additional uncertainty beyond the
usual docking error and should be read as directional signal, not validated affinity
predictions.

## Triangulated conclusion

1. **On-target engagement of both TTBK1 and MAO-B is corroborated by an independent,
   higher-level-of-theory method** (20 ns MD + MM-GBSA), strengthening confidence in the
   docking-stage dual-target hypothesis beyond what docking alone could support. MAO-B
   engagement is the more robust of the two by both methods.
2. **Isoform/paralog selectivity remains unresolved and is the dominant risk for this
   series.** The docking-only anti-target data indicate this flavonol scaffold does not
   discriminate TTBK1 from TTBK2 (0/15 shortlist candidates favorable in Phase 6b) and
   only weakly/inconsistently discriminates MAO-B from MAO-A. MD/MM-GBSA was not extended
   to the anti-targets in this phase (single-GPU budget was allocated to validating
   on-target engagement first); resolving selectivity definitively would require the same
   MM-GBSA treatment applied to TTBK2 and MAO-A complexes, or biochemical isoform-selectivity
   assays.
3. **Net assessment:** this candidate should be reported as a validated **dual on-target
   engager** (TTBK1 and MAO-B, cross-validated by docking and MD/MM-GBSA) with an
   **open, unresolved isoform-selectivity liability** rather than as a selective dual-target
   hit. This distinction should be stated explicitly in the manuscript to avoid
   overselling the selectivity claim.

Data: `selectivity_triangulation.json`. Figures: `mmgbsa_comparison.png`,
`selectivity_triangulation.png`.
