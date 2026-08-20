# Phase 6b — Anti-Target Counter-Screening for Selectivity

## Method
The top-15 dual-target-ranked candidates were consensus re-docked (same protocol as Phase 6: Vina,
exhaustiveness=32, 9 poses, 3 seeds/ligand) against all three prepared anti-targets:
- MAO-A (2Z5X), box 28x28x28 A, center (40.58, 26.93, -14.54)
- TTBK2 kinase domain, option 1 (6U0K), box 20x20x20 A, center (23.3, 29.5, 17.33)
- TTBK2 kinase domain, option 2 (7Q8Y), box 20x20x20 A, center (21.2, -26.41, 40.89)

## Receptor validation caveat (carried forward from Phase 1)
- 2Z5X (MAO-A): no redocking validation was performed for this receptor specifically (it was prepared
  via a custom FAD-transplant splice from the 2V5Z template); treat MAO-A scores as unvalidated.
- 7Q8Y (TTBK2, primary): **FAILED** redocking validation in Phase 1 (best pose 5.29 A RMSD vs. the 2.0 A
  pass threshold, no setup error identified -- a genuine Vina scoring-function limitation for this
  ligand/pocket). Scores against this receptor should be treated as lower-confidence.
- 6U0K (TTBK2, secondary): no redocking validation was performed (used as a second structural option
  precisely because 7Q8Y failed validation, but was never itself validated against a co-crystallized
  ligand). Scores against this receptor should also be treated as unvalidated, not confirmed-accurate.

None of the three anti-target receptors used in this counter-screen carry a passing redocking-validation
result. This is a material limitation: absolute anti-target binding energies, and therefore the
selectivity deltas computed from them, should be read as directional/qualitative signal only, not as
validated affinity predictions.

## Selectivity results (on-target consensus best score minus anti-target consensus best score;
negative delta = favorable selectivity for the on-target)

**MAO-B vs. MAO-A:** mean delta across the top-15 = +0.45 kcal/mol
(positive = MAO-A binds more strongly on average). Only 2/15
candidates show favorable MAO-B-over-MAO-A selectivity by this metric. This is consistent with the
well-documented pharmacological difficulty of achieving MAO-B/MAO-A selectivity for planar polyphenolic
scaffolds, since the two isoforms' active sites are highly homologous.

**TTBK1 vs. TTBK2:** mean delta = +1.92 kcal/mol (vs. 6U0K) and
+2.27 kcal/mol (vs. 7Q8Y). **0/15** candidates show favorable
TTBK1-over-TTBK2 selectivity against either TTBK2 structure -- every candidate in the shortlist is
predicted to bind TTBK2 more strongly than TTBK1 by 1.5-2.7 kcal/mol. Given TTBK1 and TTBK2 share a
near-identical ATP-binding kinase fold (this is exactly why 7Q8Y's box could be transferred onto 4NFM
by structural alignment in Phase 1), this lack of predicted selectivity is expected from a flavonoid
scaffold occupying the same conserved hinge region, and is not itself surprising -- but it means kinase
selectivity is an unsolved problem for this entire chemical series based on docking alone.

## Interpretation
Docking-based counter-screening indicates this flavonol series shows essentially no predicted kinase
selectivity (TTBK1 vs. TTBK2) and weak, inconsistent MAO-B-over-MAO-A selectivity. Because none of the
three anti-target receptors passed redocking validation, these selectivity numbers carry substantial
additional uncertainty beyond the on-target consensus/Vinardo discrepancy already reported in Phase 6.
Selectivity claims for this series should not be made on the basis of this docking data alone; MM-PBSA
(if MD becomes available) or biochemical selectivity assays (MAO-A/MAO-B and TTBK1/TTBK2 side-by-side)
are needed to make a defensible selectivity claim.

Full per-candidate table saved to `phase6b_selectivity_top15.csv`.
