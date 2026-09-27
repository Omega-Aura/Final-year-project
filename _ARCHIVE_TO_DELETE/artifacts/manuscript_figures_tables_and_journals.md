# Manuscript Figure & Table Index

Companion to `manuscript_dualtarget_flavonol.md`. Every entry is a saved project artifact.

## Main figures

| # | Artifact | Section | Content |
|---|---|---|---|
| 1 | `rl_v2_score_progression.png` | 3.2 | Composite reward trajectory across the 32-step BBB-aware RL run |
| 2 | `phase6_vinardo_crosscheck_scatter.png` | 3.3 | Vina-consensus vs. Vinardo scatter, top-15, both targets |
| 3 | `pose_primary_TTBK1.png` | 3.4 | Lead pose in the TTBK1 ATP pocket (4NFM), PyMOL render |
| 4 | `pose_primary_MAOB.png` | 3.4 | Lead pose in the MAO-B substrate cavity (2V5Z), FAD-adjacent |
| 5 | `md_trajectory_analysis.png` | 3.5 | 4-panel: protein/ligand RMSD and per-residue RMSF, both systems, 20 ns |
| 6 | `mmgbsa_comparison.png` | 3.6 | Total MM-GBSA dG per target + per-component decomposition |
| 7 | `selectivity_triangulation.png` | 3.7 | On- vs. anti-target docking margins + cross-method on-target confirmation |

## Supplementary figures

| # | Artifact | Content |
|---|---|---|
| S1 | `pose_benzylMeA_TTBK1.png` / `pose_benzylMeA_MAOB.png` | Vinardo-preferred regioisomer A poses, both targets |
| S2 | `pose_benzylMeB_TTBK1.png` / `pose_benzylMeB_MAOB.png` | Regioisomer B poses, both targets |

## Main tables

| # | Source artifact | Section | Content |
|---|---|---|---|
| 1 | `phase0_species_selection.csv` | 2.1 | Six-species ranking: flavonoid coverage x literature gap |
| 2 | `phase1_summary.json` | 2.2 | Receptor preparation status and redocking validation RMSDs |
| 3 | (in text) | 3.1-3.2 | Filtering cascade attrition, both generative campaigns |
| 4 | `phase6_top15_with_vinardo.csv` | 3.3 | Top-15 consensus + Vinardo scores, both targets |
| 5 | `phase5_top3_interaction_fingerprints.csv` | 3.4 | ProLIF contact residues, 3 leads x 2 targets |
| 6 | `mmgbsa_TTBK1_4NFM_results.dat` / `mmgbsa_MAOB_2V5Z_results.dat` | 3.6 | MM-GBSA component energies |
| 7 | `phase6b_selectivity_top15.csv` | 3.7 | Anti-target margins, all 15 candidates x 3 anti-target receptors |

## Supplementary data

| Artifact | Content |
|---|---|
| `phase0_flavonoid_library.csv` | 24-entry natural-product library, SMILES + PMID provenance + confidence tier |
| `phase0_reference_inhibitors.csv` | Reference inhibitor measured potencies (ChEMBL) |
| `phase0_target_citations.csv` | PubMed query strings and raw hit counts underlying the novelty assessment |
| `phase0_structures.csv` | PDB structure selection with resolution and ligand annotations |
| `rl_v2_shortlist_56.csv` | Full 56-candidate all-filters-pass shortlist with descriptors |
| `phase6_consensus_shortlist.csv` | Consensus docking, all 56 candidates, both on-targets |
| `mmgbsa_TTBK1_4NFM_results.csv` / `mmgbsa_MAOB_2V5Z_results.csv` | Per-frame MM-GBSA energies (n = 100 each) |
| `selectivity_triangulation.json` | Machine-readable triangulation record for the lead |

---

# Target Journal Recommendations

The manuscript's defining characteristic for journal fit is that **its two most solid
findings are negative** (the scaffold BBB ceiling; scoring-function non-robustness) and its
positive finding is explicitly hypothesis-generating with an unresolved selectivity
liability and no experimental validation. Venues that require wet-lab confirmation of a
computational lead are not viable, and venues that reward "we found a hit" framing would
require overselling. The recommendations below are ordered by fit.

## 1. *Journal of Cheminformatics* (Springer Nature) — recommended primary

**Fit rationale.** Open access, cheminformatics-methods scope, and an established
willingness to publish rigorous negative and cautionary methodological results. The
Vinardo cross-check finding — that fine-grained docking rank order within a narrow
high-scoring slice is not robust to scoring-function choice, negatively correlated at MAO-B
— is a genuine contribution to virtual-screening practice and is squarely in scope. The
BBB TPSA-floor result is a reusable constraint for anyone designing CNS-directed
flavonoids. REINVENT4 is published in this journal, which helps the generative-design
framing land.

**Positioning.** Lead with the methodology: a BBB-constrained generative pipeline with
built-in anti-target counter-screening and cross-scoring-function validation, demonstrated
on a novel dual-target hypothesis. The compound is the case study, not the headline.

**Required before submission.** No wet-lab work needed. Complete the reference-14
bibliographic details; consider adding the anti-target MM-GBSA calculations (Section 4.5),
which would meaningfully strengthen the selectivity argument at modest compute cost.

## 2. *Molecules* (MDPI), section "Computational and Theoretical Chemistry" or "Medicinal Chemistry"

**Fit rationale.** Publishes natural-product CADD studies at this level of computational
depth without requiring experimental validation, and both cited flavonoid/MAO-B precedents
(ref. 13) appear there, so the audience is directly addressable. Rapid review. The
*Evolvulus alsinoides* phytochemistry angle fits the journal's natural-products readership.

**Caveat.** This is the venue where the risk of the paper being read as a routine
"docking + MD of a natural product hit" study is highest, since many such papers appear
there. The negative findings and the selectivity liability must be foregrounded in the
abstract and title, not deferred to the discussion, or the contribution will be
misread as one more virtual-screening hit report.

## 3. *Frontiers in Chemistry* / *Frontiers in Molecular Biosciences*, computational section

**Fit rationale.** Explicitly accepts hypothesis-generating computational work and
negative results; open review adds transparency that suits a paper whose main claim is
methodological caution. Good venue for the "here is a target-pair hypothesis and here is
exactly how far the computation supports it" framing.

## 4. *International Journal of Molecular Sciences* (MDPI) — fallback

**Fit rationale.** Broad scope, accepts purely computational MTDL studies. Lower
specificity of audience than the above; use if the methodological framing does not land at
option 1.

## Not recommended

- **Journal of Medicinal Chemistry / European Journal of Medicinal Chemistry.** Both
  effectively require synthesis and biochemical assay data for a new chemical series. A
  purely computational lead with unresolved isoform selectivity will not clear review.
- **ACS Chemical Neuroscience.** Would require cellular or in vivo CNS data to support the
  neurodegeneration framing.
- **Any journal where the dual-target rationale must be presented as established.** The
  PubMed result is zero records for a TTBK1 + MAO-B strategy; a venue that expects a
  validated target-pair rationale in the introduction is a poor fit and would pressure the
  manuscript toward overclaiming.

## Title options

1. "A blood-brain-barrier ceiling on the flavonol scaffold, and what survives a
   scoring-function change: generative dual-target design against TTBK1 and MAO-B"
   *(methodology-forward; best fit for option 1)*
2. "Dual engagement without isoform selectivity: computational design and counter-screening
   of 7-deoxyflavonol candidates against TTBK1 and MAO-B"
   *(finding-forward, honest about the liability; best fit for options 2-3)*
3. "Generative design of BBB-permeant 7-deoxyflavonols as dual TTBK1/MAO-B candidates:
   consensus docking, molecular dynamics, MM-GBSA and anti-target counter-screening"
   *(descriptive/conventional; safest for option 4)*
