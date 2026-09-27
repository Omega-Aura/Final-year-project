# Step 0 — Source Species and Flavonoid Library

## What we did

Selected a plant source and assembled a curated flavonoid library from it, then used that library
as the seed chemical space for reinforcement-learning ligand generation.

[`flavonoid_library.csv`](flavonoid_library.csv) holds 24 compounds from *Evolvulus alsinoides*,
each with: compound name, species, PubChem CID, canonical SMILES, scaffold class, glycoside status
and its aglycone (name, CID, SMILES), a **confidence rating**, a **literature source (PMID)**, and
an RDKit parse check.

## Why we did it

The study is a natural-product-led design campaign, so the chemical space has to come from a real
species with real isolation literature rather than from a generic flavonoid enumeration. Recording
the **aglycone** alongside each glycoside matters because glycosides are large and polar — they
generally fail BBB gating, while their aglycones often pass, and both target families are
intracranial.

Two columns exist specifically to stop unverifiable entries propagating: `confidence` (e.g.
"HIGH — direct isolation paper, species-specific") and `source_reference` (a PMID). A compound
whose provenance is only a review or a database entry is marked as such rather than treated as
equivalent evidence.

`rdkit_parse_ok` is a guard: a SMILES that does not parse cannot be silently carried forward as an
empty molecule.

## Reference

- **Species selection rationale and phytochemistry survey** —
  [`prior_phase/phase0_phytochemistry_report.md`](prior_phase/)
- **Per-compound literature** — the `source_reference` (PMID) column
- **Supporting tables** — `prior_phase/phase0_species_selection.csv`,
  `prior_phase/species_flavonoid_summary.csv`, `prior_phase/compound_classified.csv`
- **Reference inhibitors** — `prior_phase/phase0_reference_inhibitors.csv`, carried forward as
  [`../01_smiles/references.csv`](../01_smiles/references.csv)
- **Generator** — REINVENT4 LibInvent, staged reinforcement learning

## Inputs and parameters

[`reinvent4_output/`](reinvent4_output/) holds the RL campaigns run against this library:

| Item | Role |
|---|---|
| `flavonol_7deoxy_scaffold.smi` | the scaffold the generator decorates |
| `rl_staged_learning_v2.toml` | campaign 2 configuration |
| `campaign1_first_failed/` | first campaign — **did not produce usable candidates**; kept as a negative result |
| `campaign2_v2/` | second campaign — score progression, filtered output, combined run table |

The final campaign (v3, BBB-aware reward) lives in [`../generation/`](../generation/); its
checkpoint is the one retained. The two superseded campaign checkpoints have been moved to
deleted on 27 Sep 2026 (see [`../provenance/`](../provenance/)) — their CSV outputs and
score-progression plots remain here, so the campaign
history is still readable without the 95 MB model weights.

## Analysis performed

Generated molecules were passed through the drug-likeness and BBB cascade in
[step 8](../08_analysis/) and docked in [step 4](../04_docking/). The campaign progression
(`rl_*_score_progression.png`) shows reward improving across staged learning.

## Final result

A 24-compound literature-sourced library with full provenance, and a flavonol scaffold that seeded
three RL campaigns. Campaign 1 failed to yield candidates passing the gates; campaign 2 and the
BBB-aware v3 campaign produced the pool that became the 56-candidate shortlist
([`../01_smiles/candidates_56.csv`](../01_smiles/candidates_56.csv)).

The kept failure matters: campaign 1's configuration used a docking-only reward, and it produced
**zero** candidates passing the literal BBB gate. That is what motivated the BBB-aware reward in
later campaigns — see `../08_analysis/prior_phase/phase5_shortlist_report.md`.

## Relevant files

| Path | Role |
|---|---|
| [`flavonoid_library.csv`](flavonoid_library.csv) | **the curated library** — 24 compounds with provenance |
| `reinvent4_output/flavonol_7deoxy_scaffold.smi` | generator scaffold |
| `reinvent4_output/rl_staged_learning_v2.toml` | campaign 2 config |
| `reinvent4_output/campaign2_v2/` | campaign 2 outputs and score progression |
| `reinvent4_output/campaign1_first_failed/` | campaign 1, retained as a negative result |
| [`prior_phase/`](prior_phase/) | phytochemistry report, species selection, compound classification |
