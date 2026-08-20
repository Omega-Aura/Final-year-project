# Phase 0 — Target Validation Report: TTBK1 + MAO-B Dual-Target Flavonoid CADD Study

## 1. Disease biology summary

**TTBK1** phosphorylates tau at disease-relevant epitopes (including Ser422) and is genetically
and biochemically implicated in Alzheimer's disease and other tauopathies. A growing but smaller
literature also links TTBK1 (and its paralog TTBK2) to TDP-43 phosphorylation and proteinopathy
relevant to ALS/FTLD-TDP. TTBK1 chemical probe development remains at an early, mostly preclinical
stage; the only well-characterized inhibitor series with public measured potency values found in
this search is Halkina et al. 2021 (J Med Chem 64:6358, PMID 33944571) from Biogen, describing
brain-penetrant azaindazole/pyrrolopyridine TTBK1 inhibitors that reduce tau pSer422 in vivo.

**MAO-B** is a validated, clinically established Parkinson's disease drug target. It oxidatively
deaminates dopamine and other monoamines, generating H2O2 as a byproduct; this contributes to
oxidative stress and dopaminergic neurodegeneration. MAO-B is also upregulated in reactive
astrocytes (astrogliosis) in neurodegenerative brain tissue. Three approved/clinical MAO-B
inhibitors (safinamide, selegiline, rasagiline) have extensive measured potency data in ChEMBL.

**Dual TTBK1 + MAO-B targeting — literature assessment (IMPORTANT, HONEST FINDING):**
PubMed searches for a combined "flavonoid + TTBK1" literature and for a generic "dual
TTBK1+MAO-B neurodegeneration" strategy returned **zero hits** in this search (see
`phase0_target_citations.csv` methodology and raw counts in `handoff/pubmed_searches.json`).
By contrast, "flavonoid + MAO-B" searches returned multiple hits — flavonoids as MAO-B
inhibitors is a reasonably well-established sub-field. TTBK1 and MAO-B are not co-implicated
in any single validated disease mechanism in the retrieved literature: TTBK1 is tau/tauopathy
biology (Alzheimer's-spectrum), MAO-B is dopaminergic/oxidative-stress biology (Parkinson's-
spectrum). Both are neurodegeneration-relevant CNS targets, and there is broader precedent for
multi-target-directed ligands (MTDLs) combining kinase inhibition and MAO-B inhibition or
antioxidant activity in neurodegenerative disease drug design generally, but **no literature was
found that proposes or validates TTBK1+MAO-B specifically as a coherent dual-target strategy.**
This project's dual-target framing should be presented explicitly as a **novel, unvalidated
hypothesis** motivated by (a) both targets being licensed drug-discovery targets in
neurodegeneration, (b) flavonoids' documented promiscuity across both kinase and MAO
pharmacology, and (c) the exploratory, hypothesis-generating nature of this study — not as an
extension of an established combination-therapy rationale.

Full citation list (27 references, PMID+DOI+title+year+relevance) saved as
`phase0_target_citations.csv`. Four references named in the project brief were checked: Halkina
2021 J Med Chem confirmed exactly; Bashore 2023 Sci Rep confirmed exactly; Ahamad 2024
Pharmaceuticals confirmed exactly; the brief's "Sharma 2020 Molecules" reference resolved via
PubMed citation-lookup to Chaurasiya et al. 2020, Molecules 25:5358 (same journal/volume/pages/
topic — O-methylated flavonoids vs. MAO-A/MAO-B — but a different first-author name than stated
in the brief). This discrepancy is flagged rather than silently corrected.

## 2. Reference inhibitors (`phase0_reference_inhibitors.csv`)

| Target | Compound | Potency (measured) | Source |
|---|---|---|---|
| MAO-B | Safinamide | IC50 = 7.67 nM | ChEMBL CHEMBL2380255, J Med Chem 2013 |
| MAO-B | Selegiline | IC50 = 7.0 nM | ChEMBL CHEMBL1141467, J Med Chem 2008 |
| MAO-B | Rasagiline | IC50 = 4.4 nM | ChEMBL CHEMBL1141467, J Med Chem 2008 |
| TTBK1 | Halkina 2021 compound 31 | biochemical IC50 = 2.7 nM; cellular IC50 = 315 nM | ChEMBL CHEMBL5280210, doc CHEMBL5236545 = Halkina et al. J Med Chem 2021 (PMID 33944571) |

All four compounds have PubChem CIDs, canonical SMILES, and ChEMBL-sourced measured potency
values with assay descriptions and source documents traceable in the CSV. No well-characterized
TTBK1 inhibitor with public measured potency other than the Halkina series was found — this is
consistent with TTBK1 being an early-stage, non-clinical drug target with essentially one
public high-quality chemical series.

## 3. Structures (`phase0_structures.csv`, files in `structures/`)

### On-targets
- **MAO-B: use 2V5Z, not 2V60.** The project brief specified PDB 2V60 as "MAO-B + safinamide +
  FAD." This is **incorrect** — 2V60 does contain FAD, but its co-crystallized inhibitor is a
  coumarin-4-carbaldehyde analog, not safinamide (both structures are from the same Binda et al.
  2007 J Med Chem paper on human MAO-B complexes with selective noncovalent inhibitors). The
  actual safinamide-bound structure is **2V5Z** (1.60 Å, R-free 0.087; also better resolution
  than 2V60's 2.00 Å / R-free 0.146). Both were downloaded; 2V5Z is recommended for the
  safinamide reference pose and native redocking validation. FAD (HET code FAD) is confirmed
  present, covalently linked, in both. MAO-B active site (safinamide, 4.5 Å) residues: Tyr60,
  Pro102/104, Trp119, Leu164/167/171, Phe168, Cys172, Ile198/199, Gly205, Gln206, Ile316,
  Tyr326, Phe343, Tyr398, Tyr435 — consistent with the published aromatic-cage/substrate-cavity
  architecture. No missing residues near the active site in either structure (only disordered
  N-/C-termini).

- **TTBK1: 4NFM confirmed as the kinase domain, but it is APO (no inhibitor).** 4NFM
  (Homo sapiens, UniProt Q5TCY1, 2.12 Å, R-free 0.182) contains only glycerol (cryoprotectant),
  not a genuine ligand. The catalytic Lys38 (VAIK motif) and DFG motif (Asp176-Phe177-Gly178)
  are both fully resolved, so the ATP pocket is well-defined structurally even without a bound
  ligand. Disordered regions are limited to residues 12–21 (pre-kinase-domain N-terminus) and
  315–343 (C-terminal tail beyond the kinase domain) — neither overlaps the catalytic cleft.
  **Constraint for Phase 1:** because 4NFM has no co-crystal ligand, native-ligand redocking
  validation cannot be performed directly on this structure. Recommended workaround: (a) define
  the grid box centered on the ATP pocket using the resolved Lys38/DFG176 landmarks (pocket
  residues identified by proximity search: Ile48, Gly41, Tyr49, Ile40, Ala61, Lys63, Leu62,
  Glu47, Glu50, Lys39, Gly42), and (b) validate the docking/scoring protocol by redocking known
  ligands into the homologous TTBK2 structures (6U0K, 7Q8Y — see below), which share high
  kinase-domain identity with TTBK1, then transfer the validated protocol parameters to 4NFM.

### Anti-targets
- **MAO-A: 2Z5X confirmed human (UniProt P21397), harmine-bound.** 2.20 Å, R-free 0.191, FAD
  present, fully ordered (no REMARK 465 missing residues). Active site (harmine, 4.5 Å):
  Tyr69, Ile180, Asn181, Ile207, Phe208, Gln215, Cys323, Ile325, Ile335/336/337, Met350,
  Phe352, Tyr407, Tyr444 — the aromatic cage (Tyr407/Tyr444) and gating Phe208 match published
  MAO-A active-site architecture and are appropriately positioned for a selectivity
  counter-screen against the MAO-B aromatic cage (Tyr398/Tyr435).

- **TTBK2: real experimental human structures exist — no AlphaFold fallback needed.**
  Contrary to the brief's assumption, two inhibitor-bound human TTBK2 kinase-domain crystal
  structures were found and verified against UniProt Q6IQ55: **6U0K** (1.74 Å, R-free 0.233,
  bound to inhibitor DTQ) and **7Q8Y** (1.60 Å, R-free 0.080, bound to inhibitor 9IV). Both
  have essentially complete active sites (only minor N-/C-terminal disorder). **7Q8Y is
  recommended as the primary TTBK2 anti-target** given its higher resolution and much lower
  R-free. ATP-pocket residues (9IV contact, 4.5 Å): Ile27, Gly28/29, Ile35, Ala48, Lys50,
  Cys78, Met94, Gln95, Leu96, Gln97, Asn100, Ala102, Asp103, Ser145, Leu162 — directly usable
  as the TTBK2 redocking/grid-box reference, and as the TTBK1-vs-TTBK2 selectivity-relevant
  residue set once mapped onto 4NFM by kinase-domain alignment.

## 4. Constraints and flags carried into later phases
1. Use **2V5Z**, not 2V60, as the MAO-B/safinamide reference structure for Phase 1 receptor
   prep and redocking validation. 2V60 is retained as a downloaded file but is NOT a safinamide
   complex.
2. **4NFM (TTBK1) has no co-crystal ligand.** Native redocking validation must be done via a
   surrogate (TTBK2 6U0K/7Q8Y cross-validation) before trusting docking scores on 4NFM; grid box
   must be defined from the Lys38/DFG176 ATP-pocket landmarks, not from a ligand centroid.
3. **TTBK2 anti-target: use 7Q8Y** (or 6U0K as a second pose) — real experimental structures are
   available; no AlphaFold model is needed for this target, simplifying the Phase 5 selectivity
   counter-screen and avoiding AlphaFold-model-specific caveats (e.g., pLDDT-based pocket
   uncertainty).
4. **Dual TTBK1+MAO-B targeting has no direct literature precedent** — this must be stated
   plainly in the manuscript introduction/limitations as a hypothesis-generating combination,
   not a validated multi-target strategy. Avoid manuscript language implying an established
   dual-target consensus.
5. One project-brief citation (the "Sharma 2020" MAO flavonoid reference) resolves to
   Chaurasiya et al. 2020 Molecules 25:5358 by DOI/journal/topic match, not "Sharma" — verify
   author attribution before citing in the manuscript.
6. TTBK1 (4NFM) and TTBK2 (6U0K/7Q8Y) structures cover different, non-overlapping residue
   ranges/numbering conventions (4NFM: ~res 12-343 of Q5TCY1; TTBK2 structures: ~res 0-299 of
   Q6IQ55) — Phase 5 selectivity analysis will need an explicit sequence alignment between the
   two kinase domains to map equivalent pocket residues, since direct residue-number comparison
   is not valid.
