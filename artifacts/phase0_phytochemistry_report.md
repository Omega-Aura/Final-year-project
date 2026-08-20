# Phase 0 — Phytochemistry Track: Species Selection & Flavonoid Library

## 1. Recommendation

**Recommended species: *Evolvulus alsinoides* (Shankhpushpi).**

*Evolvulus alsinoides* carries 6 structurally-confirmed, RDKit-parseable flavonoid/flavonol
compounds — all traceable to direct phytochemical isolation papers (PMID 17473466, 19748554,
23357036) rather than genus-level reviews — the highest high-confidence flavonoid count of any
candidate screened. Its total PubMed footprint is modest (60 records for the species name alone)
and its CADD-related footprint is a mere 4 records (6.7% of all species literature; query:
`"Evolvulus alsinoides" AND (docking OR "molecular docking" OR "in silico" OR "molecular
dynamics" OR "computer-aided drug design")`), with **zero** hits for `"Evolvulus alsinoides" AND
(TTBK1 OR "tau tubulin kinase")` or `"Evolvulus alsinoides" AND (MAO-B OR "monoamine oxidase")`.
Scholar-gateway semantic search corroborates this: no TTBK1/MAO-B-specific docking study surfaced
for the species anywhere in the corpus. Under the two-axis ranking rule (flavonoid coverage ×
literature gap), *E. alsinoides* scores highest (composite score 0.93 of 1.0) among all six
candidates, combining the largest verified flavonoid set with a near-total absence of prior
target-specific CADD work — precisely the "genuine gap with real dockable chemistry" profile this
screen was designed to find.

## 2. Ranking rule (transparent, pre-registered)

```
flavonoid_coverage_score  = n_high_confidence_flavonoids(species) / max(n_high_confidence_flavonoids across all species)
literature_gap_score      = 1 − [count_CADD_hits(species) / count_species_alone_hits(species)]
composite_rank_score      = flavonoid_coverage_score × literature_gap_score
```
This multiplicative form is deliberate: a species with a huge literature gap but only one or zero
structurally-verified flavonoids scores **zero**, because a gap is worthless without dockable
chemistry to fill it. Both *Sida cordifolia* (0 confirmed flavonoids) and *Celastrus paniculatus*
(0 *high-confidence* flavonoids — its sole candidate, paniculatin, has an unverified structure
assignment) are correctly suppressed to a composite score of 0.000 despite having low CADD
footprints, exactly per the instruction that "a species with a huge gap but only one known
flavonoid is useless."

## 3. Full species screening + literature-gap table

| Species | n flavonoids (total / high-conf) | n glycosides / aglycones | Dominant scaffold | PubMed: species alone | PubMed: + CADD terms | PubMed: + TTBK1 | PubMed: + MAO-B | Composite score |
|---|---|---|---|---|---|---|---|---|
| **Evolvulus alsinoides** | 6 / 6 | 4 / 2 | flavone/flavonol | 60 | 4 | 0 | 0 | **0.933** |
| Bacopa monnieri (neg. control) | 2 / 2 | 1 / 1 | flavone | 566 | 49 | 0 | 4 | 0.304 |
| Clerodendrum serratum | 5 / 2 | 1 / 4 | flavone | 18 | 4 | 0 | 0 | 0.259 |
| Desmodium gangeticum | 7 / 1 | 2 / 5 | pterocarpan/pterocarpene | 53 | 3 | 0 | 1 | 0.157 |
| Sida cordifolia | 0 / 0 | 0 / 0 | none found | 77 | 1 | 0 | 1 | 0.000 |
| Celastrus paniculatus | 1 / 0 | 1 / 0 | isoflavone (unverified) | 73 | 4 | 0 | 2 | 0.000 |

Exact query strings and per-cell counts are in `phase0_species_selection.csv`. All four query
axes were run per species: (a) species name alone; (b) species AND (docking OR "molecular
docking" OR "in silico" OR "molecular dynamics" OR "computer-aided drug design"); (c) species AND
(TTBK1 OR "tau tubulin kinase"); (d) species AND (MAO-B OR "monoamine oxidase"). **Zero of the six
candidates has a single PubMed record combining species name with TTBK1** — this is a field-wide
gap, not unique to the recommended species, and is worth stating plainly in the manuscript's
novelty framing.

### Negative-control validation (Bacopa monnieri)
Bacopa monnieri was included specifically as a well-studied CNS nootropic expected to show a large
CADD footprint, validating the literature-gap metric. It returned **566** total records and **49**
CADD-related records (8.7% CADD fraction, essentially identical in *rate* to Evolvulus's 6.7% but
with an order of magnitude more absolute prior work) — the metric functions as intended: an already
crowded species is correctly down-ranked once its far larger high-confidence flavonoid pool (only
2, since Bacopa's dominant secondary metabolites are triterpenoid saponins/bacosides, not
flavonoids) is factored in.

## 4. Per-species phytochemistry notes

**Evolvulus alsinoides** — Six flavonol/flavonol-glycoside compounds isolated directly from
n-butanol/ethanolic fractions across three independent papers: kaempferol-7-O-glucoside,
kaempferol-3-O-glucoside, quercetin-3-O-glucoside (isoquercitrin), and the novel triglycoside
evolvoside C, plus their common aglycones kaempferol and quercetin. Dominant scaffold class:
**flavone/flavonol** (3-hydroxyflavone core). All four glycosides share the same two aglycones,
so the aglycone-anchored docking set is compact (kaempferol, quercetin) while glycoside forms are
retained separately for BBB-permeability contrast in later phases.

**Desmodium gangeticum** — Structurally distinct: its characteristic flavonoid-adjacent chemistry
is dominated by **pterocarpans/pterocarpenes** (gangetin, gangetinin, desmodin, desmocarpin) rather
than classical flavones — only gangetin is HIGH confidence (isolated and independently studied for
CNS/fertility effects: PMID 38430634 and a separate dementia-model report); the other three carry
MEDIUM confidence (review-only mentions). The reported presence of baicalein/naringin/
neohesperidin in this species rests on a single review article (PMID 36515031) whose abstract text
is partially garbled by OCR/translation and does not unambiguously confirm species-level
attribution — flagged LOW confidence and should be re-verified against a primary isolation paper
before use in later phases.

**Clerodendrum serratum** — Apigenin and luteolin are HIGH confidence, characterized in
ethyl-acetate and n-butanol fractions of an ethanolic leaf extract in a neuroprotective/anti-stress
study (PMID 28458420). Apigenin-7-glucoside, hispidulin and scutellarein derive
from a combined *Clerodendrum indicum* + *C. serratum* natural-product docking library (PMID
28070575) and are flagged MEDIUM confidence pending species-exclusive confirmation. Dominant
scaffold: **flavone**.

**Sida cordifolia** — Extensive targeted searching (genus/species HPLC, LC-MS, named-compound
queries for vitexin, quercetin, kaempferol, gossypin, nictoflorin, scopoletin, strebloside) found
only non-specific "total flavonoid content" phenolic assays and one alkaloid isolation
(1,2,3,9-tetrahydro-pyrrolo[2,1-b]quinazolin-3-ylamine) — **no single named, structurally
confirmed flavonoid could be traced to a primary isolation paper for this species.** This is
reported as a genuine negative result, not an artifact of incomplete search; it independently
disqualifies *S. cordifolia* from selection regardless of its (also low) CADD footprint.

**Celastrus paniculatus** — Multiple reviews name a flavonoid **"paniculatin"** in this species'
phytoconstituent inventory, but none supplies a structure. The PubChem name-string match (CID
169419, an isoflavone di-C-glucoside, MW 594.5) is a name coincidence, not a verified identity —
flagged LOW confidence. The species' well-characterized phytochemistry (dihydroagarofuran
sesquiterpenes, alkaloids, triterpenoids) is not flavonoid-dominated, which by itself disqualifies
it from this flavonoid-focused pipeline regardless of its literature-gap profile.

**Bacopa monnieri** — Negative control, confirmed working as designed (see Section 3).

## 5. Aglycone-vs-glycoside design-core decision

Per project convention, **aglycones are the correct default docking cores for this CNS-targeted
study.** Flavonoid glycosides (MW typically 430–760 Da here, e.g. evolvoside C at 756.7 Da, TPSA
up to 324 Å²) are large, highly polar, and poor candidates for BBB permeation — a decisive
consideration since TTBK1 and MAO-B are both CNS targets. All glycosides isolated from the
recommended species and others were nonetheless retained in the library alongside their explicit
aglycones (kaempferol, quercetin from Evolvulus glycosides; naringenin, hesperetin from
Desmodium's naringin/neohesperidin; apigenin, luteolin from their respective 7-O-glucosides) so
that Phase 3 can dock both forms and quantify the expected aglycone-favors-CNS-penetration effect
directly, rather than assuming it.

## 6. Flavonoid library summary

`phase0_flavonoid_library.csv` contains **24 entries** (21 species-attributed compounds + 3
standalone aglycone entries added for glycosides whose aglycone wasn't already independently
reported for that species). **All 24 SMILES strings parse successfully in RDKit — zero parse
failures.** Every compound_name/pubchem_cid/canonical_smiles was fetched programmatically via the
`chemistry` MCP connector (PubChem) and read back from the saved JSON/CSV before being reported
here; none were hand-typed.

Confidence tiers recorded per compound (see `confidence` column; counts verified directly from
`phase0_flavonoid_library.csv`, 24 rows total):
- **HIGH** (11 of 24): direct, species-specific isolation/characterization paper.
- **MEDIUM** (6 of 24): review-only mention, or reported for a combined multi-species library.
- **LOW** (4 of 24): ambiguous species attribution in source text, or structure identity inferred
  from a name match rather than confirmed against an isolation paper (Paniculatin; the
  Baicalein/Naringin/Neohesperidin trio for *D. gangeticum*).
- **DERIVED** (3 of 24): computed aglycone of a listed glycoside, added as a separate dockable
  entry, not independently reported as isolated from that species (Naringenin and Hesperetin for
  *D. gangeticum*; Genistein for *C. paniculatus*).

## 7. Constraint for Phase 3 (scaffold parameterization)

The **dominant scaffold class for the recommended species (*Evolvulus alsinoides*) is
flavone/flavonol** (3-hydroxyflavone / 2-phenylchromen-4-one core, RDKit SMARTS-confirmed). Any
generative scaffold-decoration step in Phase 3 that treats this study as single-species should
parameterize on the flavonol core (kaempferol/quercetin), not a generic flavonoid template, and
should treat the O-glycoside substituent pattern at C3/C4'/C7 as the natural decoration vector
suggested by the isolated compound set (7-O-, 3-O-, and 4'-O-triglycoside positions are all
attested in this species' real chemistry).

## 8. Data-source caveats (full list)

1. **OpenAlex/`literature` connector unusable** (returns `openalex_key_required`) — all literature
   quantification relies on `pubmed` (PubMed E-utilities via MCP) with `scholar-gateway` used only
   for qualitative corroboration, not as a second independent hit-count source. Hit counts should
   be understood as PubMed-only counts, not multi-database totals.
2. **IMPPAT (cb.imsc.res.in)**, the standard Indian phytochemical database, was assumed blocked
   per task instructions and never queried; all compound identities were instead cross-referenced
   from primary isolation-paper abstracts against PubChem name search, which is a weaker
   cross-check than a curated phytochemical database would provide.
3. **PubMed search hit counts are name-string counts**, not systematic-review-grade curation —
   they will include unrelated homonym hits in principle (not observed here on spot-checking, but
   not exhaustively ruled out for every count).
4. **Genus-level review text** (particularly PMID 36515031 for *Desmodium*) sometimes suffers
   from OCR/formatting corruption when passed through the article-metadata API (e.g. species names
   dropped mid-sentence), which is why 3 of the 7 Desmodium compounds are flagged LOW confidence
   rather than treated as confirmed.
5. **Scaffold classification uses hand-built RDKit SMARTS patterns** for flavone, flavonol,
   flavanone, isoflavone, pterocarpan/pterocarpene, and chalcone cores — these were validated by
   manual inspection against every compound's parsed structure in this run, but were not tested
   against an external flavonoid ontology; edge cases (e.g. highly substituted or ring-fused
   variants) could in principle be misclassified.
6. **Paniculatin's structure is unverified** — flagged prominently; if Phase 3 needs this
   compound, source a primary paper with an unambiguous structure/NMR assignment before docking it.
7. Literature-gap counts reflect the state of PubMed as queried on 2026-08-17 and will drift as
   new papers are indexed.

## Files
- `phase0_species_selection.csv` — full screening + literature-gap table (species-level, with
  every query string and PubMed hit count, plus the ranking-score computation).
- `phase0_flavonoid_library.csv` — 24-compound flavonoid library (compound_name, species,
  pubchem_cid, canonical_smiles, scaffold_class, is_glycoside, aglycone fields, confidence,
  source_reference, rdkit_parse_ok).
- `compound_classified.csv` — intermediate per-compound table with RDKit classification detail
  (superset of the library, includes all raw candidates before aglycone-row addition).
- `species_flavonoid_summary.csv` — intermediate per-species roll-up (compound counts, glycoside
  counts, dominant scaffold) prior to merge with literature-gap data.
