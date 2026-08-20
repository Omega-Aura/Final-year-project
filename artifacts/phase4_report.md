# Phase 4 Filtering Report — RL-Optimized Flavonol Library

## Pipeline
Input: 1,280 scored molecules from the 40-step staged-learning RL run (REINVENT4 LibInvent,
dual-target TTBK1+MAO-B docking + SA-score reward), across 5 checkpointed chunks.

Filters applied in order (see `filtering/run_filters.py`):
1. RDKit sanitization + canonical dedup, novelty check against the 23-compound
   seed flavonoid library (phase0_flavonoid_library.csv, incl. aglycones)
2. Lipinski Ro5 + extended druglikeness descriptors
3. BOILED-Egg BBB/GI permeability gate (Daina & Zoete 2016 ellipses, digitized from pyBOILEDegg,
   computed from RDKit WLOGP/TPSA) — MANDATORY per project brief since both targets (TTBK1, MAO-B)
   are CNS-relevant (TTBK1 in neurons, MAO-B predominantly in glia)
4. PAINS (A/B/C) + BRENK reactive-group alerts (RDKit FilterCatalog)

## Attrition
| Stage                          | Count |
|---------------------------------|------:|
| Input (RL-scored rows)          | 1,280 |
| Invalid SMILES                  |    12 |
| Duplicate of seed library        |     0 |
| Internal duplicate               |   602 |
| Unique valid novel               |   666 |
| Lipinski pass                    |   587 |
| GI-absorption (GIA) ellipse pass |   198 |
| PAINS/BRENK alert-free           |   217 |
| **BBB ellipse pass**             | **0** |

## Key finding: scaffold-level BBB ceiling
Minimum TPSA across all 666 unique candidates is **90.9 Ų** — set by the flavonol core's three
mandatory phenolic/carbonyl oxygens (3-OH, 4-C=O, 5-OH) plus the pyranone ring oxygen. The
BOILED-Egg BBB ellipse requires substantially lower TPSA (roughly <79 Ų at typical WLogP), so
**no R-group decoration on the B-ring can bring a flavonol into the BBB-permeant region** — this
was already flagged during the Phase 4 smoke-test (seed compounds kaempferol/quercetin and the
literature TTBK1/2 inhibitor reference also sit outside the ellipse) and is now confirmed against
the full RL-optimized library of 666 unique, novel, non-PAINS candidates.

Per your decision to keep the BBB gate as a literal, mandatory hard filter (not advisory), the
formal Phase 4 output under that rule is a **null shortlist (0/666)**.

## Reporting-only candidate set (BBB gate excluded)
For reference/reporting, the top candidates passing Lipinski + GI-absorption + alert-free
(i.e., everything except the BBB gate) are dominated by aggressively fluorinated B-ring
decorations reaching Score > 0.85–0.91, TPSA pinned at the 90.9 Ų floor, TTBK1/MAO-B raw docking
scores of 8.5–9.0 / 10–12.8 kcal/mol. 90 molecules pass Lipinski+GIA+alert-free; 62 of those also
score > 0.7 on the RL composite reward.
