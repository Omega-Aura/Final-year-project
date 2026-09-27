# Filtering — BOILED-Egg / BBB Permeability Gate

## What we did

Implemented the BOILED-Egg blood–brain-barrier and gastrointestinal-absorption model used as a hard
gate in the filtering cascade and as a live reward component during ligand generation.

## Why we did it

Both target families are intracranial: TTBK1/TTBK2 are tau kinases and MAO-A/MAO-B are brain
monoamine oxidases. A compound that cannot cross the blood–brain barrier is not a candidate no matter
how well it docks, so BBB permeability is applied as a **gate before potency**, not as a
nice-to-have afterwards.

The same model is used in two places on purpose:

- as a **pass/fail gate** in the cascade ([`../08_analysis/`](../08_analysis/)), and
- as a **continuous signed distance** inside the RL reward ([`../generation/`](../generation/)).

A binary gate gives the generator no gradient — everything fails and there is nothing to learn from.
The signed distance tells it *how far* outside the permeable region a molecule sits, which is what
made the BBB-aware campaign work after the docking-only campaign produced zero passing molecules.

## Reference

- **Model** — BOILED-Egg: permeability predicted from the position of a molecule in
  **WLogP vs TPSA** space. The "yolk" region corresponds to predicted BBB permeability, the "white"
  to predicted GI absorption.
- **Descriptors** — RDKit `WLogP` and `TPSA`

## Inputs and parameters

[`boiled_egg_coords.py`](boiled_egg_coords.py) defines the ellipse geometry and the point-in-region
test, and returns a signed distance rather than only a boolean.

Consumers:

| Consumer | Uses it for |
|---|---|
| [`../scripts/filter_cascade.py`](../scripts/filter_cascade.py) | `BBB_pass` / `GIA_pass` columns in the cascade |
| [`../scripts/bbb_score.py`](../scripts/bbb_score.py) | standalone scoring |
| [`../generation/scoring/bbb_score.py`](../generation/scoring/) | `BBB_signed_dist` reward endpoint |

## Analysis performed

Applied to the full 405-molecule pool and to every generated set. Results are the `BBB_pass` and
`GIA_pass` columns in `../08_analysis/filter_cascade_*.csv`.

## Final result

The gate is what defines the project's central early finding: **the native flavonol library has an
absolute BBB ceiling** — the natural glycosides are too large and too polar to pass, which is why
the aglycones and then the deoxygenated scaffold were pursued (manuscript §3.1–3.2).

All 56 shortlisted candidates pass both the BBB and GI gates.

## Relevant files

| Path | Role |
|---|---|
| [`boiled_egg_coords.py`](boiled_egg_coords.py) | the model — region geometry and signed-distance test |
| [`prior_phase/bbb_score.py`](prior_phase/) | the phase 0–8 scoring variant |
| `../08_analysis/filter_cascade_*.csv` | the gate outcomes |
