# Step 4 — Molecular Docking

## What we did

Docked the 56 shortlisted candidates, the reference inhibitor set, and each receptor's native
ligand into every target site with AutoDock Vina, in triplicate (three independent seeds), across
several deliberately varied receptor preparations.

77 run directories, ~5,000 output files. The naming scheme is the experiment design:

```
04_docking/<RECEPTOR><variant>_<LIGAND SET>_seed<NN>/
```

| Component | Values | Meaning |
|---|---|---|
| `<RECEPTOR>` | `7JXX`, `7Q8Y`, `7Q8V`, `2V5Z`, `2Z5X`, `4BTK` | target structure |
| `<variant>` | *(none)* | crystallographic water shell **retained** |
| | `dry` | all waters removed |
| | `noFAD` | FAD cofactor removed (control) |
| | `brg` | bridging-water variant |
| `<LIGAND SET>` | `candidates_56` | the shortlisted candidate set |
| | `references` | known inhibitors, for anchoring |
| | `native_<LIG>` | the receptor's own crystal ligand (validation, see [step 5](../05_validation/)) |
| | `cand_003` | the lead alone, for focused re-runs |
| `seed<NN>` | `11`, `22`, `33` | independent Vina seeds |

## Why we did it

Docking is the primary screen: it produces the score that ranks candidates and the pose that
everything downstream depends on. Three design choices need justification:

**Why three seeds.** Vina's search is stochastic. A single run's score is one sample, and
seed-to-seed spread on this system runs up to ~0.1 kcal/mol. Consensus across three seeds with a
reported SD makes the number defensible, and makes it obvious when a difference is smaller than
the noise.

**Why wet *and* dry variants.** This is the most important control in the directory. An early
TTBK1-vs-TTBK2 comparison was run with receptors that happened to retain **different numbers of
crystallographic waters** — five in one site, zero in the other. That asymmetry alone produced an
apparent selectivity margin. The `dry` variants exist so the comparison can be made
water-symmetric, and the paired wet/dry runs quantify exactly how much of any margin is a water
artifact rather than chemistry.

**Why a `noFAD` control.** MAO-A and MAO-B carry a covalently linked FAD cofactor that forms one
wall of the substrate cavity. Docking into the site with FAD removed tests whether candidate
scores depend on the cofactor being present — if they do not, the pose is probably not in the real
cavity.

## Reference

- **Engine** — AutoDock Vina
- **Scoring** — default Vina function; Vinardo used as an independent cross-check in
  [step 8](../08_analysis/)
- **Ligand preparation** — Meeko / `obabel`, protonation via `dimorphite-dl` (see
  [step 2](../02_ligands/))
- **Receptor preparation** — [`../scripts/prep_receptor.sh`](../scripts/prep_receptor.sh); grid
  boxes defined per receptor in `03_receptors/<PDB>/box.json`
- **Pocket landmarks** — recorded in
  [`../03_receptors/receptor_characterization.md`](../03_receptors/receptor_characterization.md)

## Inputs and parameters

```
--exhaustiveness 32
--num_modes 9
--seed {11, 22, 33}
--receptor 03_receptors/<PDB>/receptor.pdbqt
--center/--size  from 03_receptors/<PDB>/box.json
```

Driver: [`../scripts/dock.sh`](../scripts/dock.sh). Batch runners:
[`../scripts/run_week2_redock.sh`](../scripts/run_week2_redock.sh) (resumable production re-dock),
[`../scripts/run_water_test.sh`](../scripts/run_water_test.sh),
[`../scripts/run_selectivity.sh`](../scripts/run_selectivity.sh),
[`../scripts/run_minwater_probe.sh`](../scripts/run_minwater_probe.sh).

Two preparation traps are handled in the driver and worth knowing about:

- The ligand list is **not** a blind glob of `02_ligands/pdbqt/*.pdbqt` — that directory holds
  more ligands than any single run should use.
- A truncated SDF is detected by checking for the record terminator, which is the last thing
  written. Vina reports a failed preparation as `BEST None 1000000000.00` — a failure that reads
  like a score unless it is caught.

## Analysis performed

Score collection and consensus: [`../scripts/collect_results.py`](../scripts/collect_results.py).
Interaction profiling: [`../scripts/run_prolif_v2.py`](../scripts/run_prolif_v2.py). Ranking,
consensus and selectivity margins are in [step 8](../08_analysis/); protocol validation is in
[step 5](../05_validation/).

## Final result

- Every receptor passes native-ligand redocking (0.27–1.45 Å) — see [step 5](../05_validation/).
- With the water shell made **symmetric** between the two kinases, TTBK2 is favoured by
  ~1.6 kcal/mol and **0 of 56 candidates** are favourable toward TTBK1.
- On the validated receptor pair the margin is **+0.12 kcal/mol** — but that comparison is
  confounded by the 0-vs-5 water asymmetry, which is why the symmetric number above is the one to
  use.
- The 9IV calibration in [step 5](../05_validation/) shows the protocol carries a
  **0.861 kcal/mol systematic bias toward TTBK2** on a case where experiment says there is no
  preference. Margins below ~0.9 kcal/mol are therefore not interpretable as selectivity.
- `cand_003` is the lead carried into MD: **#8/56 on 7JXX (−8.39)** and **#2/56 on 2V5Z
  (−11.41)**. Both ranks are over the 56 candidates, excluding the native reference ligand that
  shares the file.

  *Corrected 2026-09-27:* this line previously attributed the −11.41 to **7Q8Y**. It belongs to
  **2V5Z**. `consensus_new.csv` has no `7Q8Y`/`cand_003` row at all, and the mislabel mattered
  because it credited the lead's strongest score to the TTBK2 anti-target rather than to MAO-B,
  the intended target — reversing what the number says about selectivity.

The docking stage's honest summary: poses are reliable, rankings are usable, and the
TTBK1/TTBK2 selectivity margin is at or below the protocol's demonstrated bias.

## Relevant files

| Path | Role |
|---|---|
| `<RUN>/<ligand>.log` | Vina log per ligand, all 9 modes with scores |
| `<RUN>/<ligand>_out.pdbqt` | docked poses |
| `<RUN>/<ligand>_out.sdf`, `out.sdf` | poses converted for RMSD / profiling |
| [`../scripts/dock.sh`](../scripts/dock.sh) | the docking call itself |
| [`../scripts/collect_results.py`](../scripts/collect_results.py) | log → score table |
| [`prior_phase/`](prior_phase/) | phase 0–8 docking: baseline scores, interaction fingerprints, rendered poses, and the earlier receptor set (`4NFM`, `6U0K`, `2V60`) |

`prior_phase/` is kept separate because those runs used different receptor structures. Do not pool
its scores with the runs above.
