# Ligand Generation — REINVENT4 Staged Reinforcement Learning (campaign v3)

## What we did

Ran the final generative campaign: REINVENT4 LibInvent decorating a 7-deoxyflavonol scaffold under
a **four-component reward** that optimises docking score against both targets *simultaneously* with
synthetic accessibility and BBB permeability.

This is campaign **v3**, the one behind the final candidate set. Earlier campaigns are in
[`../00_library/reinvent4_output/`](../00_library/) (`campaign1_first_failed`, `campaign2_v2`).

## Why we did it

**Why generate at all.** The native flavonol library has an absolute BBB ceiling — the natural
compounds cannot reach an intracranial target. Scaffold deoxygenation plus generative decoration is
what converts a null result into candidates (manuscript §3.1–3.2).

**Why a BBB-aware reward.** Campaign 1 used a docking-only reward and produced **zero** molecules
passing the literal BBB gate. Optimising potency alone reliably produces large, polar, potent
molecules that cannot cross into the CNS. `BBB_signed_dist` puts permeability inside the objective
rather than applying it as a post-hoc filter that everything fails.

**Why dual-target in one reward.** The study's premise is a single compound engaging both target
families. Scoring TTBK1 and MAO-B with equal weight under a **geometric mean** means a molecule
cannot win by excelling at one target and failing the other — the geometric mean collapses when any
component is near zero. An arithmetic mean would have permitted exactly that trade-off.

**Why a diversity filter.** Without `IdenticalMurckoScaffold` bucketing, staged learning converges
onto one scaffold and reports a high score for what is effectively a single molecule.

## Reference

- **Generator** — REINVENT4 LibInvent, `staged_learning`
- **Prior / agent** — `REINVENT4/priors/libinvent.prior`
- **Learning strategy** — DAP, σ = 128, rate 1e-4
- **Diversity filter** — `IdenticalMurckoScaffold`, bucket size 25, minscore 0.4
- **Docking inside the reward** — AutoDock Vina via
  [`scoring/dock_score.py`](scoring/dock_score.py)
- **BBB component** — BOILED-Egg signed distance,
  [`scoring/bbb_score.py`](scoring/bbb_score.py)

## Inputs and parameters

Config: [`rl_staged_learning_v3_local.toml`](rl_staged_learning_v3_local.toml) (resolved copy in
`_rl_staged_learning_v3.json`).

```
scaffold      generation/flavonol_7deoxy_scaffold.smi
batch_size    32
randomize_smiles  true
stage 1       min_steps 15, max_steps 20, max_score 0.9
aggregation   geometric_mean
```

Reward components:

| Component | Weight | Detail |
|---|---|---|
| `TTBK1_dock` | 0.35 | Vina into `03_receptors/7JXX/receptor.pdbqt`, box 18³ at (178.93, 19.64, 46.21) |
| `MAOB_dock` | 0.35 | Vina into the MAO-B site |
| `SA_score` | 0.15 | synthetic accessibility |
| `BBB_signed_dist` | 0.15 | BOILED-Egg signed distance |

Docking inside the reward loop uses reduced settings (`--exhaustiveness 4 --num_modes 5`) for
throughput. **Those scores are not the reported scores** — every candidate is re-docked at
`--exhaustiveness 32 --num_modes 9` across three seeds in [step 4](../04_docking/).

## Analysis performed

Generated molecules → filtering cascade and consensus re-docking
([step 8](../08_analysis/), [step 4](../04_docking/)). Reward progression is in
`tb_logs_0/` (TensorBoard) and the `rl_v3_local*.csv` summaries.

## Final result

The v3 campaign produced the pool that became
[`../01_smiles/candidates_56.csv`](../01_smiles/), including the lead `cand_003`. 49 generated
molecules out-scored the lead **on the RL objective**
([`rl_v3_beats_cand003.csv`](rl_v3_beats_cand003.csv)) — a statement about the reward function, not
about binding; those were re-docked and re-filtered before any claim was made.

A convergent result worth noting: top candidates independently acquired **CF3 / CHF2**
substituents. The generator found that pattern; it was not designed in.

## Relevant files

| Path | Role |
|---|---|
| [`rl_staged_learning_v3_local.toml`](rl_staged_learning_v3_local.toml) | **the campaign definition** |
| `rl_stage1_v3.chkpt` | trained agent checkpoint (95 MB) — retained; resumes or re-samples the campaign |
| [`rl_v3_local_1.csv`](rl_v3_local_1.csv) | per-step sampled molecules with scores |
| [`rl_v3_dedup.csv`](rl_v3_dedup.csv) | deduplicated output |
| [`rl_v3_beats_cand003.csv`](rl_v3_beats_cand003.csv) | molecules out-scoring the lead on the RL objective |
| `rl_v3_local.log` | run log |
| [`scoring/dock_score.py`](scoring/dock_score.py), [`scoring/bbb_score.py`](scoring/bbb_score.py) | reward endpoints |
| `flavonol_7deoxy_scaffold.smi` | the decorated scaffold |
| `tb_logs_0/` | TensorBoard reward progression |

`REINVENT4/` at the repository root is a local third-party clone (untracked, not vendored). The
config references `REINVENT4/priors/libinvent.prior`, so the clone must be present to re-run.
