# system_TTBK2_p2_r2 — velocity replicate of `system_TTBK2_p2`

**This directory contains no new system build.** `complex.prmtop` and `complex.inpcrd` were
copied byte-for-byte from [`../system_TTBK2_p2/`](../system_TTBK2_p2/), so the starting
coordinates, topology, force field, solvation and ion count are *identical* to run 1. Verified at
creation:

```
complex.prmtop  37a7ee4461f2495043e08a24dc7491d1e4bc8ceb2960001e143820e2eae7e61c
complex.inpcrd  6e2621ac63e31dcad189c8ef4453b2c5c96bfc27a06ee784a00aec6679798728
```

Re-check with:

```
sha256sum ../system_TTBK2_p2/complex.prmtop complex.prmtop
sha256sum ../system_TTBK2_p2/complex.inpcrd complex.inpcrd
```

**The only difference is the random velocity seed.** `run_md_system.py` calls
`setVelocitiesToTemperature(300 K)` without an explicit seed, and the LangevinMiddle integrator
likewise draws its own, so re-running the same inputs produces an independent trajectory from the
same starting structure. Same mechanism as `system_TTBK1_r2`, `system_MAOA_p3_r2` and
`system_MAOB_p2_r2`.

There is **no cofactor restraint here**, unlike the MAO replicates: TTBK2 has no FAD, so this run
uses `run_md_system.py` rather than `run_md_restrained.py`. That matters when comparing replicate
spreads across arms — the MAO spreads (0.24–0.26 kcal/mol) are damped by a restraint that this
system does not have, which is exactly why a TTBK2 spread had to be measured rather than borrowed.

## Why this run exists

`system_TTBK2_p2` is the **only on-pose TTBK2 trajectory in the project**. Its ΔG of
−30.95 kcal/mol is the single number the entire TTBK1-vs-TTBK2 comparison rests on, and as of
2026-09-27 it carried **no error bar of its own**. The comparison was reported as "no resolvable
discrimination" on the strength of TTBK1's replicate spread (2.08 kcal/mol) exceeding the
1.45 kcal/mol ΔΔG — but that borrows one arm's uncertainty to judge a difference involving the
other, which is the same arm-asymmetry this project has been bitten by four times
(`../../LOGBOOK.md`, entry K §4.3 of the manuscript).

The two TTBK2 pose-1 runs *do* form a replicate pair, and their spread is 6.89 kcal/mol — but
both of those trajectories left the docked pose (3.72 Å and 5.69 Å), so that spread describes
structures that were never the complex. It is the project's evidence that the reported SEM fails
worst when the pose is least stable; it is **not** a usable error bar for an on-pose energy.

This run supplies the missing quantity: the run-to-run spread of an **on-pose TTBK2** trajectory,
measured rather than inherited. Three outcomes are possible and all three are informative:

| Outcome | What it means |
|---|---|
| Spread ≪ 1.45 | The TTBK ΔΔG may be larger than the noise after all — the "not supportable" verdict would need revisiting against the calibration floor (~1.3 kcal/mol), which would then become the binding constraint. |
| Spread ≈ 2 (like TTBK1) | Confirms the published verdict, now with a symmetric error bar on both arms instead of one. |
| Spread ≫ 2 | Strengthens it further, and extends the pose-stability/SEM relationship with a fifth data point. |

No build inputs (`build.leap`, `7Q8Y_prot.pdb`, ligand parameters) are copied here on purpose:
`tleap` was **not** re-run, and shipping its inputs would imply it had been. They are in
[`../system_TTBK2_p2/`](../system_TTBK2_p2/).
