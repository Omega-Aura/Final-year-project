# system_MAOA_p3_r2 — velocity replicate of `system_MAOA_p3`

**This directory contains no new system build.** `complex.prmtop` and `complex.inpcrd` were
copied byte-for-byte from [`../system_MAOA_p3/`](../system_MAOA_p3/), so the starting coordinates, topology, force
field, solvation, ion count and FAD restraint targets are *identical* to run 1. Verify with:

```
sha256sum ../system_MAOA_p3/complex.prmtop complex.prmtop
sha256sum ../system_MAOA_p3/complex.inpcrd complex.inpcrd
```

**The only difference is the random velocity seed.** `run_md_restrained.py` calls
`setVelocitiesToTemperature(300 K)` without an explicit seed, and the LangevinMiddle integrator
likewise draws its own, so re-running the same inputs produces an independent trajectory from the
same starting structure. This is the same way `system_TTBK1_r2` was produced.

**Why this run exists.** Every MAO conclusion as of 2026-09-27 rests on one trajectory per pose.
The MAO ΔΔG direction (MAO-B favoured across all six on-pose pairings) is robust, but its
*magnitude* had no error bar, because the only measured replicate spreads in this project come
from the TTBK systems: 2.08 kcal/mol for TTBK1 and 6.90 for TTBK2, both 10–30x the SEM that
`MMPBSA.py` reports. Pose spread (0.44–1.13 kcal/mol here) is a different and demonstrably
smaller source of variance and cannot substitute for it. This replicate measures the MAO arms'
own spread directly. See [`../../07_mmgbsa/README.md`](../../07_mmgbsa/README.md).

No build inputs (`build.leap`, protein PDB, ligand/FAD parameters) are copied here on purpose:
`tleap` was **not** re-run, and shipping its inputs would imply it had been. They are in
[`../system_MAOA_p3/`](../system_MAOA_p3/).
