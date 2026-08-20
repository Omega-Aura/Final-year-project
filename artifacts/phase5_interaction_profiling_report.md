# Phase 5 — Interaction Profiling & Pose Rendering

## Scope
Re-docked the three lead flavonoid candidates from the Phase 6 consensus ranking against both
validated on-target receptors (TTBK1 = PDB 4NFM, MAO-B = PDB 2V5Z) with AutoDock Vina at high
exhaustiveness, reproducing the Phase 6 consensus binding energies. Best poses were converted to
RDKit molecules using meeko's `PDBQTMolecule`/`RDKitMolCreate` (bond topology taken from the SMILES
embedded in the PDBQT REMARK lines, not from nearest-neighbor atom matching — see Methods note),
profiled with ProLIF, and rendered in PyMOL.

Candidates:
- `cand_primary_CF3` — CF3-substituted 7-deoxyflavonol regioisomer (Phase 6 primary recommended pick;
  most consistent under both Vina and Vinardo scoring)
- `cand_benzylMe_A`, `cand_benzylMe_B` — benzyl/methyl regioisomers that ranked #1 and #3 under the
  Vina+Vinardo combined average (flagged in the Phase 6 report as possibly an averaging artifact)

## Interaction fingerprints (ProLIF)

| Candidate | Receptor | # contacts | Key residues |
|---|---|---|---|
| cand_primary_CF3 | TTBK1 (4NFM) | 10 | ILE40, ILE48, GLN110, GLY111, ASN113 (H-bond acceptor), LEU175 |
| cand_primary_CF3 | MAO-B (2V5Z) | 14 | PHE168, LEU171, CYS172, ILE198 (H-bond donor), ILE199, GLN206, ILE316, TYR326, PHE343, TYR435 |
| cand_benzylMe_A | TTBK1 (4NFM) | 12 | ILE40, ILE48, GLN110, GLY111, ASN113, ASN159, LEU175, ASP176 |
| cand_benzylMe_A | MAO-B (2V5Z) | 19 | TYR60, LEU164, PHE168, LEU171, CYS172, ILE198, ILE199, GLN206, ILE316, TYR326, PHE343, TYR435 |
| cand_benzylMe_B | TTBK1 (4NFM) | 12 | GLY41, ILE48, GLN110, LYS156 (H-bond donor), SER158, LEU175, ASP176 |
| cand_benzylMe_B | MAO-B (2V5Z) | 17 | PHE168, LEU171, CYS172, TYR188, ILE199, GLN206, ILE316, TYR326, PHE343, TYR398, TYR435 |

Observations:
- All three candidates engage the MAO-B aromatic cage (TYR326/PHE343/TYR435, adjacent to the FAD
  cofactor site) that the pipeline has targeted since receptor-prep validation. The two benzyl/methyl
  regioisomers make more total contacts at both targets (12/17-19) than the CF3 candidate (10/14),
  consistent with their larger, more hydrophobic substituents filling more of each pocket.
- Only `cand_benzylMe_B` registers a scored contact with TYR398, MAO-B's catalytic tyrosine; the
  other two candidates pack near it in the rendered poses but ProLIF does not score it as a formal
  contact for those poses.
- At TTBK1, all three candidates contact the ASN113/GLN110/LEU175 region of the ATP pocket, with each
  candidate forming a distinct single hydrogen bond (ASN113 acceptor for CF3, LYS156 donor for
  benzylMe_B) — benzylMe_A shows no scored H-bond at TTBK1 in this profiling pass.
- These fingerprints do not by themselves establish on-target vs. anti-target selectivity; that
  question was already addressed (largely negatively) in the Phase 6b counter-screening report, whose
  caveats (unvalidated anti-target receptor preparations) still apply.

## Rendered poses
Six pose figures were generated (one per candidate x target), each showing the receptor as a
transparent cartoon, the ligand as yellow sticks, and the ProLIF-identified contact residues as
teal sticks with labels:

- `pose_primary_TTBK1.png`, `pose_primary_MAOB.png`
- `pose_benzylMeA_TTBK1.png`, `pose_benzylMeA_MAOB.png`
- `pose_benzylMeB_TTBK1.png`, `pose_benzylMeB_MAOB.png`

## Methods note — pose-reconstruction bug found and fixed
The first rendering pass produced visibly wrong ligand connectivity (a tangled, non-chemical
wireframe). Root cause: the existing `pdbqt_utils.build_posed_mol` helper matched pose atoms to the
RDKit template by nearest-neighbor distance, but compared the docked pose's pocket coordinates
against the template's original, undocked conformer coordinates — two unrelated coordinate frames,
so the atom-identity assignment was effectively arbitrary. This was fixed by switching to meeko's own
`PDBQTMolecule` / `RDKitMolCreate.from_pdbqt_mol`, which rebuilds each pose directly from the SMILES
meeko embeds in the PDBQT file's REMARK lines, giving one correctly-bonded conformer per pose. The
fix was verified by round-tripping the canonical SMILES of a reconstructed pose against the intended
candidate structure and confirming sane (in-pocket) coordinates before regenerating the ProLIF
fingerprints and renders. The corrected fingerprints were consistent with the pre-fix run (the bug
affected only 3D geometry used for rendering, not the profiling), but all downstream renders and the
final fingerprint table in this report use the corrected poses.

## Status
Plan step "Render pose and interaction figures" (Phase 5 — Precision validation & interaction
profiling) is complete. MD-dependent phases (system building, production MD, MM-PBSA) remain blocked
pending an attached GPU compute target.
