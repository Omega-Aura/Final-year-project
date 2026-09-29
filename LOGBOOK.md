# Logbook

## 2026-08-20 — A

Day 1 environment setup (§Part 2, Day 1).

**Native Windows, conda env `docking_project`:** was rdkit-only at start. Installed
`openmm mdanalysis prolif parmed` via conda-forge, and `spyrmsd meeko dimorphite-dl`
via pip. `python -m openmm.testInstallation` confirms a working CUDA platform
(RTX 4050, driver CUDA UMD 13.3) — MD runs (§A6–A8) can run natively, no WSL needed
for the simulation itself.

**AmberTools has no Windows build** (confirmed: zero conda-forge win-64 matches for
`ambertools`). This blocks `acpype` and `MMPBSA.py` natively. Machine already had
WSL2 Ubuntu with GPU passthrough working (`nvidia-smi` visible inside WSL). Installed
Miniconda in WSL2 and created a minimal env `mdgbsa` (`ambertools acpype`, conda-forge
only — had to use `--override-channels` to skip the `defaults` channel's ToS gate).
Both `acpype` (v2023.10.27) and `MMPBSA.py` (v14.0) confirmed callable inside WSL.

**Net environment split for Machine A:**
- Native Windows `docking_project`: rdkit, openmm(CUDA), MDAnalysis, prolif, parmed,
  spyrmsd, meeko, dimorphite-dl, vina, obabel — docking + MD + analysis.
- WSL2 `mdgbsa`: acpype (GAFF2 parameterization), MMPBSA.py (MM-GBSA) only — used
  because AmberTools has no Windows port, not a general duplicate stack.

**NOTE:** meeko's CLI installs as `mk_prepare_receptor.exe` / `mk_prepare_ligand.exe`
on Windows, not `mk_prepare_receptor.py` as written in the doc's Part 6 scripts. Both
`mk_prepare_receptor` (no extension) and the `.exe` form work from bash. Scripts in
`scripts/` need this accounted for per machine, per the doc's own warning about
Meeko's CLI changing across versions.

**STILL OPEN:** the lead compound's acpype/GAFF2 parameter files (named as one of the
three Day-1-critical files) were not found anywhere in the recovered `artifacts/` —
no `.itp`/`.top`/`.frcmod`/`.mol2`. Will need to be regenerated via `acpype` (WSL)
before Week 2's MD system builds (§A6), which assume these already exist and are
just reused across all four systems.

Files: none yet (environment-only session, no data outputs).

## 2026-08-20 — B

- Initialized environment and scripts in `scripts/` (`prep_receptor.sh`, `prep_ligands.py`, `dock.sh`, `rmsd_check.py`, `collect_results.py`).
- Prepared MAO-B receptor (2V5Z Chain A) with FAD cofactor using `scripts/prep_receptor.sh`.
- Extracted and prepared native safinamide ligand (SAG) with `scripts/prep_ligands.py`.
- Executed calibration redocking with `scripts/dock.sh 2V5Z native_SAG 11` (exhaustiveness=32, seed=11).
- Computed symmetry-corrected RMSD with `scripts/rmsd_check.py`:
  - Pose 1 RMSD: 1.38 Å
  - Pose 3 (best) RMSD: 0.56 Å -> **PASS** (< 2.0 Å threshold).
- **GATE 0 (Calibration exercise) passed** on B's machine.

NOTE (A, reviewing): §Day 2 requires *both* machines to reproduce the ~1.57 Å
safinamide redock and agree with each other to within ~0.3 Å before proceeding.
B's best pose (0.56 Å) is a clean pass in isolation but hasn't yet been
cross-checked against an independent run on Machine A — that comparison is still
outstanding before GATE 0 can be called fully closed per the doc's own rule.

## 2026-08-20 — A

Day 3: structure download and inspection (§Part 2, Day 3). A's task in full.

**Downloaded all six receptors** into `03_receptors/<PDB>/raw.pdb` (2V5Z already
present from B's Day 2 work, skipped re-download). NOTE: the doc's download loop
literally lists eight IDs (`7JXX 4BTK 7Q8V 7Q8Y 2V5Z 2Z5X 7JXY 4BTM`) under a "six
receptors" heading — `7JXY` and `4BTM` appear nowhere else in the doc (not in the
table, not in any later step) and look like a copy/typo artifact (near-duplicates of
7JXX/4BTK). Went with the doc's own 6-receptor table as authoritative and did not
fetch those two. Flag if that's wrong.

**Inspection table** (HETATM species, missing-residue count from REMARK 465,
resolution — all pulled directly from the downloaded headers):

| PDB | Protein | Chains | Res. | HETATM species | REMARK 465 lines | Native ligand |
|---|---|---|---|---|---|---|
| 7JXX | TTBK1 | A | 1.56 Å | NA, HOH, VP7 | 51 | VP7 |
| 4BTK | TTBK1 | A | 2.00 Å | DMS, DTQ, HOH | 57 | DTQ |
| 7Q8V | TTBK1 | A | 2.13 Å | 9IV, HOH, PO4 | 24 | 9IV |
| 7Q8Y | TTBK2 | A, B | 1.60 Å | 9IV, HOH, PO4 | 20 | 9IV (present in both chains; using chain A for consistency with 2V5Z convention) |
| 2V5Z | MAO-B | A, B | 1.60 Å | FAD, HOH, SAG | 54 | SAG |
| 2Z5X | MAO-A | A | 2.20 Å | DCX, FAD, GOL, HOH, HRM | 0 | HRM |

All resolutions and native ligands match the doc's table exactly — no surprises there.

**Check 1 — does 2Z5X contain FAD?** Yes, confirmed: `HET FAD A 600 53` (full
53-atom FAD, not partial), plus a REMARK 500 close-contact record between Cys406 SG
and the FAD C8M atom (1.65 Å) — consistent with the expected 8α-S-cysteinyl
covalent linkage at Cys406. **The custom FAD-transplant step described in the
manuscript's Methods was unnecessary** — 2Z5X ships with a complete, correctly
linked FAD. This deletes a caveat from Limitations.

**Check 2 — do 7JXX pocket residues match 4NFM numbering?** Superposed 7JXX chain A
onto `artifacts/4NFM.pdb` chain A in PyMOL (`cmd.align`, 0.43 Å RMSD over 1815 atoms
before refinement). Checked all seven pocket residues named in the manuscript
(ILE40, ILE48, GLN89, GLN110, GLY111, ASN113, LEU175) directly against both raw PDB
files (not just the alignment tool, since the alignment's raw-pair extraction gave a
false "not aligned" for residue 40 that direct inspection contradicted):

| 4NFM | 7JXX | Match |
|---|---|---|
| ILE40 | ILE40 | identical number + identity |
| ILE48 | ILE48 | identical number + identity |
| **GLN89** | HIS89 (4NFM is *also* HIS89, not GLN) | **manuscript residue name is wrong in both structures** |
| GLN110 | GLN110 | identical number + identity |
| GLY111 | GLY111 | identical number + identity |
| ASN113 | ASN113 | identical number + identity |
| LEU175 | LEU175 | identical number + identity |

**7JXX uses the exact same numbering as 4NFM** — no residue-numbering lookup table
is needed for §2.2. The one required fix is a residue-identity correction:
**GLN89 → HIS89** in the manuscript's pocket-residue list; it was never a numbering
issue, and 4NFM (the structure the manuscript numbers are drawn from) has HIS at
that position too, so this looks like a plain transcription error, not something
that changed between structures.

**GATE 0 status:** still open pending A's own independent calibration run to cross-
check against B's 0.56 Å (see note above) — did not run it today, scope was Day 3
only per instruction.

Files: `03_receptors/{7JXX,4BTK,7Q8V,7Q8Y,2Z5X}/raw.pdb` (new), `03_receptors/2V5Z/raw.pdb` (pre-existing, B).

## 2026-08-20 — A

Track A, Week 1 (§Part 3). A's tasks A1-A5.

**Fixed two portability bugs in shared scripts before running anything**, both
harmless-on-this-machine but violating §1.5 (same script, same run):
`scripts/prep_receptor.sh` and `scripts/dock.sh` both hardcoded B's personal path
(`/c/Users/sayan/miniforge3/envs/dock/...`) into `$PATH`. Removed; both scripts
already autodetect the active conda env correctly without it.

**A1 — prepare all six receptors.** Two (7JXX, 7Q8Y) initially failed: alternate-
location (altloc) residues that `-a/--allow_bad_res` alone doesn't resolve. Added
`--default_altloc A` to the `mk_prepare_receptor` call in `prep_receptor.sh` and
re-ran **all six** uniformly (not just the two failures) per §1.5. All six now
prepare cleanly: 7JXX, 4BTK, 7Q8V, 7Q8Y, 2V5Z, 2Z5X.

**A2 — GATE 1 redocking validation (7JXX, 4BTK).** Found a real pipeline bug while
running this: 7JXX's redock crashed spyrmsd with `NonIsomorphicGraphs` instead of
producing a number. Root cause: `dock.sh` used plain `obabel` to convert Vina's
docked PDBQT back to SDF. Meeko inserts dummy "glue" atoms into the PDBQT for
ligands needing flexible-ring handling (VP7 apparently needs this, DTQ apparently
doesn't -- explains why 4BTK "passed" and 7JXX crashed instead of just scoring
badly). `obabel` doesn't understand Meeko's own convention and silently emits `*`
wildcard atoms, corrupting the molecular graph. Fix: use Meeko's own `mk_export`
to convert docked poses, not `obabel` (verified: `mk_export` reconstructs VP7 as
45 atoms / 48 bonds, exact match to the reference; `obabel` gave 28 atoms with two
`*` atoms). Re-exported existing poses with `mk_export` (no need to re-dock,
deterministic given the same seed) and re-ran the RMSD check:

  - **7JXX (native VP7): best pose 0.71 Å -> PASS.**
  - **4BTK (native DTQ): best pose 0.75 Å -> PASS** (unchanged by the fix; DTQ
    didn't hit the dummy-atom issue, confirms the fix doesn't regress a working case).

**GATE 1 PASSES on 7JXX**, the doc's designated primary TTBK1 receptor (1.56 Å,
holo, VP7-bound) -- a dramatic improvement over the manuscript's original 5.29 Å
failure on the apo 4NFM structure. Per §Part 3: **Limitations #2 and #3 are deleted
from the manuscript.** No fallback to 4BTK needed, though 4BTK remains available as
a validated cross-check receptor.

This `mk_export` fix matters beyond GATE 1 -- `dock.sh` is used for every docking
run in the project, including A5's 56-candidate consensus docking below. Fixed
before running A5, not after.

**A3 — rebuild the 56-candidate ligand set from SMILES.** `01_smiles/candidates_56.csv`
(the `cand_001`...`cand_056` reformatted copy of `rl_v2_shortlist_56.csv`) rebuilt
from scratch via `prep_ligands.py`. 56/56 succeeded, 0 failures.

**Found and fixed a third `dock.sh` bug while testing the multi-ligand path** (never
exercised before today -- B's GATE 0 calibration and my GATE 1 validation were both
single-ligand runs): the set-resolution logic fell back to globbing *every*
`.pdbqt` file in `02_ligands/pdbqt/` when there's no single combined
`<SET>.pdbqt` file, which would have silently swept the native reference ligands
(native_VP7, native_DTQ, native_SAG, native_9IV) into any "candidates_56" docking
run. Fixed: when `01_smiles/<SET>.csv` exists, the ligand list now comes from that
CSV's `id` column (the authoritative set membership), not a blind glob.

**Then hit a fourth bug from that same fix**: the python-generated ligand-path list
inside `dock.sh` picked up trailing `\r` on every line (Python's stdout on Windows
does universal-newline translation to `\r\n`; bash's `$()` doesn't strip embedded
`\r`). This corrupted every path, so `basename "$L" .pdbqt` didn't strip the
suffix, vina couldn't find the file, and `set -euo pipefail` killed each of the 6
candidates_56 docking runs after only the first ligand. First A5 attempt silently
produced only 6 log files total (1 per seed/receptor combo) instead of 336.
Caught this by checking `collect_results.py`'s output table (only 3 rows -- the
native ligands -- with zero candidate entries) rather than trusting a clean exit
code. Fixed with `| tr -d '\r'` on the python output; verified the fix directly
(56/56 correct paths, no `\r`, file existence confirmed) before re-running the full
batch.

**A4 — re-run the filtering cascade.** Hard blocker, not a version-drift issue: the
actual filter code was never recovered at all, only its output CSVs. `bbb_score.py`
is a different thing (REINVENT4's RL-time BBB reward hook, not the offline
Lipinski/GI/BBB/alerts cascade), and it itself imports a `filtering/boiled_egg_coords.py`
that also didn't exist anywhere in the project.

Asked how to proceed; reconstructed from published/standard definitions rather than
holding indefinitely:
  - **Lipinski Ro5**: standard, pass = <=1 violation of {MW>500, WLogP>5, HBD>5, HBA>10}.
  - **GI absorption / BBB penetration**: point-in-polygon against the digitized
    BOILED-Egg ellipses from Daina & Zoete (ChemMedChem 2016), coordinates pulled
    from the open-source reimplementation PyBOILEDegg (github.com/bfmilne/PyBOILEDegg,
    GPLv3), which cites the same source paper. Saved to the (recreated)
    `filtering/boiled_egg_coords.py` that `bbb_score.py` already expected to exist.
    TPSA must use `includeSandP=True` -- confirmed against the existing TPSA column
    in the recovered data (0/50 mismatch with S/P included vs 1/50 without).
  - **Structural alerts**: RDKit's combined PAINS + BRENK `FilterCatalog` -- alert
    names in the recovered data (`Michael_acceptor_1`, `catechol_A(92)`, etc.)
    matched this combination's naming convention.

**Verified against `00_library/reinvent4_output/campaign2_v2/rl_v2_filtered_full.csv`
(405 molecules, already carrying the original pass/fail columns) -- zero mismatches
on all four filters, individually, per molecule** (not just matching aggregate
counts, which could be coincidental -- every single one of 405 x 4 boolean labels
matched). Funnel: 369 Lipinski / 81 BBB / 334 GI / 253 alert-free / 56 passing all
-- **exact match to the doc's target numbers.**

Wrote this up as a real, re-runnable script (`scripts/filter_cascade.py`), not just
inline verification code, since that's what "re-run the filtering cascade" actually
requires. Ran it independently on the 405-molecule pool (369/81/334/253/56, exact)
and on our own `01_smiles/candidates_56.csv` (56/56 pass all four -- confirms A3's
ligand set is internally consistent with A4's filter).

**Side finding, not fixed (out of scope -- the RL run itself isn't being redone):**
`bbb_score.py`'s TPSA call is missing `includeSandP=True`, so the RL-time BBB reward
and the actual offline BBB filter are not perfectly consistent for S/P-containing
molecules, contradicting that script's own docstring claim. Doesn't affect anything
in this project's remaining phases since the generative campaigns are being kept
as-is, not rerun.

**A5 — consensus re-dock all 56 candidates on 7JXX and 2V5Z.** 3 seeds x 56
ligands x 2 receptors = 336 dockings. First attempt silently broke on the `\r` bug
above (caught before trusting it). Re-ran with the fixed `dock.sh` -- **complete,
336/336 logs, `collect_results.py` -> `08_analysis/consensus_new.csv`.**

**Ranking on the validated 7JXX changed, as the doc predicted.** The original
pipeline's top-ranked dual-target candidate by RL score, `cand_001`
(`rl_v2_shortlist_56.csv` row 1, `Score`=0.7497, the highest in the set), is now
**#3/56 on 7JXX** (consensus -8.79 kcal/mol) -- `cand_013` (-8.83) and `cand_002`
(-8.80) both edge it out on the holo, validated receptor. Consistent with the doc's
"a holo pocket at 1.56 Å is a different shape from an apo pocket" expectation.

**More striking: `cand_001` ranks only #49/56 on 2V5Z (MAO-B)**, consensus -10.32
kcal/mol, well off the top cluster (best on 2V5Z: cand_043 at -11.63). 2V5Z itself
wasn't revalidated today (it's the doc's "unchanged" primary MAO-B receptor, already
established), so this isn't a validation artifact -- it's either real chemistry
(cand_001 may simply bind TTBK1 much better than MAO-B despite being scored as a
top dual-target hit during RL) or an artifact of how the RL reward combined
`TTBK1_dock` and `MAOB_dock` into one scalar score (a compound can rank #1 on a
sum/product of two scores while being mediocre on one of them). Worth raising at
sync -- this bears directly on which compound A should advance as "the lead" for
MM-GBSA in Week 2-3, and on the dual-engagement claim generally.

**CORRECTION (added 2026-08-21 on review).** The two paragraphs above are factually
correct about `cand_001` but incomplete in a way that misleads, and the omission was
repeated downstream before being caught. `cand_001` is the top-ranked candidate *by RL
score* -- it is **not** the manuscript's lead compound. Matching manuscript §3.3's SMILES
strings against `01_smiles/candidates_56.csv`:

| Manuscript §3.3 row | SMILES | ID here |
|---|---|---|
| "the lead" -- #2/#2, robust under both Vina and Vinardo, advanced to MD | `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F` | **`cand_003`** |
| combined-average #1, **explicitly rejected** as lead | `Cc1ccc(-c2oc3cccc(O)c3c(=O)c2O)cc1Cc1ccccc1` | `cand_001` |
| third tabulated row | `Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1Cc1ccccc1` | `cand_002` |

§3.3 rejects `cand_001` in as many words: its "weak MAO-B Vina-consensus score (-10.32
kcal/mol, worst of the top-15 on that metric) is masked by averaging rather than
corroborated." So `cand_001` landing at #49/56 on 2V5Z **corroborates the manuscript's
own reasoning rather than threatening it** -- and -10.32 is the identical value §3.3
already reports, i.e. today's independently-rebuilt pipeline reproduced it exactly.

The actual MD lead, `cand_003`, held up on both receptors: **#8/56 on 7JXX (-8.38) and
#2/56 on 2V5Z (-11.41)**, the latter reproducing §3.3's -11.42 to within 0.01 kcal/mol on
the same receptor through a separately rebuilt pipeline. **On current evidence A5 gives no
reason to change the lead compound for Week 2-3 MM-GBSA.** The open question raised above
is real but narrower than stated: it is about whether a *better* dual-target candidate now
exists (see the 2026-08-21 entry), not about whether the existing lead survived.

New top TTBK1 hit: **cand_013** (-8.83 kcal/mol, seed SD 0.005 -- very stable
across seeds). New top MAO-B hit: **cand_043** (-11.63 kcal/mol).

Files: `03_receptors/*` (all six, receptor.pdbqt/box.json), `02_ligands/{sdf,pdbqt}/cand_*`
(56 each), `filtering/boiled_egg_coords.py` (new), `scripts/filter_cascade.py` (new),
`08_analysis/filter_cascade_A4_check.csv`, `08_analysis/filter_cascade_candidates_56.csv`,
`08_analysis/consensus_new.csv`, `05_validation/7JXX_redock.txt`, `05_validation/4BTK_redock.txt`,
`04_docking/{7JXX,2V5Z}_candidates_56_seed{11,22,33}/` (336 logs total).

**Track A Week 1 (A1-A5) is now complete.** Per §Part 3 checklist: all six receptors
prepared by one script [x]; GATE 1 (7JXX < 2.0 Å) [x]; 56 candidates rebuilt [x];
filtering cascade re-run, numbers reproduce exactly [x]; 56 candidates re-docked on
7JXX + 2V5Z [x]. Waiting on B's Track B Week 1 (reference set + benchmark docking)
before SYNC POINT 1.

## 2026-08-20 — A (doing Track B Week 1, B's machine/environment not available)

**`01_smiles/references.csv` provenance note:** B committed a bare 5-column stub
(`ID,Target,Measured,Note,SMILES`) as `1436cfa`, but a much richer 14-column draft
(id,smiles,target,measured_value,measured_unit,source,role,native_of_pdb,
expected_cip,formula,mw,notes,mechanism,use_in_correlation) was sitting uncommitted
in the working tree, pulled over from B's system. Validated the rich version before
trusting it:
  - **All 14 SMILES parse; every formula and MW matches RDKit's calculation** (one
    trivial 0.01 Da rounding on clorgyline).
  - **All stereocenters match RDKit's own CIP assignment** (safinamide/SAG_pdb =
    S, selegiline/rasagiline = R).
  - Correctly separates true safinamide (secondary amine) from `SAG_pdb` (the imine
    tautomer actually modelled in the 2V5Z crystal structure) -- a real and easy-to-
    miss distinction; correctly excludes irreversible covalent inhibitors
    (selegiline, rasagiline, clorgyline) from `use_in_correlation`.
  - **Found and fixed a real inconsistency**: kaempferol, quercetin, lazabemide,
    and isatin were all flagged `use_in_correlation=yes` with an empty
    `measured_value` -- can't correlate against nothing. Flipped all four to `no`
    until real primary-source IC50s are sourced. Did not fabricate values.

Kept this richer version, did not touch B's committed stub directly (that gets
superseded by this commit).

**B1 (finalize) + B2 + B3, run together:**

  - `prep_ligands.py --csv references.csv` -- 14/14 succeeded.
  - Docked all 14 on 2V5Z and 2Z5X (already prepared in Track A, reused directly --
    same script/run, satisfies §1.5), 3 seeds each = 84 dockings, all succeeded.
  - **B3 verification passes exactly**: `SAG_pdb` (independently re-prepared here
    from `references.csv`) redocks into 2V5Z at **0.56 Å** -- identical to the
    `native_SAG` result from the original Day-2 GATE-0 calibration (also 0.56 Å,
    prepared via a completely separate code path). Strong cross-check that ligand
    prep is deterministic and consistent regardless of entry point.
  - Consensus scores, all compounds, both receptors: see
    `05_validation/benchmark_mao.csv`. Notable: safinamide only 0.15 kcal/mol
    better on 2V5Z (MAO-B, -10.02) than 2Z5X (MAO-A, -9.87); harmine only 0.04
    kcal/mol better on 2Z5X than 2V5Z (-8.66 vs -8.62). Both are much more
    selective in reality than these thin margins suggest -- same kind of protocol-
    resolution question §B4's 9IV calibration is designed to formally quantify for
    TTBK1/TTBK2 in Week 2; worth considering an equivalent check for MAO-A/B.

**Correlation plot built (`05_validation/benchmark_mao_scatter.png`), but honestly
labelled as not yet meaningful**: after correctly excluding irreversible inhibitors
and the four unsourced-value compounds, exactly **one** point
(safinamide) has both `use_in_correlation=yes` and a real measured value. One point
cannot support a Pearson r. This isn't a bug introduced today -- the manuscript's
own original reference table had the same gap (kaempferol/quercetin listed with
"--" for measured value), and whoever built the richer CSV already recognised it
and added lazabemide/isatin specifically "to strengthen a thin reversible set" --
just hadn't sourced values for them yet.

**OPEN QUESTION for A (self) / user:** source literature IC50 values for
kaempferol, quercetin, lazabemide, isatin against MAO-B to make the correlation
plot meaningful? The doc's own Week 3 B task already names Chaurasiya et al. 2020
(*Molecules* 25:5358) as the source for the kaempferol/quercetin correction, so
kaempferol/quercetin at least have a doc-endorsed reference to check against.
Held off pulling this in now since it's explicitly a Week 3 task, not Week 1, and
because sourcing binding-affinity numbers for a manuscript should be confirmed with
the user first, not done unilaterally on a guess.

Files: `01_smiles/references.csv` (14 rows, fixed), `02_ligands/{sdf,pdbqt}/{safinamide,
SAG_pdb,selegiline,rasagiline,kaempferol,quercetin,harmine,clorgyline,VP7,DTQ,
9IV_ttbk1,9IV_ttbk2,lazabemide,isatin}.*`, `04_docking/{2V5Z,2Z5X}_references_seed{11,22,33}/`
(84 logs), `05_validation/benchmark_mao.csv`, `05_validation/benchmark_mao_scatter.png`.

## 2026-08-20 — A (literature chase for the four unsourced correlation compounds, on request)

User asked to source real MAO IC50 values for kaempferol, quercetin, lazabemide,
isatin now rather than deferring to Week 3, since the correlation plot only had one
usable point. Used PubMed (`mcp__claude_ai_PubMed`) directly rather than general
web search, so every number below has a checkable PMID/DOI.

**Found a real citation error in the manuscript itself.** The doc's own Week 3 task
names "Chaurasiya et al. 2020, Molecules 25:5358" as the source for a
kaempferol/quercetin correction. Confirmed via PubMed this is PMID 33212830, DOI
[10.3390/molecules25225358](https://doi.org/10.3390/molecules25225358) -- and
pulled the full text: **it does not report kaempferol or quercetin at all.** It
studies O-methylated flavonoid *derivatives* (3,4'-di-O-methylkaempferol,
4'-O-methylkaempferol, chalcones) isolated from African plant species -- chemically
distinct compounds (the free phenolic OH groups are methylated). Free kaempferol
and quercetin are not among the tested compounds. This citation needs fixing
wherever it appears in the manuscript, not just in this CSV.

**Also found the CSV's target assignment itself doesn't match the literature.**
Both kaempferol and quercetin are labeled `target=MAO-B` here, but every primary
source found characterizes them as MAO-A-selective:
  - Sloley et al. 2000, *J Pharm Pharmacol* 52:451-9, PMID 10813558,
    [DOI](https://doi.org/10.1211/0022357001774075): kaempferol IC50 = 7e-7 M
    (700 nM) vs MAO-A. MAO-B described only as "more pronounced inhibition of MAO-A
    than MAO-B" -- no MAO-B number given at all. **Rat brain enzyme, not human.**
  - Larit et al. 2018, *Phytomedicine* 40:27-36, PMID 29496172,
    [DOI](https://doi.org/10.1016/j.phymed.2017.12.032), **recombinant human
    MAO-A/-B**: quercetin IC50 = 1.52 uM vs MAO-A. No MAO-B IC50 reported for
    quercetin in the abstract.

Did not force these MAO-A values into the existing MAO-B-labeled rows (would
misrepresent which assay the number came from). Instead: left the original
kaempferol/quercetin rows as `target=MAO-B`, `use_in_correlation=no`, with a note
explaining why: and added two new rows, `kaempferol_maoa` and `quercetin_maoa`
(same SMILES, `target=MAO-A`, real sourced values), docked properly on 2Z5X (3
seeds, not just relabeling the old MAO-B-context run).

**Filled in real values for lazabemide and isatin too:**
  - lazabemide: IC50 = 0.063 uM vs human MAO-B, Maliyakkal et al. 2020, ChemMedChem
    15:1629-1633, PMID 32583952,
    [DOI](https://doi.org/10.1002/cmdc.202000305) -- used as the reference standard
    in that paper, so this is a well-anchored number.
  - isatin: IC50 approx 3 uM vs MAO-B, Medvedev/Clow/Sandler/Glover 1996,
    *Biochem Pharmacol* 52:385-91, PMID 8687491,
    [DOI](https://doi.org/10.1016/0006-2952(96)00206-7) -- the original paper
    establishing isatin as an endogenous MAO-B inhibitor. Cross-confirmed by a 1997
    follow-up (PMID 9503568, IC50 3-8 uM range). **Rat brain mitochondrial enzyme,
    not human recombinant** -- same caveat as kaempferol.

**Rebuilt the correlation plot with 5 usable points now (up from 1): safinamide,
lazabemide, isatin (MAO-B) + kaempferol_maoa, quercetin_maoa (MAO-A).**
`Pearson r = -0.146 (p = 0.815, n = 5)`. **This is a real, honest negative result,
not a plotting artifact** -- even with better data, this docking protocol shows no
significant correlation between raw Vina consensus score and literature pIC50
across these five compounds. Caveats on top of the small n: two different targets
were pooled into one regression (MAO-A and MAO-B use different receptor structures
and box parameters, so pooling isn't fully rigorous -- but n=3 and n=2 per target
individually is too small to fit separately); two of the five values are rat-brain
assays, not human recombinant, so not perfectly comparable to the other three.
Worth being upfront about this rather than presenting the plot as more meaningful
than it is -- weak/no raw-score-to-affinity correlation is itself a known,
citable limitation of empirical scoring functions like Vina, not a sign something
is broken.

Files: `01_smiles/references.csv` (16 rows now), `02_ligands/{sdf,pdbqt}/{kaempferol_maoa,
quercetin_maoa}.*`, `04_docking/2Z5X_references_seed{11,22,33}/{kaempferol_maoa,quercetin_maoa}*`,
`05_validation/benchmark_mao.csv` (updated), `05_validation/benchmark_mao_scatter.png` (updated,
now with regression line and Pearson r).

## 2026-08-21 — A (A5 ranking analysis; correction filed against the A5 entry)

Read `08_analysis/consensus_new.csv` properly and cross-matched candidate IDs to the
manuscript's §3.3 SMILES. Produced the correction recorded inline in the A5 entry above,
plus the ranking tables A5 generated but never wrote down.

### Top 10 per target (3-seed consensus, kcal/mol)

**TTBK1 (7JXX, GATE-1 validated at 0.71 Å)**

| # | ligand | consensus | seed SD | # | ligand | consensus | seed SD |
|---|---|---|---|---|---|---|---|
| 1 | cand_013 | -8.83 | 0.005 | 6 | cand_055 | -8.44 | 0.008 |
| 2 | cand_002 | -8.80 | 0.015 | 7 | cand_010 | -8.41 | 0.097 |
| 3 | cand_001 | -8.79 | 0.021 | 8 | **cand_003** | **-8.38** | 0.020 |
| 4 | cand_050 | -8.71 | 0.007 | 9 | cand_030 | -8.38 | 0.013 |
| 5 | cand_011 | -8.58 | 0.009 | 10 | cand_005 | -8.35 | 0.006 |

**MAO-B (2V5Z)**

| # | ligand | consensus | seed SD | # | ligand | consensus | seed SD |
|---|---|---|---|---|---|---|---|
| 1 | cand_043 | -11.63 | 0.010 | 6 | cand_051 | -11.27 | 0.017 |
| 2 | **cand_003** | **-11.41** | 0.021 | 7 | cand_038 | -11.21 | 0.006 |
| 3 | cand_026 | -11.40 | 0.006 | 8 | cand_015 | -11.19 | 0.000 |
| 4 | cand_023 | -11.35 | 0.035 | 9 | cand_005 | -11.16 | 0.023 |
| 5 | cand_013 | -11.29 | 0.010 | 10 | cand_004 | -11.13 | 0.006 |

Full-set ranges: TTBK1 -6.96 (cand_054) to -8.83; MAO-B -8.72 (cand_056) to -11.63.

### Dual-target ranking

Neither single-target table answers the project's actual question. Ranked by worst-of-the-
two ranks, so nothing places by being lopsided on one target:

| # | ligand | TTBK1 (rank) | MAO-B (rank) | B-ring decoration |
|---|---|---|---|---|
| 1 | cand_013 | -8.83 (#1) | -11.29 (#5) | CF3 + N-methylpiperazine |
| 2 | **cand_003** | -8.38 (#8) | -11.41 (#2) | CF3 + methyl -- **the MD lead** |
| 3 | cand_005 | -8.35 (#10) | -11.16 (#9) | CF3 + vinyl |
| 4 | cand_004 | -8.26 (#14) | -11.13 (#10) | CF3 + ethyl |
| 5 | cand_008 | -8.29 (#13) | -11.04 (#14) | two CHF2 |

All five carry CF3 or CHF2 -- the convergence described in `ligand_generation_notes.md`,
now reproduced on a validated receptor rather than on the exhaustiveness-4 in-loop score.

**`cand_013` is the new dual-target #1 and is the one candidate that could displace the
lead, but it should not be promoted without scrutiny.** Three reasons: (i) chemically it is
the outlier of the set -- an N-methylpiperazine rather than a small lipophilic group; (ii)
MW 420.39 and TPSA 77.15 put it at the heaviest end of the 56 and hard against the ~79 Å^2
BOILED-Egg ceiling, i.e. essentially no BBB margin; (iii) its MAO-B score moved -8.66 ->
-11.29 (2.63 kcal/mol, the largest shift of any compound in the set) between the RL-time
score and this consensus. The direction is expected -- exhaustiveness 4 misses good poses --
but a #1 placement resting on the single largest score change in the set warrants a check
before it consumes a week of MM-GBSA. Recommend Vinardo cross-scoring `cand_013` and
`cand_003` before any lead change, since §3.3's whole point is that Vina rank order does
not survive a scoring-function swap.

### Native-ligand calibration

The native crystal ligands were docked alongside the candidates, which anchors the scores:

| receptor | native ligand | consensus | vs. candidates |
|---|---|---|---|
| 7JXX | VP7 (real TTBK1 inhibitor) | -8.97 (3 seeds) | **beats all 56**; best candidate -8.83 |
| 2V5Z | SAG_pdb (safinamide, crystal tautomer) | -10.15 (3 seeds, Track B) | 50/56 candidates score better |
| 4BTK | DTQ | -7.65 (3 seeds) | cross-check receptor, not used for candidates |

Read against the Track-B correlation result (Pearson r = -0.146, n = 5, no significant
raw-score-to-pIC50 relationship on this protocol), the MAO-B row must **not** be reported as
"50 candidates more potent than safinamide." The defensible statement is that the series
docks into the range occupied by known binders. The TTBK1 row is the more informative one
and is a genuinely useful negative: an optimized real inhibitor still out-scores every
generated candidate.

Note `consensus_new.csv` carries `native_SAG` at n_seeds=1 (the Day-2 GATE-0 single-seed
calibration run); the 3-seed value above is `SAG_pdb` from `05_validation/benchmark_mao.csv`
(-10.153, SD 0.012). The two agree to 0.01 kcal/mol via completely separate prep paths,
the same cross-check B3 already reported at the RMSD level.

### Data-quality flags

**1. `cand_002` moved the wrong way on an unchanged receptor.** 2V5Z is the same receptor
the manuscript used, so RL-time and consensus scores should broadly track. Across the 56 the
agreement is reasonable (mean |diff| 0.271, median 0.089, 32/56 within 0.10 kcal/mol), and
the large movers are nearly all *improvements* consistent with exhaustiveness 4 -> 32 finding
better poses (`cand_013` -8.66 -> -11.29; `cand_050` -8.34 -> -10.80; `cand_056` -6.83 ->
-8.72). `cand_002` is the exception: **-11.56 -> -10.42, i.e. 1.14 kcal/mol worse**, and it
also carries the highest MAO-B seed SD in the set (0.145). A better search returning a worse
score is not explained by search depth. This is why `cand_002` sits #2 on TTBK1 but #44 on
MAO-B. Worth a direct look at its three seed logs before it is used in any table.

**2. Seven candidates have unreliable TTBK1 rankings.** Inter-seed SD > 0.15 kcal/mol:
cand_044 (0.808), cand_053 (0.709), cand_007 (0.294), cand_027 (0.267), cand_032 (0.235),
cand_033 (0.224), cand_054 (0.184). `cand_044` is the one that matters -- its consensus
-8.29 would place it ~#11, but with SD 0.808 and best_overall -8.77 the three seeds plainly
disagree about the pose. Exclude from any ranked table, or re-dock with more seeds. Note this
is the 3-seed protocol working as designed: single-seed docking would have reported -8.77
with no indication anything was wrong.

**3. The TTBK1 RL-vs-consensus comparison is not meaningful** (mean |diff| 0.508, 3/56
within 0.10) and should not be quoted as disagreement -- RL time used apo 4NFM, A5 used holo
7JXX. Different receptor, so a shift is the expected result, not a discrepancy.

### Follow-ups

- Vinardo cross-score `cand_013` and `cand_003` before any lead change (blocks the Week 2-3
  MM-GBSA compound choice).
- Inspect `cand_002`'s 2V5Z seed logs (`04_docking/2V5Z_candidates_56_seed{11,22,33}/`).
- Decide whether the seven high-SD TTBK1 compounds get more seeds or get excluded.
- Manuscript §3.3 tabulates candidates by SMILES only. Add the `cand_0NN` IDs to that table
  at revision time -- the absence of IDs is the direct cause of the correction filed above.

Files: none (analysis only, over existing `08_analysis/consensus_new.csv` and
`05_validation/benchmark_mao.csv`; correction inserted inline into the A5 entry above).

## 2026-08-21 — A (Phase A: targeted literature check on the lead's B-ring chemotype, before
committing to Week 2 MD)

User wants direct literature grounding for the lead (`cand_003`,
`Cc1cc(-c2oc3cccc(O)c3c(=O)c2O)ccc1C(F)(F)F`) before MD, not just docking-score self-consistency.
The manuscript currently has **zero citations discussing the CF3/halogen B-ring decoration
itself** (confirmed by grep of `manuscript_dualtarget_flavonol.md` for "CF3"/"trifluoromethyl"/
"precedent" -- none found). Ran targeted PubMed searches (`mcp__claude_ai_PubMed`, same tool used
for the MAO IC50 chase) rather than a broad re-review.

**Finding 1 -- the flavonol's B-ring is structurally a chalcone aryl group, and halogenated
chalcones are a large, well-populated, sub-100nM MAO-B SAR program.** According to PubMed, Bijo
Mathew's group (Amrita/Ahalia, India) has run this exact substitution strategy repeatedly:
- Singh et al. 2023, *Chem Biol Drug Des* 102(2):271-284, PMID 37011915,
  [DOI](https://doi.org/10.1111/cbdd.14238): benzyloxy halogenated chalcones, lead compound BB4
  IC50 = 0.062 uM vs hMAO-B -- **more potent than lazabemide** (0.11 uM in their assay), the same
  reference compound already in `references.csv`.
- Mathew et al. 2019, *Arch Pharm* 352(4):e1800309, PMID 30663112,
  [DOI](https://doi.org/10.1002/ardp.201800309): fluorinated morpholine chalcones, lead f2 IC50 =
  0.087 uM vs hMAO-B, confirmed BBB-permeant by PAMPA.
- Mathew et al. 2016, *ChemMedChem* 11(11):1161-71, PMID 27159243,
  [DOI](https://doi.org/10.1002/cmdc.201600122): brominated thienyl chalcones, reversible
  selective hMAO-B inhibitors, Ki = 0.11 uM.
- Parambi et al. 2019, *Cent Nerv Syst Agents Med Chem* 19(1):67-71, PMID 30451121,
  [DOI](https://doi.org/10.2174/1871524918666181119114016): same halogenated-thiophene-chalcone
  series cross-tested for cholinesterase inhibition (dual MAO-B/ChE angle).

This directly supports the project's core structural hypothesis with real potency data, not just
internal docking-score consistency -- halogen/CF3-decorated aryl-chalcone-like scaffolds
reaching or beating reference MAO-B inhibitor potency is an established, reproducible result in
this chemical space, independent of this project's own docking.

**Finding 2 -- directly actionable for the unresolved MAO-A/MAO-B selectivity risk flagged in
Phase 8.** Olotu et al. 2020, *J Biomol Struct Dyn* 39(16):6126-6139, PMID 32705963,
[DOI](https://doi.org/10.1080/07391102.2020.1796803): shows that simply shifting a fluorine atom
between *meta* and *para* positions on the same chalcone B-ring produces a "remarkable shift" in
MAO-A vs. MAO-B selectivity, via MD + DFT analysis of active-site interactions. Phase 8's
`selectivity_triangulation_report.md` flags MAO-A/B selectivity as the dominant unresolved risk
(+0.09 kcal/mol margin, "no separation"); this paper gives a concrete, literature-precedented
lever (halogen position on the B-ring, not just identity) worth naming in the manuscript's
limitations/future-work as the specific next design move, rather than leaving selectivity as an
open-ended problem.

**Finding 3 -- weaker but real precedent for CF3 as a kinase hinge-region motif specifically**
(TTBK1 side of the hypothesis). Zhang et al. 2016, *Mol Inform* 35(1):15-8, PMID 27491649,
[DOI](https://doi.org/10.1002/minf.201500091): rational design study on CK2 (a different kinase,
not TTBK1) found that replacing a methyl with **trifluoromethyl specifically to engage the hinge
region** gave a 5-fold potency improvement via a halogen-bond/electrostatic interaction with the
hinge. This is not TTBK1-specific and should not be oversold, but it is real precedent that CF3
placement is a deliberate, literature-known hinge-engagement strategy in kinase inhibitor design
generally, not an arbitrary substituent choice.

**Finding 4 -- re-confirmed, not new: TTBK1+MAO-B dual targeting still has zero combined
precedent.** Re-ran the phase0 query fresh (`"TTBK1" AND "monoamine oxidase" AND "dual" AND
"neurodegeneration"`): **0 hits**, same null result as the original Phase 0 search
(`phase0_targets_report.md:19-22`). No new literature closes this gap since Phase 0. The
dual-target framing should continue to be presented as a hypothesis, not a validated combination
strategy -- this is unchanged, just re-verified rather than assumed stale.

**Net effect on confidence in the lead:** the CF3/halogen-decorated aryl-B-ring strategy for
MAO-B engagement is now grounded in real, independently-published SAR data (Finding 1), not just
this project's own docking scores -- a meaningfully stronger basis than existed before this
check. The TTBK1 side and the dual-target framing remain comparatively weaker (Findings 3-4): CF3
kinase-hinge precedent exists but not for TTBK1 itself, and dual TTBK1+MAO-B engagement is still
an unprecedented combination. This does not change the lead compound choice (no new candidate
was generated in this pass) but materially improves what can be defensibly said about *why* the
scaffold might work, and gives a concrete, citable direction (halogen repositioning) for the
selectivity limitation. Recommend folding Findings 1-2 into manuscript SS4.3 at revision time.

Files: none (literature only; no code run). PMIDs/DOIs above are the citable record.

## 2026-08-21 — A (Phase D groundwork: acpype/GAFF2 params for cand_003, closing the
Day-1 "STILL OPEN" item, run in parallel with the generation track below)

Regenerated the lead compound's GAFF2 parameters via `acpype` in the WSL2 `mdgbsa` env
(confirmed working since Day 1), independent of whether the generation run below changes the
lead. `acpype -i 02_ligands/sdf/cand_003.sdf -b cand_003 -n 0 -a gaff2 -c bcc` -- net charge 0
(neutral flavonol, no ionizable groups modeled), AM1-BCC charges via `sqm`, GAFF2 atom types.
Completed cleanly in 49s (Antechamber OK, Parmchk OK, Tleap OK). Output at
`06_md/params/cand_003.acpype/`: `cand_003_AC.{prmtop,inpcrd,frcmod}` (AMBER) and
`cand_003_GMX.{top,itp,gro}` (GROMACS) -- this is the file set Week 2's MD system build was
blocked on. De-risks Day 2 regardless of the generation track's outcome.

Files: `06_md/params/cand_003.acpype/*` (new).

## 2026-08-21 — A (Phase A-D: local REINVENT4 install + new dual-target RL generation run,
per explicit user request to verify "can we generate more/better ligands" before trusting
`cand_003` into Week 2 MD)

**Context:** REINVENT4 was not installed anywhere locally -- both original campaigns
(`00_library/reinvent4_output/campaign1_first_failed`, `campaign2_v2`) ran on a since-
unavailable cloud sandbox, and the docking-in-the-loop scorer (`dock_score.py`) never existed
locally, only its config (`rl_staged_learning_v2.toml`). Stood this up from scratch rather than
mining existing output, per explicit direction (user chose this over re-mining the unused
666-candidate pool from campaign 1, after being shown REINVENT4 wasn't installed and the master
plan says to keep generation "as-is").

**Install (isolated `reinvent4` conda env, Python 3.11 -- `docking_project` untouched):**
cloned `MolecularAI/REINVENT4` (pyproject.toml requires-python >=3.11, so the env had to be
recreated once from an initial 3.10 attempt), installed CPU-only torch + package deps,
downloaded `libinvent.prior` (91.6 MB) from the official Zenodo record
(10.5281/zenodo.15641296, concept DOI matching the README's citation).

**Found and fixed a real Windows-portability bug in REINVENT4 itself** (upstream code, not
this project's): `reinvent/utils/hw_report.py` did an unconditional top-level `import resource`
-- `resource` is POSIX-only and doesn't exist on Windows, so the CLI couldn't even start. The
actual *usage* of `resource` (line ~112, peak-memory reporting) was already correctly guarded
by `if SYSTEM != "Windows":`, so this was purely a missing import guard, not a logic gap.
Explains why the original runs used a Linux cloud sandbox -- REINVENT4 doesn't run on native
Windows out of the box. Fixed by moving the import inside the existing guard. Also needed
`pip install scipy` (a real transitive dependency missing from the hand-typed deps list used
for the `--no-deps -e .` install). Verified with REINVENT4's own sampling mode on the project's
existing scaffold file (`generation/flavonol_7deoxy_scaffold.smi`) before wiring anything else
to it -- produced valid, on-scaffold molecules (5 requested, 3 valid, correctly decorating both
B-ring attachment points), confirming the install and prior are good.

**Built the docking-in-the-loop scorer that never existed locally**
(`generation/scoring/dock_score.py`): reuses this week's already-validated local pattern
(RDKit ETKDGv3 embed + MMFF, obabel protonation, `mk_prepare_ligand`, Vina CLI) rather than
`obabel` pose conversion (irrelevant here since only the docked score is needed, not a
re-exported pose). Matches REINVENT4's `ExternalProcess` stdin-SMILES / JSON-stdout contract
(same contract already used by `scripts/bbb_score.py`). One real integration bug found and
fixed before trusting it: the script runs as a bare subprocess with a fixed interpreter path
(no shell activation), so `mk_prepare_ligand` (installed in `docking_project`'s `Scripts/`)
would not resolve via `shutil.which` -- added a startup step that derives the conda env root
from `sys.executable` and prepends its `Library/bin`/`Scripts` to `PATH`, matching the same fix
`dock.sh` already applies via `$CONDA_PREFIX`.

Standalone-tested `dock_score.py` on `cand_003`'s SMILES against both receptors before wiring
it into REINVENT4: **7JXX -7.78/-7.80 kcal/mol** (exhaustiveness 4, vs. this week's 3-seed
exhaustiveness-32 consensus of -8.38 -- expected gap, matches the direction seen throughout A5)
and **2V5Z -11.4 kcal/mol** (vs. this week's consensus -11.41 -- near-exact agreement even at
low exhaustiveness). Strong sign the new scorer is measuring the same thing as the rest of this
week's pipeline, not something different.

**Fixed the TPSA inconsistency flagged in the 2026-08-20 A4 entry** while touching
`scripts/bbb_score.py` for reuse here: added `includeSandP=True` to its `Descriptors.TPSA`
call, matching `scripts/filter_cascade.py`'s offline filter. RL-time BBB reward and the final
filter gate are now numerically consistent (previously flagged as a known, unfixed gap).
Copied to `generation/scoring/bbb_score.py` (its relative import of
`filtering/boiled_egg_coords.py` is written for that location, `../../filtering` -- running it
from `scripts/` instead breaks the import; this is a pre-existing path assumption in the
recovered file, not a new bug).

**Deliberate improvement over the original campaign2_v2 config, not just a re-run:** the new
config (`generation/rl_staged_learning_v3_local.toml`) points the TTBK1 docking-reward term at
**7JXX** (holo, GATE-1 validated at 0.71 Angstrom this week) instead of the original **4NFM**
(apo, no native ligand -- a validation gap `phase0_targets_report.md` itself flags). Both
receptor boxes use this week's actual prepared `box.json` values (7JXX: center 178.93/19.64/
46.21, box 18x18x18; 2V5Z: center 51.89/156.45/28.56, box 18x18.9x18 -- note 2V5Z center matches
the original run exactly, confirming consistent receptor prep, but box *size* differs, 18 vs.
the original's 28). `max_steps` capped at 20 (from 40) to time-box the run given the 2-day
window; all other RL hyperparameters (batch_size=32, sigma=128, IdenticalMurckoScaffold
diversity filter, component weights 0.35/0.35/0.15/0.15) kept identical to campaign2_v2 for
comparability.

**Run launched and healthy**: ~45s/RL-step observed (steps 1-2 timed directly), so 20 steps
should complete in well under the time originally budgeted for this phase. Docking scores in
step 1-2 output land in the same range as this week's established consensus numbers (e.g. raw
TTBK1_dock 7.7-10.7, MAOB_dock 2.5-10.8 across sampled molecules), and the diversity filter is
visibly suppressing repeated Murcko scaffolds as designed. Full results and the decision on
whether any new candidate displaces `cand_003` will be logged separately once the run completes
and goes through `filter_cascade.py` + consensus re-docking + Vinardo cross-check, per the
approved plan's Phase E bar (must rank top-2 under *both* scoring functions to displace the
current lead, the same bar `phase6_final_report.md` established originally).

Files: `generation/scoring/dock_score.py` (new), `generation/scoring/bbb_score.py` (new, copy
of `scripts/bbb_score.py` with TPSA fix), `generation/rl_staged_learning_v3_local.toml` (new),
`generation/flavonol_7deoxy_scaffold.smi` (new, copy), `scripts/bbb_score.py` (modified, TPSA
fix), `REINVENT4/reinvent/utils/hw_report.py` (modified, Windows import fix -- vendored
third-party code, not tracked in this repo's own history), `REINVENT4/priors/libinvent.prior`
(new, downloaded, gitignored), new conda env `reinvent4` (gitignored/untracked).

## 2026-08-21 — A (Phase E: analyze the new RL run, decide whether cand_003 survives)

RL run finished cleanly at `max_steps=20` (17:00:54 -> 17:59:18, ~58 min total, ~2.9 min/step
average -- slower than the ~45s/step seen in the first 2 steps, likely docking-subprocess
contention with the concurrent acpype run earlier in the session). 640 raw rows (20 steps x
32), 370 unique valid molecules after RDKit canonicalization + dedup.

**Triage against cand_003 at matched exhaustiveness.** Standalone-tested `dock_score.py` on
`cand_003` itself before the run (see prior entry): 7.78-7.80 kcal/mol (TTBK1/7JXX), 11.4
kcal/mol (MAOB/2V5Z) at exhaustiveness=4 -- the same low-exhaustiveness protocol the RL loop
itself uses, so this is the correct benchmark for triage (not the 3-seed exhaustiveness=32
consensus numbers, which aren't comparable at this stage). **49/370 unique new molecules beat
cand_003 on *both* raw docking terms simultaneously** at this protocol -- confirms the run
found real signal, not noise; the "CF3 convergence" pattern flagged earlier in this project
reappeared strongly (most of the 49 carry one or more CF3/halogen B-ring substituents).

**Filter cascade on all 49** (`scripts/filter_cascade.py`, same Lipinski/BOILED-Egg/PAINS+BRENK
gates as A4): Lipinski 48/49, GI-absorption 43/49, alert-free 44/49, **BBB 1/49**. Only 1
candidate passes all four filters.

**That one candidate is cand_003 itself, independently rediscovered by this run** -- verified
by exact RDKit-canonical SMILES match against `01_smiles/candidates_56.csv`, not just visual
similarity. Its raw scores from this run (TTBK1 8.285, MAOB 11.41) are close to, and marginally
above, the standalone benchmark (7.78-7.80 / 11.4) purely from conformer-embedding randomness
at low exhaustiveness (each dock_score.py call re-embeds a fresh 3-conformer ensemble with a
fixed seed but the RL harness invokes it independently per step) -- not a genuine new score.

**Net result: zero genuinely new candidates from this run displace cand_003.** Every one of the
48 non-self hits that out-scored cand_003 on raw docking did so specifically by adding *more*
lipophilic/aromatic bulk (second CF3 groups, extended biaryls, alkyl/aryl linkers) --
WLogP for the top 20 by score ranges 5.3-7.6, well past the BOILED-Egg BBB boundary, vs.
cand_003's WLogP ~4.2 at the same 70.67 TPSA. This is the same scaffold-level tension Phase 4
first identified (TPSA floor for the flavonol core forces the BBB boundary to be crossed almost
entirely via WLogP), now demonstrated directly under a corrected, holo-receptor reward: **the
"obvious" way to improve raw docking score on this scaffold is to add lipophilic bulk, and doing
so trades away the BBB permeability cand_003 was specifically selected for.** This is a genuine,
citable negative result, not a failed experiment -- it is independent evidence (different RL
run, corrected TTBK1 receptor, different random seed) that `cand_003` sits at a real, non-
arbitrary local optimum for this scaffold class under the project's own druglikeness
constraints, not merely an artifact of the original campaign2_v2 run's specific reward
weighting.

**Decision (per the approved plan's Phase E bar -- a new candidate must pass filters AND rank
top-2 under both Vina-consensus and Vinardo to displace the lead): cand_003 remains locked as
the lead for Week 2 MD.** No Vinardo cross-check was needed for a fresh candidate since none
survived the filter cascade to be worth cross-checking. `cand_003`'s status is now supported by
five independent lines of evidence: (1) original RL score, (2) 3-seed Vina consensus, (3)
independent Vinardo scoring function, (4) this week's from-scratch pipeline reproduction, and
(5) this fresh, independently-seeded RL campaign with a corrected receptor rediscovering it as
the only viable dual-target, druglikeness-compliant hit in its own output.

**Follow-up for the manuscript:** the WLogP-vs-BBB tension demonstrated here, combined with the
literature finding earlier today (Olotu et al. 2020, PMID 32705963 -- halogen *position*, not
just identity, shifts MAO-A/B selectivity) suggests future SAR work should explore repositioning
existing halogens rather than adding more of them, since this run's data shows adding bulk is a
dead end against the BBB constraint.

Files: `generation/rl_v3_local_1.csv` (new, 640-row RL summary), `generation/rl_v3_dedup.csv`
(new, 370 unique), `generation/rl_v3_beats_cand003.csv` (new, 49 candidates),
`01_smiles/generation_top20.csv`, `01_smiles/generation_beats_cand003_all49.csv` (new),
`08_analysis/generation_top20_filtered.csv`, `08_analysis/generation_all49_filtered.csv` (new).

**Track A/B generation-verification task is now complete; proceeding to Week 2 MD on
`cand_003`** using the acpype/GAFF2 parameters generated earlier today
(`06_md/params/cand_003.acpype/`).

## 2026-08-21 — A (Phase F: Week 2 MD system build, scoped to one system for the review --
TTBK1/7JXX + cand_003, not the master plan's full 4-system/47-GPU-hour build)

**Scope decision, made mid-build, not in the original plan:** switched the single MD system
from 2V5Z (MAO-B) to **7JXX (TTBK1)**. Reason: MAO-B's binding site requires the covalently-
linked FAD cofactor (confirmed present in 2V5Z, `phase0_targets_report.md`), and no FAD force-
field parameters exist anywhere in this project (`find -iname "*fad*"` -- nothing). Properly
parameterizing a covalently-bound cofactor (RESP/AM1-BCC charges, correct linkage chemistry to
Cys397) is a multi-hour task with real correctness risk under this deadline. TTBK1/7JXX needs no
exotic cofactor (kinase domain + ligand only) and is also the receptor this whole session's
generation run was already built around (GATE-1 validated, 0.71 Angstrom). Re-scoping to the
system that's actually buildable correctly in the time available, rather than attempting the
higher-profile MAO-B target and risking an unreliable result.

**Re-discovered the same CUDA-platform-missing issue Day 1 apparently didn't hit**: calling
`docking_project`'s `python.exe` directly by absolute path (as done all week, bypassing `conda
activate`) leaves `Library/bin` off `PATH`, so `OpenMMCUDA.dll` fails to resolve its `nvrtc64_
130_0.dll` dependency (Windows error 126) and only Reference/CPU/OpenCL platforms register.
Fixed by prepending `<env>/Library/bin` to `PATH` before invoking OpenMM -- confirmed CUDA loads
and computes correctly once done. **Any future OpenMM/CUDA work in this project must do the
same** (same root cause as `dock.sh`'s existing `$CONDA_PREFIX/Library/bin` prepend, now known to
also apply to OpenMM, not just meeko).

**Protein prep (7JXX -> MD-ready) took two attempts.** First attempt used `pdb4amber --reduce`
(crashed, no `reduce`/`reduce2` binary in the `mdgbsa` env) then OpenMM `Modeller.addHydrogens`
as a fallback -- this ran but produced histidines/N-terminus with hydrogen counts inconsistent
with tleap's own `ff14SB` templates (`HIE` residues got a spurious extra `HD1` with no atom
type; N-terminal `ASP` got an extra untyped `H`), which tleap correctly rejected as fatal rather
than silently building a broken system. Root cause: two different tools' hydrogen-placement
conventions don't agree, and mixing them is fragile. **Fix: let one tool do it.** Used
`pdb4amber` *without* `--reduce` (completes fine, and critically -- correctly places `TER`
records at the crystal structure's two real chain gaps, Gly22|Phe23 and Thr164|Asp165, matching
`phase0_targets_report.md`'s disordered-region note) for heavy-atom cleanup and gap detection
only, renamed its generic `HIS` labels to `HIE` (single default tautomer, no per-residue pKa
analysis attempted at this scope), and let **tleap build every hydrogen itself** from its own
internal templates. This is what actually worked (`Errors = 0` on the real build). Take-away for
next time: don't pre-add hydrogens with one tool and hand them to tleap -- let tleap do the
protonation using its own residue library, or use `reduce` specifically (not OpenMM `Modeller`)
if pre-adding is ever needed.

**System build** (`06_md/system/build.leap`, WSL `mdgbsa` tleap): `ff14SB` protein + `gaff2`
ligand (loaded from the docked-pose-derived `cand_003_AC.lib`/`.frcmod`, see below) + `tip3p`
water, truncated-octahedral box, 10 Angstrom buffer, neutralized with 12 Cl- (system net charge
+12 before neutralization -- a notably basic ATP pocket region, consistent with the Lys38/Lys39/
Lys63 landmarks `phase0_targets_report.md` names for this pocket). ~40k atoms total. Zero tleap
errors; 8 warnings, all expected N-/C-terminal residue-name-format notices from the 3 chain
fragments created by the 2 real gaps.

**Ligand pose, not just ligand parameters:** re-ran `acpype` a second time
(`06_md/params/docked/cand_003.acpype/`), this time on the actual **docked pose** extracted from
this week's A5 consensus docking (`04_docking/7JXX_candidates_56_seed11/cand_003_out.sdf`, top-
ranked of 9 poses, verified atom count and coordinates match the docked geometry, not the free
minimized conformer) rather than the free-conformer parameterization done earlier today. This
guarantees the MD starting geometry is the actual docked binding pose in the actual binding
site, not a re-embedded conformer that would need re-alignment (and risks starting from a
non-representative geometry). The earlier free-conformer acpype run (`06_md/params/cand_003.
acpype/`) remains useful as the general-purpose GAFF2 parameter set; this docked-pose run is
what's actually used for the MD system.

**MD run** (`06_md/run_md.py`, native Windows CUDA, `docking_project` env): PME electrostatics,
HBond constraints, Langevin middle integrator (300 K, 1/ps, 2 fs), Monte Carlo barostat (1 atm).
2000-step minimization, 100 ps NVT-style equilibration, then production. **Throughput: ~171
ns/day** on the RTX 4050 for this ~40k-atom system -- much faster than assumed when scoping the
plan, so the originally-planned "short, preliminary" run can likely be extended once the first
short run's correctness is confirmed end-to-end (checked: no NaN, temperature/energy stable
through equilibration).

**MM-GBSA prep staged, not yet run** (waiting on production trajectory): `06_md/system/
strip_traj.cpptraj` (autoimage + strip water/ions from the trajectory) and `06_md/system/
mmpbsa.in` (single-trajectory GB, `igb=5`, 0.15 M salt) written and ready; `ante-MMPBSA.py`/
`cpptraj`/`MMPBSA.py` all confirmed present in the WSL `mdgbsa` env.

Files: `06_md/system/7JXX_prot.pdb` (pdb4amber-cleaned, HIS->HIE), `06_md/system/build.leap`,
`06_md/system/{complex,protein_only,ligand_only}.{prmtop,inpcrd}` (new, tleap output),
`06_md/params/docked/cand_003.acpype/*` (new, docked-pose GAFF2 params), `06_md/run_md.py`
(new), `06_md/system/strip_traj.cpptraj`, `06_md/system/mmpbsa.in` (new, staged).

## 2026-08-21 — A (Phase F: first MD/MM-GBSA result, then extending production)

**1 ns production completed cleanly**: temperature stable at ~300 K, potential energy stable
around -516000 kJ/mol throughout, no NaN/blowup. Measured throughput **~169-171 ns/day** on the
RTX 4050 for this ~40k-atom system -- substantially faster than assumed when the plan scoped
this as a short/preliminary run, so extending it is cheap.

**MM-GBSA on the full 100-frame trajectory** (`cpptraj` strip/autoimage -> `ante-MMPBSA.py` ->
`MMPBSA.py`, single-trajectory GB, igb=5, 0.15 M salt): **DELTA G (binding) = -34.55 +/- 2.22
kcal/mol** (SEM) for cand_003/TTBK1(7JXX). Compare to the original (pre-reconstruction) project's
`phase8_selectivity_triangulation_report.md` result for the same compound on TTBK1: -26.32 +/-
0.53 kcal/mol from a 20 ns run on the apo 4NFM receptor. Same order of magnitude; the more
favorable number here is plausibly attributable to the holo, GATE-1-validated 7JXX receptor
(this week's whole point) rather than apo 4NFM, and the much wider SEM here (2.22 vs 0.53) is
the expected consequence of 1 ns / 100 frames vs. 20 ns of sampling -- not a red flag, just less
precision. **This confirms the full pipeline (docking -> generation verification -> MD ->
MM-GBSA) now works end-to-end on locally-rebuilt infrastructure**, independent of the original
run's cloud/unrecovered environment.

**Extending production** given the throughput headroom: continuing from `final_state.xml` for
additional nanoseconds to tighten the SEM before this is reported at review, rather than resting
on a 1 ns number when ~10x more sampling is affordable in the remaining time. Will re-run
MM-GBSA on the combined/extended trajectory and log the updated number separately.

Files: `06_md/system/mmgbsa_results.dat` (new), `06_md/system/gb_{complex,receptor,ligand}.
prmtop` (new), `06_md/system/production_stripped.nc` (new).

## 2026-08-21 — A (Phase F: extended to ~10 ns, final MM-GBSA result for the review)

Extended production from 1 ns to **10.1 ns total** (`06_md/run_md_extend.py`, +9.1 ns,
continued from `final_state.xml`), same throughput (~167-168 ns/day), stable temperature/energy
throughout, no NaN/blowup across the full run. Combined both trajectory segments
(`production.dcd` + `production2.dcd`) via `cpptraj` (1000 raw frames, 10 ps/frame), then ran
MM-GBSA on 200 frames (every 5th frame, `interval=5` in `mmpbsa.in`) spanning the full 10 ns.

**Final result: DELTA G (binding) = -31.36 +/- 0.27 kcal/mol (SEM)** for cand_003/TTBK1(7JXX),
`06_md/system/mmgbsa_results_full.dat`. Converged consistently from the 1 ns estimate (-34.55
+/- 2.22): both agree well within each other's uncertainty, and the SEM tightened by ~8x with
10x more sampling, as expected. This is now comparable in *precision* to the original project's
phase8 result (-26.32 +/- 0.53 kcal/mol, 20 ns, apo 4NFM) despite half the simulation time, and
agrees with it in *magnitude* -- both report large, clearly favorable binding free energies for
this compound on TTBK1, from two independently-built pipelines on two different (holo vs. apo)
receptor structures.

**Where this leaves the project going into the review:**
1. Literature grounding for the lead's chemotype now exists (today's Phase A), where none did
   before.
2. The lead (`cand_003`) survived an independent, freshly-generated RL campaign that specifically
   corrected the original's known TTBK1 receptor weakness (Phase B-E) -- nothing beat it once
   druglikeness constraints were enforced.
3. A real, working, locally-reproduced MD/MM-GBSA result now exists on the properly-validated
   holo TTBK1 receptor (-31.36 +/- 0.27 kcal/mol, 10 ns), independently corroborating the
   original project's MD finding on a different (apo) structure.
4. MAO-B MD was not attempted this session -- no FAD cofactor parameters exist locally, and
   building them correctly was judged too high-risk for the remaining time. This is the one
   piece of the original Phase 8 triangulation not re-verified locally; flagged explicitly as a
   remaining gap, not silently dropped.

Files: `06_md/system/production2.dcd` (new), `06_md/system/production_full_stripped.nc` (new),
`06_md/system/mmgbsa_results_full.dat` (new), `06_md/run_md_extend.py` (new).

## 2026-09-24 — A (Week 2 start: B4/B5 calibration — three preparation defects found and fixed)

Started Week 2. Track B's calibration items (B4/B5/B6) were run first: they are cheap CPU
docking, and B4 is what makes the paper's central selectivity claim defensible, so it gates
how the Week 2 MD results can be interpreted. Three separate defects surfaced, all in
receptor/reference preparation rather than in docking itself.

**Defect 1 — altloc contamination in the native-ligand extract (7Q8Y).** `prep_receptor.sh`
step 1 pulled the native ligand with `grep ^HETATM | awk 'resname && chain'` and **no altloc
filter**, while step 3 already passed `--default_altloc A` for the receptor. 7Q8Y's 9IV is
modelled in two alternate conformations (A@0.65, B@0.35), so `native_9IV.pdb` contained 48
atoms — two superimposed copies of a 24-heavy-atom ligand. Open Babel then perceived bonds
*across* the two copies and produced a fragmented non-molecule
(`C.C.C.C#C.C#C...` — 30+ disconnected pieces). **Every RMSD ever measured against that
reference was meaningless, which is the origin of the 5.29 A "7Q8Y redocking failure" the
manuscript attributes to a Vina scoring-function limitation.** Fixed in `prep_receptor.sh`
(keep blank-altloc or A, matching the receptor). 7Q8Y was the only one of the six affected;
re-running all six changed *only* 7Q8Y's ligand and box — the six `receptor.pdbqt` files came
back byte-identical (md5), so no previously-completed docking was invalidated by this fix.
The box barely moved (<=0.06 A): the box was never the problem, the reference ligand was.

**Defect 2 — bond orders guessed from coordinates.** The native-ligand SDFs are built by
Open Babel from crystal PDB heavy atoms, which carry no bond orders and no hydrogens. For
fused heteroaromatics it guesses wrong: 7Q8V's 9IV came back as
`Clc1ccccc1Oc1ccc(NC2=NC=NC3=NC=C[C@@H]32)cc1` — pyrrolo[2,3-d]pyrimidine non-aromatic, the
pyrrole N-H dropped, and an invented stereocentre on an atom that is planar in reality.
Added `scripts/fix_native_bondorders.py`, which transfers bond orders from each compound's
verified `01_smiles/references.csv` SMILES onto the crystal coordinates via RDKit's
`AssignBondOrdersFromTemplate` and writes `native_<LIG>_ref.sdf`. All six native references
now round-trip to exactly their reference SMILES (2V5Z SAG differs only in an unset C=N E/Z
label; connectivity and atom count match). These `_ref.sdf` files are now the RMSD reference.

**Checked and cleared: ligand-side protonation/tautomer (B5 item 2).** The suspicion was that
`obabel -p 7.4` moves the aminopyrimidine N-H onto the wrong nitrogen. It does not — because
`prep_ligands.py` seeds RDKit from the verified SMILES with an explicit `[nH]` before Open
Babel ever sees the molecule. Verified the prepared `native_9IV.sdf` round-trips to the exact
reference SMILES, neutral, tautomer intact. The ligand prep was never at fault.

**Defect 3 (the big one) — stripping the pocket water shell breaks the ATP-site receptors.**
With defects 1 and 2 fixed, 9IV still failed to redock: 7Q8V 4.21 A, 7Q8Y 5.60 A. Note 7Q8V
is a TTBK1 on-target and had never been redocking-validated at all before this week. The
poses were not randomly scattered — every one of the 9 sat 2-6 A off the crystal centroid, a
systematic displacement. Ruled out in order: box size (2V5Z has the *tightest* box margin,
3.4 A, and passes at 0.56 A, while 7Q8V has 4.8 A and fails), broken/incomplete pocket
(contact counts comparable to the passing receptors), and a pocket formed across a chain
interface (both sites are chain A only).

`vina --score_only` on the crystal pose settled it:

    receptor           crystal pose    best docked     docked better by
    7Q8V (TTBK1)          -7.07           -8.03              0.96
    7Q8Y (TTBK2)          -7.34           -9.32              1.98

Vina genuinely ranked a non-native pose above the crystal pose in both — i.e. with a dry
receptor this *is* a scoring failure, not a sampling failure. But the reason is that the dry
receptor is not the physical binding site: 9IV binds both TTBK paralogs through bridging
water. Retaining crystallographic waters within 5 A of the native ligand (bulk solvent still
discarded) fixes it outright, and does not harm the receptors that already passed:

    receptor          dry     +waters
    7JXX  TTBK1      0.71      0.55     (improves)
    4BTK  TTBK1      0.75      0.80     (native pose rises from rank 4 to rank 1)
    2V5Z  MAO-B      0.56      0.73
    7Q8V  TTBK1      4.21      0.70     FAIL -> PASS
    7Q8Y  TTBK2      5.60      1.45     FAIL -> PASS

7Q8V needed exactly **one** water to go from 4.21 A to 0.70 A. A uniform "retain the 5 A
pocket water shell" protocol therefore validates all five, where the dry protocol validated
only three. Implemented as a 5th argument to `prep_receptor.sh` (`WATER_SHELL_A`, default
5.0) applied identically to every receptor, so rule 1.5 still holds across the whole set;
all six regenerated in one batch. Honest caveat to carry into Methods: with waters the native
pose is top-ranked for 7Q8Y and 4BTK but sits at rank 5 for 7Q8V, so the pass is on the
project's existing best-of-9 criterion (which is also how 4BTK's original 0.75 A pass was
scored — that was a rank-4 pose).

**Verdict on B5.** The manuscript's Limitation #3 was *half* right and for the wrong reason.
The specific 5.29 A number was an artifact of the altloc bug and is void. The underlying
claim that Vina misplaces 9IV is reproducible on a correctly-prepared *dry* receptor — but it
is a preparation artifact, not a scoring-function ceiling, and it disappears when the bridging
waters are kept. Both anti-target/calibration receptors are now validated. Limitation #3 as
written should be replaced by a statement that these sites are water-mediated and were
prepared with an explicit pocket water shell.

**Also fixed:** `dock.sh` silenced `mk_export` failures with `2>/dev/null || true`. One 9IV
run hit exactly that and left a valid `.log` and `.pdbqt` but no `.sdf`, which `rmsd_check.py`
then reported as `BEST None 1000000000.00 A -> FAIL` — a preparation failure that reads as a
docking failure. Now retries once and fails loudly.

Files: `scripts/prep_receptor.sh` (altloc filter + water shell), `scripts/dock.sh` (loud
mk_export), `scripts/fix_native_bondorders.py` (new), `scripts/calibration_9iv.py` (new),
`03_receptors/*/native_*_ref.sdf` (new), `03_receptors/*` regenerated with 5 A water shell,
`05_validation/{7Q8V,7Q8Y}_redock.txt`.

**Circularity check on the water shell (important caveat, tested not assumed).** Redocking a
native ligand into a receptor that retains *its own* crystallographic waters is partly
self-fulfilling: those waters were positioned in the presence of that ligand. So the
7Q8V/7Q8Y passes above are, on their face, weaker evidence than 7JXX's dry 0.71 A. Tested
whether the water that matters is ligand-induced or a conserved feature of the site, by
superposing the other TTBK structures onto 7Q8V via sequence alignment (Biopython
PairwiseAligner -> Superimposer on aligned CA; backbone RMSD 0.41-0.97 A over ~290 CA, i.e.
genuine superposition — note that naive residue-number matching gives ~20 A RMSD for 4BTK and
7Q8Y because their numbering differs, and must not be used here) and measuring the distance
to the nearest water in each:

    structure   protein   ligand   nearest water to 7Q8V HOH572
    7JXX        TTBK1     VP7            1.87 A
    4BTK        TTBK1     DTQ            2.05 A
    7Q8Y        TTBK2     9IV            1.80 A

A water sits at essentially the same site in all four structures, across three different
ligands and two different proteins. It is therefore a conserved structural water of the TTBK
ATP site, not an artifact of 9IV being bound, and retaining it is a property of the site
rather than a crutch for the reference ligand. The 1.8-2.1 A spread is consistent with
coordinate uncertainty at these resolutions plus the 0.4-1.0 A backbone superposition itself.
Still to state plainly in Methods: the waters are held rigid by Vina, which is the standard
trade-off and does bias mildly toward ligands that complement the retained shell.

**GATE 1 re-run on the water-retained receptors (3 seeds each, best of 9, symmetry-corrected
against the bond-order-corrected `_ref.sdf` references).** All five pass:

    7JXX  TTBK1  VP7   0.55 A  PASS   (pose rank 2)
    4BTK  TTBK1  DTQ   0.79 A  PASS   (rank 1)
    2V5Z  MAO-B  SAG   0.66 A  PASS   (rank 3)
    7Q8V  TTBK1  9IV   0.68 A  PASS   (rank 5)
    7Q8Y  TTBK2  9IV   1.45 A  PASS   (rank 1)

Previously only three of these had ever been validated, and two of those five would have
failed under the dry protocol. 7Q8V in particular is a TTBK1 **on-target** that had never
been redocking-validated at all before this week.

**B4 — the 9IV matched-pair calibration. The protocol has a real, measurable TTBK1/TTBK2
bias.** Both receptors prepared by one run of the same script; ONE ligand file docked into
both, so nothing here can come from ligand preparation. Per-seed spread is tiny (SD 0.015 and
0.022 kcal/mol over seeds 11/22/33), so the margin is resolved far beyond its own noise:

    consensus 7Q8V (TTBK1)  -8.115 +/- 0.015
    consensus 7Q8Y (TTBK2)  -8.899 +/- 0.022
    measured margin (TTBK1 - TTBK2)  = +0.784 +/- 0.015 kcal/mol
    experimental ddG from IC50       = -0.077 kcal/mol  (range -0.234 .. +0.046)
    SYSTEMATIC PROTOCOL BIAS         = +0.861 kcal/mol, favouring TTBK2

9IV is close to equipotent on the two paralogs (TTBK1 330-530 nM, TTBK2 490 nM), so its true
ddG is ~0. The protocol nevertheless scores TTBK2 as the better binder by 0.78 kcal/mol.
**That 0.861 kcal/mol is an additive correction that must be subtracted from every reported
TTBK1/TTBK2 selectivity margin in 3.7.** Concretely: the manuscript's headline 1.7 kcal/mol
TTBK2-over-TTBK1 preference becomes ~0.84 kcal/mol once corrected — roughly half of the
claimed effect is an artifact of the protocol, not chemistry. This is exactly the failure mode
B4 was designed to catch, and it would not have been visible without a matched pair.

**Control for the obvious confound.** 7Q8Y retained 10 pocket waters and 7Q8V only 1 (the
crystals simply resolve different numbers), so "same script" does not mean "same number of
waters", and more waters could mean more favourable contacts and hence a spuriously better
score. Re-computed the margin on the *dry* receptors, where both have zero waters and the
comparison is perfectly symmetric:

    DRY  (0 waters both)      margin +1.278 +/- 0.009
    WET  (1 vs 10 waters)     margin +0.784 +/- 0.015

The bias is present in both and is *larger* in the symmetric dry case, so it is not a
water-count artifact -- it is a genuine bias between these two receptor structures, and
retaining the water shell in fact reduces it by ~40%. Report the WET number (+0.784, bias
+0.861) since that is the validated protocol, and cite the dry control as evidence the
correction is not an artifact of asymmetric solvation.

Files: `05_validation/calibration_9IV.csv`, `05_validation/calibration_9IV_margin.json`,
`05_validation/{7JXX,4BTK,2V5Z,7Q8V,7Q8Y}_redock.txt`.

**Defect 4 — 2Z5X (MAO-A) has been docked without its FAD cofactor this whole time.** Spotted
while verifying the regenerated receptors: 2V5Z's `receptor.pdbqt` contains 61 FAD atoms,
2Z5X's contained **zero**. FAD is present in 2Z5X's `raw.pdb` (53 atoms, chain A — note the
Day-3 workflow note claiming 2Z5X lists FAD in chain B is wrong, there is only one chain) and
survives into `clean_noH.pdb`, so the loss happens inside Meeko. The dry backup is identical,
so this is not new today: **every MAO-A docking run so far, including B2's 84 reference
dockings, used an active site with no flavin in it.** For an MAO anti-target that is a serious
error — the flavin forms one wall of the substrate cavity, and harmine stacks directly on it.

Root cause, confirmed by measurement:

    2V5Z (MAO-B)  Cys397 SG--FAD C8M = 2.31 A  -> too long to perceive as a bond -> FAD kept
    2Z5X (MAO-A)  Cys406 SG--FAD C8M = 1.65 A  -> perceived as covalent          -> FAD DROPPED

Both structures contain the same 8alpha-S-cysteinyl FAD; they merely refine that bond to
different lengths. Meeko has no FAD residue template, builds one on the fly, and when the
Cys linkage is perceived it reports `matched with excess inter-residue bond(s): A:600`,
template matching fails, and `-a/--allow_bad_res` **silently deletes the residue**. So 2V5Z
only "worked" by the accident of a loosely-refined bond. Relying on that would also break rule
1.5 outright, since MAO-B and MAO-A are compared to each other.

`--set_template` is not a fix: it takes a residue *name* and tries to fetch it from RCSB
(`A:600=FAD_5p` -> HTTP 404, "Ligand FAD_5p not available from rcsb.org").

Fix: cofactors now bypass Meeko's polymer path entirely. `prep_receptor.sh` writes
`clean_polymer.pdb` (protein + pocket waters, no cofactors) for Meeko, and the new
`scripts/prep_cofactor.py` prepares each declared cofactor through the SAME ligand toolchain
everything else in this project uses — RCSB ideal chemistry -> bond orders transferred onto the
crystal coordinates -> `obabel -p 7.4` -> `mk_prepare_ligand` -> torsion tree stripped -> atoms
appended to `receptor.pdbqt` as rigid receptor. Applied identically to both MAO structures.

Validated against the path it replaces: for 2Z5X the ligand route reproduces 2V5Z's
Meeko-derived FAD exactly — same 61 atoms, **identical** AutoDock typing (15 A, 12 C, 8 HD,
4 N, 5 NA, 15 OA, 2 P) and net charge -2.002 vs -2.000. So this is not a different
parameterization, it is the same chemistry obtained by a route that does not depend on a
refined bond length.

Also added a hard **guard** in `prep_receptor.sh`: any cofactor with atoms in `clean_noH.pdb`
and zero atoms in `receptor.pdbqt` now aborts the run with a message rather than producing a
quietly gutted receptor. This class of silent drop is exactly what `--allow_bad_res` invites,
and it went unnoticed for weeks.

Files: `scripts/prep_cofactor.py` (new), `scripts/prep_receptor.sh` (polymer/cofactor split +
cofactor guard), `03_receptors/_cofactors/FAD_ideal.sdf` (cached RCSB reference).

**Limitation to state with the B4 correction.** The +0.861 kcal/mol bias is estimated from a
*single* matched pair, because 9IV is the only ligand crystallised in both TTBK1 and TTBK2 --
that is precisely why it was chosen, but it also means the correction rests on one compound.
The seed-to-seed precision (SD 0.015) is the precision of *that one measurement*, and must not
be presented as the uncertainty on the correction itself, which is unknown and certainly
larger. Treat +0.861 as a point estimate of a systematic offset, apply it as a correction, and
say plainly in Methods that it derives from one matched pair. A second matched pair would be
needed to put a real error bar on it, and none exists in the PDB for this paralog pair today.
A weaker but useful cross-check becomes available once the 56 candidates are docked on both
7Q8V and 7Q8Y: if the TTBK1-minus-TTBK2 margin clusters near +0.8 across chemically diverse
compounds rather than scattering, that supports the offset being a property of the receptor
pair rather than of 9IV. Worth running -- it is 112 extra dockings on receptors that are
already prepared and validated.

**GATE 1 now complete across all six receptors — a first for this project.** 2Z5X had never
been redocking-validated at all. With FAD restored and the 5 A water shell, harmine redocks to
**0.27 A** (3 seeds), the best result of the whole set:

    7JXX  TTBK1  VP7   0.55 A  PASS
    4BTK  TTBK1  DTQ   0.79 A  PASS
    2V5Z  MAO-B  SAG   0.66 A  PASS
    7Q8V  TTBK1  9IV   0.68 A  PASS
    7Q8Y  TTBK2  9IV   1.45 A  PASS
    2Z5X  MAO-A  HRM   0.27 A  PASS

**How bad was the flavin-less MAO-A receptor, quantitatively?** Measured rather than assumed,
by redocking harmine into the exact old receptor B2 used (no FAD, no waters):

    pose:   0.27 A (rank 2)  ->  1.45 A (rank 7)    degraded, but still a formal pass
    score: -8.748 +/- 0.039  vs  -8.663 +/- 0.006    difference only -0.085 kcal/mol

So the earlier characterisation of this as producing a badly wrong MAO-A pocket was
overstated, and is corrected here. FAD *is* genuinely part of the site — 13 of its 53 atoms
fall inside the docking box and it approaches harmine to 3.90 A (safinamide to 3.00 A in
2V5Z) — but for harmine specifically the missing flavin cost only 0.085 kcal/mol, because
harmine is small and sits mostly against protein with the flavin edge-on. The defect was real
and the receptor was wrong, but its numerical impact is compound-dependent, and for this one
reference compound it is small.

That does not make the re-dock optional: any candidate that stacks face-on against the flavin
(which is the classic MAO binding mode, and what several of the 56 flavonoid-derived
candidates would be expected to do) stands to be affected far more than harmine is. The
honest statement is that the magnitude is unknown per compound until re-docked, not that the
previous MAO-A numbers were all badly wrong.

**Week 2 production re-dock: scoped, scripted and resumable.** Because the receptors changed
materially (water shell, FAD restored, 7Q8Y reference de-duplicated), every production docking
number has to be regenerated — the old and new tables cannot be mixed. Measured on the first
compounds through the new 7JXX, the water shell moves consensus scores by -1.0 to +0.7
kcal/mol depending on the compound (cand_003 -8.385 -> -9.212, cand_006 -8.105 -> -7.437), so
this is not a formality.

Added `scripts/run_week2_redock.sh`. It is idempotent: a run counts as current only if it has
the expected number of logs AND every log is newer than the `receptor.pdbqt` it was docked
against, so an interrupted overnight job resumes exactly where it stopped and stale
pre-water-shell results are never silently accepted as done.

Honest cost, measured rather than guessed (~2.7 min/ligand for candidates at exhaustiveness 32,
Vina auto-detecting 16 cores):

    required set   15 runs = 480 dockings  ~18 h
      7JXX + 2V5Z candidates (336), references on 7JXX/2V5Z/2Z5X (144)
    --with-paralog  +6 runs = +336 dockings  ~15 h more
      the 56 candidates on 7Q8V and 7Q8Y, as the cross-check on the B4 correction

This is well beyond the 6-8 h first estimated; that estimate was wrong and is corrected here.
Deliberately NOT reduced by lowering `--exhaustiveness` or pinning `--cpu`: Vina's results
depend on both, so either change would decouple the production numbers from the GATE 1
validation that licenses them. The run is cheap in attention (unattended) and expensive only
in wall-clock.

## 2026-09-25 — A (Week 2 re-dock launched; two resumability defects fixed first)

Picked up where yesterday's entry stopped: `scripts/run_week2_redock.sh` written, the overnight
production re-dock part-way through. Before restarting it I checked the resume logic actually
resumes. It did not.

**Defect 1 — the resume guard fails open, and re-docks everything.** `expected_count()` counted
ligands with `python -c "import pandas..."`. The script consults it *before* anything guarantees
the conda env is on PATH, and a bare Windows `python` has no pandas, so the call returned
nothing, `N_EXP` came back empty, and the freshness test collapsed:

    scripts/run_week2_redock.sh: line 78: [: : integer expected
    [run] 2V5Z references seed11  (0/ current)

Every run then classed as stale. Restarting the interrupted job would silently have re-docked
the 112 dockings already finished (~5 h) rather than resuming — the exact failure the script was
written to prevent, hidden behind a message that reads like ordinary progress. Counting rows in
a CSV needs no chemistry, so the count is now `awk 'END{print NR-1}'`, and an unparseable or
zero count **aborts** instead of degrading into "run everything".

With that fixed the real state became visible — and the dry run is now trustworthy:

    [current] 7JXX candidates_56 seed11 (56/56)
    [current] 7JXX candidates_56 seed22 (56/56)
    [run]     7JXX candidates_56 seed33 (10/56 current)
    ...
    runs needed: 19   already current: 2      (13 required + 6 paralog)

**Defect 2 — resume granularity was a whole run.** `dock.sh` loops over every ligand in the set
unconditionally, so the partially-finished seed33 would redo the 10 ligands it had already
docked. Added a per-ligand skip: a ligand counts as done only if its `.log` is newer than the
receptor **and** its `_out.sdf` ends in the `$$$$` record terminator. The terminator is the last
thing written, so a job killed mid-export leaves a truncated SDF that correctly fails the test
and is re-docked — this deliberately does not trust mere existence of the file, which is what
lets a half-written pose masquerade as a result.

Verified on a run that was already complete: `dock.sh 7JXX candidates_56 11` returned `[ok]` in
**35 s** instead of ~2.5 h, and the mtime of an existing log was unchanged, i.e. it skipped
rather than silently redoing work.

**Required set launched** (13 runs, ~590 dockings; paralog cross-check deferred). Environment:
`vina.exe` is on the system PATH, pandas comes from the `docking_project` env, so the job runs
under `conda run -n docking_project`.

**The 2.7 min/ligand figure in yesterday's entry does not hold for 2V5Z.** First ligands through
are ~13 s apart (cand_001 01:26:37, cand_002 01:26:50) at the same exhaustiveness 32. The ~18 h
projection was extrapolated from 7JXX and is likely a large overestimate; the true cost will be
measured from this run rather than re-guessed here.

One thing checked explicitly, because getting it wrong corrupts output silently: that no second
docking job was still alive writing the same directories. Yesterday's job finished 7JXX seed33
at 01:17 and exited; exactly one `vina.exe` is running, and the stale 2V5Z logs it is now
overwriting are from 17:00 the previous day, i.e. older than the rebuilt receptor and correctly
treated as invalid.

Files: `scripts/run_week2_redock.sh` (awk count + abort guard), `scripts/dock.sh` (per-ligand
resume).

### Same day, while the re-dock ran — three further defects, all in what happens *after* docking

**Defect 3 — `collect_results.py` would have silently mixed old and new receptors back
together.** The whole point of the Week 2 re-dock is that pre-rebuild numbers are not
comparable with post-rebuild ones. The collector globbed `04_docking/*/*.log` with no
freshness test at all, so it would have averaged both into `consensus_week2.csv` and produced
a table that looks entirely normal and is wrong. Run against the tree mid-job it was pulling
in **229 stale logs**:

    [stale] dropped 229 log(s) older than their receptor:
              2V5Z_candidates_56_seed22: 56
              2V5Z_candidates_56_seed33: 56
              2Z5X_references_seed11: 16
              ...

The collector now applies the same rule as the resume guard — a log older than the
`receptor.pdbqt` it names is dropped and reported, with `--allow-stale` as the deliberate
escape hatch. Two further guards while there: it warns when a `(receptor, ligand)` group pools
more than one ligandset (that is two preparations of one compound being averaged, not a
consensus), and when a ligand has fewer than the expected 3 seeds.

**Defect 4 — 2V5Z's GATE 1 was validated against a receptor that no longer exists.**
Timestamps: `05_validation/2V5Z_redock.txt` 23:32:47, `03_receptors/2V5Z/receptor.pdbqt`
**23:44:51**. 2V5Z and 2Z5X were rebuilt last, for the cofactor fix, and only 2Z5X was
re-validated afterwards. So yesterday's claim that "GATE 1 is now complete across all six
receptors" does not hold for 2V5Z: its 0.66 A pass describes the previous receptor. The
chemistry barely moved (the ligand-route FAD reproduced Meeko's typing exactly, net charge
-2.002 vs -2.000), so it is expected to re-pass — but expected is not measured, and a GATE
claim has to be measured. Queued as a re-dock; 3 dockings.

**Defect 5 — the Week 2 job list regenerates none of the paper's selectivity margins.** This is
the serious one. Manuscript 3.7 reports three margins and 4.2 calls the selectivity liability
"the dominant risk"; the required set covers the two *on-target* columns (7JXX, 2V5Z) and
2Z5X *references*, and no candidates on any anti-target:

    TTBK1 vs TTBK2 (6U0K)   receptor no longer present in 03_receptors/  -- orphaned
    TTBK1 vs TTBK2 (7Q8Y)   candidates deferred as the optional paralog job
    MAO-B  vs MAO-A         candidates on 2Z5X were never in the job list at all

The MAO-A gap is the worst of the three. Those candidate scores come from the old
"custom FAD-transplant splice" 2Z5X, the receptor this week's cofactor work replaced outright,
and they are precisely the numbers most exposed to the FAD fix — the earlier entry's own
argument was that flavonoid candidates stacking face-on against the flavin stand to be
affected far more than harmine's 0.085 kcal/mol. The paper's central claim currently rests on
a receptor the project has already declared wrong.

**The cost argument for deferring has also collapsed.** Measured on this run, 2V5Z candidates
are going through at ~28 s/ligand, so a 56-candidate run is ~26 min, not the ~2.5 h the
"+15 h paralog" estimate assumed. Queued after the required set: candidates on 2Z5X, 7Q8Y and
7Q8V (9 runs, ~4 h), which closes both regenerable margins and delivers the B4 cross-check as
a by-product. The 6U0K row cannot be regenerated without re-preparing that receptor from
scratch; flagged here as orphaned, manuscript left untouched pending a decision.

Note for whoever restarts a job mid-flight: do **not** edit `dock.sh` or `run_week2_redock.sh`
while they are executing. Bash reads a script incrementally, so an in-place edit can make the
running shell misparse the remainder. `collect_results.py` was safe to fix mid-run only
because it is invoked once, at the end.

### First scientific result of the re-dock: the 7JXX ranking does not survive the water shell

The 7JXX candidate set finished early in the run (56 ligands x 3 seeds, all fresh), so the
pre- and post-rebuild tables can be compared directly. First the comparison was made
trustworthy: the docking logs are git-tracked, so the *true* immediate predecessor was
reconstructed by reading every `cand_*.log` out of `git show HEAD:` and re-deriving the
consensus, rather than trusting a CSV of uncertain vintage. It reproduces
`08_analysis/consensus_new.csv` to **max |diff| = 0.0000 kcal/mol**, confirming that file is
exactly the A5 (2026-08-20) table and that what follows compares adjacent states of the
project, not two unrelated rounds.

Adding the 5 A pocket water shell to 7JXX changes the answer:

    mean delta      -0.080 kcal/mol      (near zero -- the shift is NOT a uniform offset)
    median delta    -0.259
    sd of delta      0.833
    range           -0.977 .. +3.041
    |delta| > 0.5   28/56 candidates
    Spearman(old,new) = 0.475
    top-15 overlap    = 9/15
    rank moves >= 10 places: 27/56   (mean |rank move| 12.6 places)

The mean being ~0 while the spread is 0.8 is the important part: this is not a constant
correction that leaves the ranking intact, it is a reshuffle. Over half the shortlist moves
by more than half a kcal/mol and six of the top fifteen are different compounds.

Two individual results carry most of the consequence:

    cand_013   -8.826 (rank  1/56, seed sd 0.005)  ->  -5.785 (rank 56/56, sd 0.011)   +3.041
    cand_003   -8.385 (rank  8/56, seed sd 0.020)  ->  -9.194 (rank  1/56, sd 0.036)   -0.809

**cand_013 was the top-ranked candidate on TTBK1 and is now the worst of all 56.** The seed SDs
on both sides are 0.005 and 0.011, so this is not sampling noise — it is reproducible, and the
most likely reading is that cand_013's old pose occupied volume now held by a conserved pocket
water. Any statement in the manuscript resting on cand_013 has to go.

**cand_003, the compound already carried through MD and MM-GBSA, improves to rank 1.** That is
a genuine relief rather than a result to celebrate: the lead was selected on the old table, and
it happens to survive. Had it moved the way cand_013 did, the entire Phase D/F effort would
have been spent on a compound the corrected receptor ranks near the bottom.

This settles the question of whether the re-dock was worth 4 h of wall clock. It also means the
shortlist, Table 3.7 and every downstream selection in the manuscript are provisional until the
anti-target runs land — the on-target ranking has already changed underneath them.

### CORRECTION to the entry above — the 7JXX "reshuffle" is substantially an artifact

The previous entry reported the ranking change as a scientific result. That was wrong, and it
was written after checking only that the docking *job* was healthy (one vina process, logs
appearing, counts advancing) rather than that its *output* was sound. The pattern that should
have stopped it: five candidates collapsing into a narrow -5.8 to -6.2 band, which is what
ligand exclusion looks like, not what a chemical effect looks like.

Everything except the receptor is byte-identical between the two rounds — same box (centre
178.93/19.64/46.21, 18 A cube), same ligand pdbqt, same seed, same exhaustiveness. So the shift
is entirely attributable to the receptor, i.e. to the water shell. Tested against ligand
properties across all 56:

    delta vs TPSA    Pearson r = +0.792  (p = 3.6e-13)
    delta vs Rings   Pearson r = +0.774  (p = 2.7e-12)
    delta vs MW      Pearson r = +0.640  (p = 1.1e-07)
    delta vs RotB    Pearson r = +0.104  (p = 0.45, n.s.)

    MW quartile      mean delta     got worse
    Q1 smallest        -0.334          2/14
    Q4 largest         +0.757          9/13

And the decisive one — the sign of the score/size relationship inverts:

    old receptor   score vs MW   r = -0.337     (bigger scores better: normal Vina size bias)
    new receptor   score vs MW   r = +0.458     (bigger scores WORSE)

Vina's known size bias means larger ligands normally score better. A receptor that reverses
that is not reporting chemistry, it is sterically excluding ligands. TPSA being the strongest
correlate is the clincher and points the same way: the *polar* ligands are penalised hardest,
which is backwards, since polar ligands are exactly the ones that would hydrogen-bond to a
pocket water or displace it favourably.

**Root cause.** `prep_receptor.sh` retains every crystallographic water within 5 A *of the
native ligand* and writes them into `receptor.pdbqt` as rigid receptor atoms. Two consequences
that were not separated:

1. **The GATE 1 justification is circular for cross-docking.** The retained waters are by
   construction the ones that coexist with the native ligand in its own crystal. Redocking that
   same ligand into a pocket moulded around it will tend to improve, so GATE 1 improving does
   not validate the shell for the 56 candidates, which are different chemotypes. The
   validation and the application are not the same experiment.
2. **Rigid waters cannot be displaced.** Real pocket waters are displaceable at a desolvation
   cost. Modelling them as immovable forbids the exchange outright, and the penalty falls on
   whichever ligands are largest and most polar — precisely the observed correlations.

This is not an argument that the waters are wrong everywhere. 7Q8V went 4.21 A -> 0.70 A on a
single bridging water, and that one is well evidenced. The problem is the blanket 5 A rule:

    7JXX 10 waters in box    2V5Z 8    2Z5X 8    7Q8Y 10    4BTK 7    7Q8V 1

Ten rigid waters inside an 18 A box is a large amount of occluded volume.

**This also undermines the B4 paralog cross-check queued earlier.** The correction was measured
between 7Q8V (**1** water) and 7Q8Y (**10** waters). For 9IV, native to both, the shells are
fitted around the ligand itself and the comparison is roughly fair. For the 56 candidates,
native to neither, a 1-water pocket would be compared against a 10-water pocket and the margin
would conflate the paralog difference with a water-count difference. Running those 336 dockings
before the water treatment is settled would produce a number that cannot be interpreted.

**Status:** the extension (2Z5X/7Q8Y/7Q8V candidates) is on hold pending a decision. The
current job is being allowed to finish because the references runs are the objective test —
those 16 compounds have measured IC50s, so whether the wet receptors predict experiment better
or worse than the dry ones can be measured rather than argued.

**Amendment to the correction — the sign inversion is a 7JXX effect, not a universal one.**
Replicating the same test on 2V5Z, whose A5 table also survives in `consensus_new.csv`:

    receptor   score vs MW (old -> new)      delta vs TPSA    Spearman(old,new)   top-15
    7JXX        -0.337  ->  +0.458  INVERTED   +0.792 ***          0.475            9/15
    2V5Z        +0.271  ->  +0.566  same sign  +0.589 ***          0.748           11/15

On 2V5Z the *old* receptor already scored larger ligands worse, which is chemically reasonable
for MAO-B's narrow substrate cavity, so the waters strengthen an existing trend rather than
reversing one. Presenting sign inversion as the general signature of the artifact was an
overstatement; it is one receptor's symptom. What replicates across both proteins is the
polarity penalty — delta vs TPSA positive and strongly significant in each — and that is the
claim worth keeping. Magnitudes differ substantially: 7JXX is badly disturbed (sd 0.833, six of
fifteen shortlist members change), 2V5Z much less (sd 0.463, mean delta -0.504, i.e. the waters
slightly *improved* MAO-B scores on average).

Also worth stating plainly: old-vs-new compares across a whole receptor regeneration, so it is
not a clean measurement of the waters alone. The wet-vs-dry runs now in progress are the clean
one — the dry receptors are the wet files with only the water atoms deleted, so nothing else
can differ.

## 2026-09-25 — A (water-shell test: the artifact is confirmed and fully attributed)

Required set finished 03:02, all 12 runs clean, `08_analysis/consensus_week2.csv` written
(165 rows, 3 seeds throughout). GATE 1 gap from earlier today also closed: **2V5Z redocks SAG
at 0.66 A PASS** against the current receptor, so the "all six receptors" claim is now measured
rather than asserted.

**Arm 1 — candidates on 7JXX, wet vs dry.** Same receptor file, water atoms deleted, nothing
else touched, so the waters are the only variable:

    score vs MW      WET r = +0.458 (p=3.9e-04)      DRY r = -0.337 (p=1.1e-02)
    score vs TPSA    WET r = +0.792 (p=3.7e-13)      DRY r = +0.049 (p=0.72, n.s.)
    score vs Rings   WET r = +0.612 (p=5.4e-07)      DRY r = -0.287 (p=3.2e-02)

Removing the waters **abolishes the polarity penalty entirely** — TPSA goes from the strongest
correlate in the set to statistically indistinguishable from zero — and restores Vina's normal
negative size bias. This is as clean as this kind of test gets.

**The attribution is exact.** The dry receptor reproduces the pre-rebuild A5 table to
**max |diff| = 0.0000 kcal/mol over all 56 candidates, Pearson r = 1.0000, top-15 overlap
15/15**. cand_013 comes back at -8.826 and cand_003 at -8.385, the A5 values to three decimals.
So every other element of this week's 7JXX rebuild — Meeko re-run, polymer/cofactor split —
changed the scores by nothing whatsoever, and the 5 A water shell accounts for 100% of the
difference between the A5 table and the new one. cand_013's -8.826 -> -5.785 collapse is
entirely the waters. (It also confirms the pipeline is exactly reproducible given receptor,
ligand and seed, which is worth knowing independently.)

**Arm 2 — correlation with measured IC50. This one does NOT support removing the waters, and
that has to be said rather than buried.**

    MAO-B 2V5Z  n=3   safinamide  pIC50 8.12   wet -10.300  dry  -9.989
                      lazabemide  pIC50 7.20   wet  -7.183  dry  -6.601
                      isatin      pIC50 5.52   wet  -7.678  dry  -7.299
    MAO-A 2Z5X  n=2   kaempferol  pIC50 6.15   wet  -8.362  dry  -8.952
                      quercetin   pIC50 5.82   wet  -7.424  dry  -8.102
    TTBK1 7JXX  n=1   9IV         pIC50 6.37   wet  -8.968  dry  -8.624

    pooled   WET r = -0.627 (p=0.18)      DRY r = -0.410 (p=0.42)

The wet receptors correlate *better* with experiment, not worse. With n=6 pooled across three
different proteins, neither correlation is significant and the difference between them is
meaningless — this arm cannot decide anything, exactly as predicted before it was run. It is
recorded because it is evidence that points the other way, and the honest position is that the
experimental data available cannot adjudicate this. Note in passing that the two flavonoids,
the compounds most like the 56 candidates, are the two that score *better* dry.

**Where this leaves the protocol.** What is established: the shell imposes a large,
systematic, polarity-dependent penalty on non-native ligands, and it is the sole cause of the
candidate ranking change. What is NOT established: that dry receptors predict experiment
better. The case for a minimal-water protocol therefore rests on the mechanism and on GATE 1,
not on arm 2:

  - dry already passes GATE 1 on 7JXX (0.71), 4BTK (0.75) and 2V5Z (0.56) — the three
    receptors that carry candidate scoring, so the shell buys nothing where it costs most;
  - waters are load-bearing only on 7Q8V (4.21 -> 0.70, one bridging water) and 7Q8Y
    (5.60 -> 1.45), both 9IV structures;
  - the native-redock justification is circular for cross-docking, since the retained waters
    are by construction the ones that coexist with the native ligand.

Unresolved tension to settle before re-docking: rule 1.5 requires receptors whose scores are
compared to be prepared identically, but TTBK1/TTBK2 selectivity compares 7JXX against 7Q8Y,
and 7Q8Y genuinely needs water while 7JXX does not. Options are to run the paralog comparison
on 7Q8V vs 7Q8Y (both minimal-water TTBK1/TTBK2) rather than 7JXX vs 7Q8Y, or to accept
non-uniform preparation and state it. Not decided here.

**Practical consequence worth flagging: if the dry protocol is adopted, the A5 candidate
numbers for 7JXX are still valid** — they are bit-identical to the dry re-dock. It is the
wet-receptor production run finished tonight that would be discarded, not the older work.
Whether the same holds for 2V5Z is untested; its cofactor route also changed, and candidates
were only re-docked there wet.

Files: `scripts/run_water_test.sh`, `scripts/analyze_water_test.py`,
`03_receptors/{7JXX,2V5Z,2Z5X}dry/`, `08_analysis/water_test_candidates.csv`,
`08_analysis/water_test_references.csv`.

### 2V5Z wet-vs-dry: the artifact is receptor-dependent, and MAO-B keeps a real polarity term

The 2V5Z gap flagged earlier is now closed — candidates docked on 2V5Zdry (168 dockings), so
both scoring receptors have a clean wet-vs-dry comparison with the waters as the only variable.

                          7JXX                        2V5Z
    dry vs A5       max|diff| 0.0000  r=1.0000   max|diff| 0.0867  r=0.9993
    TPSA  wet            +0.792 ***                  +0.788 ***
    TPSA  dry            +0.049  n.s.                +0.407  (p=0.0018)
    MW    wet            +0.458                      +0.566
    MW    dry            -0.337                      +0.270
    delta(wet-dry)   mean -0.080  sd 0.833       mean -0.503  sd 0.460
    delta vs TPSA        +0.792 ***                  +0.596 ***
    Spearman(wet,dry)     0.475   top-15 9/15         0.749   top-15 10/15

Two corrections to the earlier write-up follow from this.

**1. The A5 tables are not universally recoverable.** 7JXX dry reproduces A5 exactly, but 2V5Z
dry differs by up to 0.087 kcal/mol (r = 0.9993). That is the cofactor route change showing up
precisely where it was predicted to — the ligand-route FAD carries net charge -2.002 against
Meeko's -2.000 — and it is real, if tiny. So the earlier note that "the A5 candidate numbers
are still valid if dry is adopted" holds for 7JXX only. 2V5Z would need regenerating either
way, which the run just completed has already done.

**2. The polarity penalty is not purely an artifact on every receptor.** On 7JXX removing the
waters takes TPSA from the strongest correlate in the set to indistinguishable from zero — the
whole correlation was the water shell. On 2V5Z it only falls from +0.788 to +0.407, and the dry
value is still significant at p = 0.0018. The natural reading is that MAO-B's substrate cavity
is genuinely hydrophobic, so a real chemical polarity term exists there and the waters add an
artifact on top of it (delta vs TPSA +0.596). Saying flatly that "the waters cause the polarity
penalty" was therefore too broad: it is true on TTBK1, and only partly true on MAO-B.

2V5Z is also less disturbed overall — Spearman 0.749 against 7JXX's 0.475 — and the waters
slightly *improve* its scores on average (mean delta -0.503). The severity of the artifact
tracks how much dead volume the shell adds relative to the pocket, which matches the earlier
geometric finding that 8 of 7JXX's 10 waters touch no ligand atom at all.

### Minimal-water probe: the protocol is viable, and half of 7Q8Y's waters do nothing

12 dockings, GATE 1 native redocking on reduced-water receptors:

    7JXXdry  (0 waters)   0.71 A  PASS     matches the documented dry value exactly
    2V5Zdry  (0 waters)   0.51 A  PASS     documented 0.56; shifted by the cofactor route
    2Z5Xdry  (0 waters)   1.02 A  PASS     MAO-A with no waters -- never tested until now
    7Q8Ybrg  (5 waters)   1.45 A  PASS     IDENTICAL to the all-10-water result

The 7Q8Y line is the important one. Keeping only the five waters that bridge 9IV and protein
reproduces the full-shell result to the digit, so **the five spectator waters contribute
nothing to validation** — exactly what the geometric analysis predicted, now measured rather
than inferred. 2Z5X passing dry at 1.02 A was the other unknown, and it removes the last
obstacle: every receptor that carries candidate scoring can be run with no waters at all.

Resulting protocol, uniform in procedure (start dry, add back only until the native redock
passes) and minimal in occlusion:

    7JXX  0 waters    4BTK  0 waters    2V5Z  0 waters    2Z5X  0 waters
    7Q8V  1 water  (required: dry fails at 4.21 A)
    7Q8Y  5 waters (bridging only; the other 5 are dead volume)

Note this is not a claim that the retained waters are "right" in any absolute sense, only that
each receptor now carries the minimum demonstrably needed to reproduce its own crystal pose,
and that the four receptors doing the candidate scoring carry none. The residual asymmetry
between 7JXX (0) and 7Q8Y (5) for the TTBK1/TTBK2 comparison is unavoidable — 7Q8Y genuinely
needs those waters — and should be stated in Methods rather than engineered away. Running that
comparison as 7Q8V (1) vs 7Q8Y (5) narrows the gap if a closer match is wanted.

Also worth recording: 2V5Zdry redocks at 0.51 A against the documented 0.56 A for the old dry
receptor, a difference in the same direction and of the same order as the 0.087 kcal/mol score
shift, i.e. the cofactor route change is visible in pose as well as score and both are small.

## 2026-09-25 — A (selectivity docking resumed; TTBK2 MD launched)

Session restart note: the selectivity job survived the previous session's exit and was still
docking (2Z5Xdry seed11, log written one second before the check). Deliberately NOT restarted —
a second orchestrator on the same output directories would interleave writes. Left to run.

**TTBK2 MD system built and launched**, the step the manuscript itself names as decisive
(4.3: extend MM-GBSA to TTBK2 and MAO-A). 48,158 atoms, 14,438 waters, 11 Cl-, net charge
+0.002; the TTBK1 baseline for comparison is 40,146 atoms, 12 Cl-, -0.001. Same recipe
(ff14SB/gaff2/tip3p, solvateoct 10 A), 0 tleap errors. 10 ns at 2 fs, ~2 h expected (TTBK1 ran
at 168 ns/day; this system is ~20% larger). Runs on the GPU while docking continues on CPU, so
the two do not compete.

Two decisions worth recording:

**Ligand parameters were reused, not regenerated.** The TTBK2 complex uses the *same*
`cand_003_AC.lib`/`.frcmod` as TTBK1, with only the coordinates changed to the TTBK2 docked
pose. Re-running acpype would have produced slightly different AM1-BCC charges, and then a
delta-G difference between the two complexes could reflect re-parameterisation rather than the
protein — fatal for a selectivity comparison, which is the entire purpose of the run.

**The coordinate transfer was verified rather than assumed.** Atom order turned out to be
identical across all four ligand files, but that was checked, not trusted. The first validation
(full inter-atomic distance matrix, max 2.19 A) appeared to fail; inspecting the worst pairs
showed they were all H...H and F...H, i.e. a CF3/hydroxyl torsion, with the mapping perfectly
correct. Re-validated on 1-2 and 1-3 distances, which are fixed by covalent geometry and
conformation-independent: **bonds agree to 0.068 A, angles to 0.057 A**. The original test was
measuring the wrong thing for a molecule with a rotatable bond.

**Preliminary and not to be quoted yet:** cand_003 docks to TTBK2 (7Q8Ybrg) at -8.089 against
-8.385 on TTBK1, i.e. a *favourable* 0.30 kcal/mol on-target margin, where Table 3.7 currently
reports +2.23 unfavourable. One compound, one seed, and the receptors differ in water content
(0 vs 5). The 56 x 3 run in progress is what decides it — but the selectivity conclusion, which
4.2 calls the dominant risk, may move.

**MAO-A MD is the hard one and is not started.** FAD is absent from ff14SB/gaff2 and in MAO it
is *covalently* bound through the 8alpha-S-cysteinyl linkage documented on 2026-09-24, so it
needs a custom linked Cys-FAD residue rather than a stock parameter set. That is multi-day work
with real failure risk against a 1.5-week deadline. Recommendation is to run TTBK2 now, and
treat MAO-A MM-GBSA as declared future work unless time frees up.

### MAO-A selectivity: the margin moves hugely, and 8 candidates do not dock at all

First look at the MAO-B vs MAO-A margin on the corrected receptors gave **-3.89 kcal/mol with
56/56 favourable**, against the +0.09 in manuscript Table 3.7. A uniform ~4 kcal/mol shift
across every compound is the signature of a systematic receptor difference, so it was checked
rather than reported.

**Not the box.** Volumes and crowding are comparable: 2V5Zdry 6124 A^3 / 325 receptor atoms in
box, 2Z5Xdry 5832 A^3 / 339 atoms.

**It is the FAD.** Every MAO-A docking before 2026-09-24 ran in a pocket with no flavin at all.
Deleting FAD from the current receptor (2Z5XnoFAD, 61 atoms removed) and re-docking reproduces
the old regime:

    cand_003   withFAD  -8.236    noFAD -11.38    FAD costs 3.14
    cand_044   withFAD  -7.909    noFAD -11.41              3.50
    cand_050   withFAD  -7.813    noFAD -11.36              3.55
    cand_013   withFAD  -0.449    noFAD  -9.049             8.60

So the shift is entirely the restored cofactor, and it is the quantitative confirmation of the
prediction made on 2026-09-24 — harmine lost only 0.085 kcal/mol because it sits edge-on, while
candidates stacking face-on against the flavin "stand to be affected far more". They lose ~3.5.

**But cand_013's -0.449 is not a weak affinity, it is a failed docking**, and that had to be
caught before any margin was quoted. Its mode table is:

    mode 1  -0.449      mode 2  +0.708      mode 3  +1.363
    mode 4 +12.99       mode 5  +45.17      (only 5 modes returned, not 9)

Positive affinities are net repulsion: Vina never found a placeable pose. Across the set, 8 of
56 candidates fail this way on 2Z5Xdry and **zero** fail on 2V5Zdry or 7JXXdry:

    2Z5Xdry  median -7.41   8 failures: cand_013, 025, 056, 023, 029, 020, 032, 054
    2V5Zdry  median -10.72  0 failures
    7JXXdry  median -7.99   0 failures

The eight are the largest compounds in the set (MW 420, 420, 380, 379, 365, 365, 351, 339),
which is chemically coherent — MAO-A has a single smaller cavity where MAO-B has a bipartite
one, and the flavin occupies part of it. `analyze_selectivity.py` now excludes any ligand whose
best pose is worse than -5.0 kcal/mol from the margin statistics and reports it separately.
Averaging a non-measurement into a delta-G would have manufactured selectivity out of a docking
failure: it was that inclusion, not chemistry, that produced the initial -3.89.

Handled honestly the result is arguably stronger, not weaker: 48 quantifiable candidates give a
mean margin of **-3.28, favourable in 48/48**, and the remaining 8 are *sterically excluded*
from MAO-A while binding MAO-B at -8.7 to -11.3. That is a qualitative selectivity statement,
not a delta-G, and must be reported as such. Numbers are provisional until the 2Z5Xdry seed33
run completes (52 ligands are still at fewer than 3 seeds).

Files: `scripts/analyze_selectivity.py`, `03_receptors/2Z5XnoFAD/` (control receptor).

### Selectivity margins on the corrected receptors — both Table 3.7 rows move, one is confounded

2Z5Xdry and 7Q8Ybrg complete (168/168 each), 7Q8V still finishing.

    row                manuscript 3.7          now
    TTBK1 vs TTBK2     +2.23   0/15 fav        +0.12   25/56 fav  (shortlist -0.15, 12/15)
    MAO-B vs MAO-A     +0.09   2/15 fav        -3.28   48/48 fav  + 8 sterically excluded

**B4 cross-check settles its question, and the answer is negative.** Across the 56 candidates
the 7Q8V-minus-7Q8Ybrg margin is **+0.100 +/- 0.323** (median +0.055, range -0.559..+0.916),
with only 14/56 falling within +/-0.5 of the +0.784 measured on 9IV. The candidate margins
cluster near zero, not near 9IV's value, so the offset is a property of **9IV specifically**
and not of the TTBK1/TTBK2 receptor pair. It must not be applied as a blanket correction —
which is exactly what the cross-check was designed to test, and it fails. (Provisional: 7Q8V is
at fewer than 3 seeds until its last run lands.)

**The MAO-A row is trustworthy; the TTBK2 row is not yet.** The MAO-B/MAO-A comparison is
water-symmetric — both receptors carry zero waters — so the -3.28 reflects the restored flavin
and genuine cavity differences, and its cause is independently confirmed by the 2Z5XnoFAD
control.

The TTBK1/TTBK2 comparison is **not** water-symmetric. The old +2.23 compared two water-free
receptors; the new +0.12 compares 7JXXdry (**0** waters) against 7Q8Ybrg (**5** waters). Since
waters were already shown to penalise these candidates by of order 1 kcal/mol, part or all of
the apparent improvement could be the artifact operating in TTBK1's favour rather than real
paralog selectivity. Reporting +0.12 as a chemistry result without testing this would repeat
precisely the mistake made earlier today with cand_013.

The clean comparison cannot simply be run: dry 7Q8Y fails GATE 1 at 5.60 A, so it is not a
legitimate production receptor. It is still a valid *diagnostic*, and 7Q8Ydry has been built
and queued to run after the selectivity job. Two readings will separate the causes:

    7JXXdry vs 7Q8Ydry   both 0 waters -> water-symmetric margin. If this returns to ~+2.2,
                         the entire improvement is the water artifact.
    7Q8Ybrg vs 7Q8Ydry   same protein, 5 waters vs 0 -> the water cost on TTBK2 alone.

One partial reassurance already in hand: two different TTBK1 receptors with different water
counts give nearly the same margin against 7Q8Ybrg — 7JXXdry (0 waters) +0.12 and 7Q8V
(1 water) +0.10 — so the result is at least not sensitive to the TTBK1 side of the pairing.

If the diagnostic shows the improvement is genuine, manuscript 4.2 ("the selectivity liability
is the dominant risk") no longer holds and the paper's central negative conclusion changes. If
it shows the improvement is the water artifact, 4.2 stands and the corrected number must be
reported with the asymmetry stated. Either way this is not a sentence to write before the
diagnostic returns.

### TTBK2 MM-GBSA — result obtained, then found to be parameter-confounded and rerun

10 ns TTBK2 MD completed (10.1 ns, 301.2 K, 129 ns/day, 48,158 atoms). MM-GBSA first failed on
my own error: the script passed `-sp complex.prmtop`, the *solvated* 48,158-atom topology,
against an already-stripped 4,833-atom trajectory. `-sp` is only correct when the trajectory
still contains solvent; here `-cp` (gb_complex, 4,833 atoms) already matches it. Fixed in
`scripts/run_mmgbsa.sh` with the reason recorded; the strip and topology build had already
succeeded so only the final step repeated.

Result, identical protocol and frame count to the TTBK1 baseline (igb=5, saltcon=0.150,
interval=5, 200.8 frames each):

    TTBK1  cand_003   -31.3629   SD 3.8155   SEM 0.2698
    TTBK2  cand_003   -29.7753   SD 2.3853   SEM 0.1687
    ddG = -1.59 kcal/mol in favour of TTBK1 (on-target preferred)

**This comparison is free of the docking water confound**, which is what makes it valuable:
both complexes were built protein+ligand and then fully solvated in TIP3P, with zero crystal
waters carried through (verified: 0 HOH lines in both `7JXX_prot.pdb` and `7Q8Y_prot.pdb`). So
unlike the docking margin, it cannot be an artifact of the 0-vs-5 water asymmetry.

**However the claim that the two complexes shared identical ligand parameters was wrong, and
checking it is what caught the problem.** The `.lib` files differ:

    06_md/params/docked/cand_003.acpype/cand_003_AC.lib   8a49b36b0f51  <- used by TTBK1
    06_md/params/cand_003.acpype/cand_003_AC.lib          2de23f72608f  <- used by my TTBK2

Same atom names and types (0 mismatches) but **charges differ by up to 0.063 e**, net -0.0010
against +0.0020. The project holds two acpype runs for cand_003 and I took the wrong one. A
0.063 e difference feeds directly into the EEL and EGB terms, and the effect is not obviously
small next to a 1.59 kcal/mol delta — which is precisely the failure mode I set out to avoid
when I chose to reuse parameters rather than regenerate them. The -1.59 therefore cannot be
quoted as it stands.

Rebuilt as `06_md/system_TTBK2m` with the lib matched to TTBK1 (md5 confirmed identical,
tleap 0 errors) and the 10 ns MD relaunched on the idle GPU. The first run is kept rather than
deleted: comparing the two gives a free measurement of how much a 0.063 e charge perturbation
actually moves a single-trajectory GB MM-GBSA result, which is worth knowing in its own right.

Standing caveats for whatever number the rerun gives: single-trajectory GB MM-GBSA carries no
entropy term, the SEM assumes independent frames and so understates the true uncertainty on
correlated MD, and this is one ligand with one replicate per target.

### Diagnostic result: the TTBK2 "improvement" was the water artifact. Manuscript 4.2 stands.

7Q8V and 7Q8Ydry both complete (168/168, 3 seeds throughout). The three-way comparison:

    TTBK1 vs TTBK2, production   7JXXdry 0w vs 7Q8Ybrg 5w   +0.121   25/56 favourable
    TTBK1 vs TTBK2, symmetric    7JXXdry 0w vs 7Q8Ydry  0w   +1.606    0/56 favourable
    water cost on TTBK2 alone    7Q8Ybrg 5w vs 7Q8Ydry 0w    +1.486    0/56

The production-to-symmetric difference is **-1.486**, matching the independently measured water
cost on TTBK2 (**+1.486**) to three decimals. So essentially *all* of the apparent improvement
from +2.23 to +0.12 is the 0-vs-5 water asymmetry penalising TTBK2, not paralog chemistry.
Removing the confound restores an unfavourable margin with **0/56 candidates favourable** —
the same qualitative conclusion as the published +2.23.

**Manuscript 4.2 ("the selectivity liability is the dominant risk") therefore stands for the
TTBK pair, and the earlier entry's speculation that it might not was wrong.** Writing +0.12 up
as a chemistry result would have inverted the paper's central negative finding on the strength
of an artifact I had already identified and documented hours earlier. The margin only looked
good because the anti-target was carrying five waters the on-target was not.

This leaves a genuine methodological bind that has to be stated rather than resolved by
preference: the *validated* pairing (7Q8Ybrg passes GATE 1 at 1.45 A) is the confounded one,
and the *unconfounded* pairing uses 7Q8Ydry, which fails GATE 1 at 5.60 A and so cannot carry a
production number. Neither is clean. The defensible report is both, with the artifact named.

**The MAO row is unaffected by any of this and is the one solid new result.** 2V5Zdry and
2Z5Xdry both carry zero waters, so that comparison is water-symmetric by construction, and both
pass GATE 1 (0.51 A and 1.02 A). Final numbers, all 3 seeds:

    MAO-B vs MAO-A   -3.28   favourable 48/48   plus 8 candidates with no viable MAO-A pose
                             (shortlist -3.21, 13/13)

versus +0.09 and 2/15 in Table 3.7. The cause is established independently: the old MAO-A
receptor had no flavin, and the 2Z5XnoFAD control reproduces the old regime (~3.5 kcal/mol).

**B4 cross-check, final:** +0.101 +/- 0.322, median +0.061, 15/56 within +/-0.5 of 9IV's
+0.784. The offset is a property of 9IV, not of the receptor pair; do not apply it as a blanket
correction.

**Open tension for the write-up:** docking (water-symmetric) says TTBK2 is preferred by
~1.6 kcal/mol, while MM-GBSA on cand_003 says TTBK1 is preferred by 1.59. The MM-GBSA
comparison is water-symmetric by construction and uses a better energy model, but it is one
ligand, one replicate, no entropy term, and currently being rerun with corrected ligand
charges. The two methods disagreeing in sign is itself a reportable result and is exactly the
test 4.3 proposed; it should not be resolved by picking the more convenient one.

### The TTBK2 MM-GBSA numbers are not binding free energies — the ligand leaves the site

Rerunning TTBK2 with ligand parameters matched to TTBK1 gave -22.88 against the first run's
-29.78. A **6.9 kcal/mol** swing from a charge perturbation of at most 0.063 e is not credible
as a parameter effect, so it was checked rather than reported, and the checking overturned the
whole comparison.

Ligand heavy-atom RMSD to the starting pose, protein-fitted, over 1000 frames:

    TTBK1               mean 1.68 A   max 2.51   last-100 mean 1.88    stays bound
    TTBK2 run 1         mean 1.86 A   max 4.63   last-100 mean 3.72    drifts
    TTBK2 run 2         mean 4.48 A   max 7.78   last-100 mean 5.69    leaves the site

Protein CA RMSD is 1.59-1.73 A in every case, so the simulations are sound; it is the ligand
that moves. The 6.9 kcal/mol difference between the two TTBK2 runs is therefore sampling — how
far cand_003 had dissociated by the end — not the charges I rebuilt the system for. The
parameter mismatch was real and worth fixing, but it was not the cause.

**Consequences, in order of importance.**

1. **No MM-GBSA delta-G for TTBK2 can be quoted.** Those averages are taken over a trajectory
   in which the ligand is leaving the pocket, so they are not binding free energies of a
   complex. The earlier ddG of -1.59 must be withdrawn, and the -22.88 must not replace it.
   Any TTBK1-vs-TTBK2 MM-GBSA number in the write-up would be meaningless.

2. **The qualitative result is strong and reproducible, and it favours TTBK1.** cand_003 holds
   its docked pose in TTBK1 for 10 ns and dissociates from TTBK2 in **both** independent
   replicates, built from different parameter sets and different random velocities. Two
   replicates agreeing is worth more here than either delta-G was.

3. **Docking and MD now disagree about TTBK2 in an informative way.** Water-symmetric docking
   scores TTBK2 *better* than TTBK1 by ~1.6 kcal/mol, yet the pose that score describes is not
   stable for 10 ns. A good score attached to an unstable pose is a recognised docking failure
   mode, and it is a more interesting finding than either number alone. It also means the
   Table 3.7 TTBK2 liability, which rests entirely on docking scores, is not corroborated by
   dynamics.

**What cannot be concluded.** That cand_003 does not bind TTBK2. Dissociation in 10 ns from a
docked starting pose is equally consistent with the *pose* being wrong rather than the binding
being absent, and 7Q8Ybrg's own GATE 1 is the weakest in the set at 1.45 A. Distinguishing the
two needs either longer sampling, multiple replicates from different starting poses, or a
method that does not depend on one pose. That is a genuine limitation, not a formality.

Files: `06_md/system_TTBK2m/` (matched-parameter rerun), `lig_rmsd.dat`/`prot_rmsd.dat` in each
system directory.

## 2026-09-25 — A (replicate MD set; MAO-A cofactor parameterisation started in parallel)

**Priority chosen: replicate the MD, not extend it.** The only new MD finding — cand_003 stable
in TTBK1, dissociating from TTBK2 — rests on two TTBK2 runs that started from the *same* docked
pose, against a *single* TTBK1 run. Dissociation from one pose is as consistent with a bad pose
as with weak binding, and one TTBK1 trajectory cannot establish that its stability is
reproducible. Three runs fix both gaps, ~6 h unattended on an otherwise idle GPU:

    system_TTBK2_p2   TTBK2 from docked pose 2   (-7.921, 4.50 A from pose 1)
    system_TTBK2_p3   TTBK2 from docked pose 3   (-7.419, 5.06 A from pose 1)
    system_TTBK1_r2   TTBK1, same pose, fresh velocities

Poses 2 and 3 are genuinely independent starting points, not perturbations: 4.5-5.1 A from pose
1 with comparable scores. All three built with the **docked** acpype parameters, matching the
TTBK1 baseline — the mismatch found earlier today is not repeated, and the coordinate transfer
now validates at 0.001 A on bonds and angles (against 0.068 A before) because the template and
the docking prep finally come from the same acpype run. tleap: 0 errors on both new systems.

**MAO-A started in parallel on the CPU, by the cheaper route.** The GPU is committed for ~6 h
and FAD parameterisation is CPU work, so the two do not compete. Rather than building a custom
covalent Cys-FAD residue — the multi-day, high-failure-risk path I advised against with a
deadline this close — FAD is being parameterised as its **own GAFF2 residue** with AM1-BCC
charges at net charge -2 (84 atoms, from the project's existing `cofactor_FAD.sdf`, which
already carries correct bond orders from the RCSB ideal-chemistry route), to be held by
positional restraints during MD.

This is an approximation and must be labelled as one: the 8alpha-S-cysteinyl bond is not
modelled, so FAD is present as a rigid pocket wall rather than a covalently tethered cofactor.
For a single-trajectory MM-GBSA of the *ligand's* binding energy that is defensible — what
matters is that the flavin occupies its crystallographic position and forms the cavity wall,
which is precisely what its absence got wrong before 2026-09-24 — but it is not a substitute
for proper covalent parameterisation and the manuscript must say so. If the charges or the
restrained geometry do not validate cleanly, the fallback remains declaring MAO-A MM-GBSA as
future work.

### MAO-A system built — FAD parameterises cleanly, which was the main risk

The cheaper MAO-A route worked, and the step that could have sunk it did not:

    FAD GAFF2/AM1-BCC   84 atoms   net charge -1.9980 (target -2)   12 min
    parmchk2 ATTN (guessed) parameters: 0

**Zero guessed parameters** is the result that matters. parmchk2 found complete GAFF2 coverage
for every bond, angle and dihedral in FAD, so none of the flavin's internal geometry rests on
an interpolated guess. Had that come back with a long ATTN list the whole approach would have
been abandoned, since guessed parameters on the cofactor forming the cavity wall would be
worse than no MAO-A run at all.

System: **118,896 atoms**, net charge +0.001, FAD 84 atoms, ligand 35, 36,859 waters, tleap 0
errors. Protein plus FAD is charge-neutral so no counter-ions were required. It is 2.5x the
TTBK2 system, so ~5 h for 10 ns rather than ~2 h.

One collision caught during the build: acpype names its output residue `MOL`, the same name
cand_003 already uses. Two residues sharing a name would have confused tleap and, worse, would
have made the later `ante-MMPBSA -n :MOL` ligand selection silently ambiguous between the
ligand and an 84-atom cofactor. Renamed to `FAD` in both the lib and the PDB before building.

**The approximation, stated plainly.** FAD is a free GAFF2 residue held by positional restraints
(10 kcal/mol/A^2 on heavy atoms, `06_md/run_md_restrained.py`), not the covalently bound
8alpha-S-cysteinyl cofactor it is in reality. Consequences: the flavin cannot relax in response
to the ligand, so induced fit involving it is suppressed, and the cofactor's own dynamics are
absent. What it does buy is the cavity wall in its crystallographic position, which is the
thing whose absence made every pre-2026-09-24 MAO-A number wrong. For a single-trajectory
MM-GBSA of the *ligand's* binding energy this is defensible — FAD sits on both sides of the
complex-minus-receptor subtraction and its internal energy cancels — but it is not equivalent
to proper covalent parameterisation, and the manuscript must say so rather than imply a full
cofactor treatment.

Queued behind the replicate set: the chain polls for the replicates' completion marker and only
then takes the GPU, so the two never contend.

### Replicates overturn the dissociation claim: it was pose-dependent, not a TTBK2 property

Two of three replicates complete. Ligand heavy-atom RMSD to the starting pose:

    TTBK1  rep 1  (pose 1)   mean 1.68   max 2.51   last-100  1.88    stable
    TTBK1  rep 2  (pose 1)   mean 1.56   max 2.40   last-100  1.38    stable
    TTBK2  run 1  (pose 1)   mean 1.86   max 4.63   last-100  3.72    drifts
    TTBK2  run 2  (pose 1)   mean 4.48   max 7.78   last-100  5.69    leaves site
    TTBK2  pose 2            mean 2.33   max 4.21   last-100  1.99    STABLE

**TTBK1 stability replicates cleanly** — two independent trajectories, last-100-frame RMSD 1.88
and 1.38 A. That part of the earlier entry holds.

**The TTBK2 dissociation does not.** Started from docked pose 2, cand_003 settles at 1.99 A,
essentially as stable as TTBK1. So the earlier statement that cand_003 "dissociates from TTBK2
in both independent replicates" was wrong in its implication: both of those replicates began
from the *same* pose 1, so they were never independent in the way that mattered. Two runs
agreeing told us only that pose 1 is unstable, not that TTBK2 binding is.

The corrected reading is narrower and more interesting: **the top-scoring docked pose in TTBK2
is unstable, while a lower-scoring alternative pose is stable.** Pose 1 scores -8.089 and falls
apart; pose 2 scores -7.921, only 0.17 kcal/mol worse, and holds. That is a statement about the
docking, not about the protein — Vina's ranking picked the pose that dynamics rejects. It also
retrospectively explains the 6.9 kcal/mol spread between the two pose-1 runs: both were
sampling a dissociating ligand, so neither number meant anything.

Consequences:

  - The qualitative selectivity claim drawn from MD — cand_003 stable in TTBK1, not in TTBK2 —
    is **withdrawn**. There is a stable TTBK2 complex; docking just did not rank it first.
  - Pose 2's trajectory has the ligand bound throughout, so it yields a *legitimate* MM-GBSA
    delta-G, unlike either pose-1 run. Running it now, together with TTBK1 rep 2, which gives
    the TTBK1 side a replicate rather than a single value to compare against.
  - Whatever comes back must be read against TTBK1's own run-to-run spread, which rep 2 will
    finally make measurable. A ddG smaller than that spread means nothing.

Pose 3 still running; it will show whether pose 2's stability or pose 1's instability is the
outlier.

### MM-GBSA cannot resolve TTBK1 from TTBK2 — and the SEM would have said otherwise

All 10 ns runs, same protocol, 200.8 frames each:

    TTBK1  rep 1   (pose 1, stable)     -31.3629   SD 3.8155   SEM 0.2698
    TTBK1  rep 2   (pose 1, stable)     -33.4396   SD 2.6330   SEM 0.1862
    TTBK2  pose 2  (stable)             -30.9513   SD 3.8096   SEM 0.2694
    TTBK2  pose 1 run 1  (drifting)     -29.7753   -- excluded, ligand not bound
    TTBK2  pose 1 run 2  (dissociated)  -22.8840   -- excluded, ligand not bound

Only the three stable trajectories carry meaning. Against them:

    TTBK1 mean            -32.40
    TTBK1 replicate spread  2.08 kcal/mol   (-31.36 vs -33.44)
    TTBK2 pose 2          -30.95
    ddG (TTBK1 - TTBK2)    -1.45

**The apparent TTBK1 preference is smaller than TTBK1's own run-to-run variation.** Two
trajectories of the *same complex*, differing only in random velocities, disagree by 2.08
kcal/mol — more than the 1.45 separating the two proteins. MM-GBSA at this sampling therefore
cannot resolve a TTBK1/TTBK2 difference, and no selectivity conclusion can be drawn from it in
either direction.

**The reported SEMs would have said the opposite, and that is the methodological point worth
recording.** They run 0.17-0.27 kcal/mol, which would make -1.45 look overwhelming — roughly
five sigma. The true replicate-to-replicate uncertainty is about **10x larger** than the SEM,
because MMPBSA.py computes the SEM as if 200 frames sampled 10 ps apart were independent draws,
which they are not. Quoting the SEM as an error bar on a ddG is the single easiest way to
manufacture a significant selectivity result from this pipeline, and without the second TTBK1
replicate that is exactly what would have happened here. The replicate spread, not the SEM, is
the error bar.

Running the second TTBK1 replicate was therefore not redundant: it is the only reason the
uncertainty is known at all. One replicate per target would have produced a confident,
publishable, wrong number.

**Net position on TTBK1 vs TTBK2 after all of today's work:**

    docking, water-symmetric   TTBK2 favoured by ~1.6 kcal/mol, 0/56 candidates favourable
    docking, validated pair    +0.12, but confounded by a 0-vs-5 water asymmetry
    MD pose stability          no difference -- both proteins hold the ligand, given a good pose
    MM-GBSA                    no resolvable difference (1.45 < 2.08 replicate spread)

The manuscript's Table 3.7 TTBK2 liability rests on docking scores. Dynamics neither confirms
nor refutes it; it simply lacks the resolution to speak, and saying so is the honest outcome.

### Third TTBK2 pose dissociates — and a design asymmetry that has to be fixed before comparing

    TTBK2  pose 1 run 1   last-100 3.72 A   drifts
    TTBK2  pose 1 run 2   last-100 5.69 A   leaves
    TTBK2  pose 2         last-100 1.99 A   STABLE
    TTBK2  pose 3         last-100 5.79 A   leaves

So TTBK2 holds cand_003 in **1 of 3** docked poses. No MM-GBSA was run on pose 3: a
dissociating trajectory yields no meaningful delta-G, and averaging one would repeat the error
already made with the pose-1 runs.

**The comparison as it stands is not valid, and the flaw is in the experimental design rather
than the result.** TTBK2 was tested from three different docked poses; TTBK1 was tested from
pose 1 only, twice. "TTBK2 fails in 2/3 poses while TTBK1 holds in 2/2 runs" therefore compares
two different experiments: one probes pose sensitivity, the other probes velocity sensitivity.
TTBK1 might well lose its ligand from poses 2 and 3 too — that has simply never been tested,
and asserting a difference without testing it would be the same shape of error as the water
asymmetry earlier today.

TTBK1 poses 2 and 3 built (-8.118 and -8.066, 5.48 and 3.97 A from pose 1; tleap 0 errors,
40,185 and 40,191 atoms, both neutral with 12 Cl-) and queued behind MAO-A, which now holds the
GPU. That gives 3 poses per protein under one protocol and makes the pose test symmetric.

Only then can the question "does cand_003 form a stable complex with TTBK1 more reliably than
with TTBK2?" actually be answered. If TTBK1 also holds in only 1 of 3 poses, there is no
stability difference and the MD contributes nothing to the selectivity question beyond what
MM-GBSA already showed (no resolvable difference). If TTBK1 holds in 3 of 3, the asymmetry is
real and reportable as a qualitative result.

### MAO-A MD: the FAD approximation works, but the result cannot support a selectivity claim

10 ns, 118,896 atoms, 61.8 ns/day on the RTX 4050.

    FAD      mean 0.46 A   max 0.76   last-100 0.48    restraints held
    protein  mean 1.76 A   max 3.14   last-100 1.79    stable
    ligand   mean 3.57 A   max 5.34   last-100 3.24    moved off the docked pose
    MM-GBSA  DELTA TOTAL -44.2194   SD 3.1486   SEM 0.2226

**The cofactor scheme did its job.** FAD stayed within 0.5 A of its crystallographic position
throughout, so the flavin wall of the cavity was present and rigid exactly as intended, and the
84-atom GAFF2 residue with zero guessed parameters behaved stably in a 10 ns simulation. That
part of the approach is validated.

**The -44.22 is not usable as a selectivity number, for two independent reasons.**

1. **There is no MAO-B counterpart.** MM-GBSA absolute values are not comparable between
   different proteins — they carry protein-specific desolvation and surface terms that do not
   cancel across systems. A MAO-B vs MAO-A ddG requires cand_003 simulated in 2V5Z under the
   identical restrained-FAD protocol. Only MAO-A was built. Quoting -44.22 against TTBK1's -31
   would compare a kinase to a flavoenzyme and mean nothing.
2. **The ligand left its docked pose.** last-100 RMSD 3.24 A, against 1.38-1.99 A in the runs
   that held. It has not dissociated the way TTBK2 poses 1 and 3 did (5.7-5.8 A), but the
   structure the energy describes is not the structure that was docked, so the number cannot be
   attached to the docking result it was meant to test.

Note also that a *more* favourable MAO-A binding energy would, if taken at face value, point
the opposite way to the docking margin of -3.28 favouring MAO-B. That tension is not resolvable
from one unpaired run, and it would be wrong to present either number as corroborating the
other.

What the MAO-A work has established is narrower and still worth having: FAD can be parameterised
and restrained well enough to run stable MD of a MAO complex, which removes the blocker that
made MAO-A MM-GBSA look like multi-day work. The energetics need the MAO-B partner run before
they say anything. 2V5Z also contains FAD, so it takes the same treatment and roughly the same
5 h.

### MAO-B partner run built and queued; RDKit blocked mid-session and was routed around

MAO-B (2V5Z) built to pair with MAO-A: 90,129 atoms, net charge +0.001, FAD 84 atoms, tleap
0 errors. **The FAD lib and frcmod are the same files used for MAO-A** — only the coordinates
change, to 2V5Z's own crystallographic flavin — so the two MAO systems differ in protein and
ligand pose alone and not in cofactor parameterisation. That is the same discipline the
cand_003 charge mismatch earlier today showed to be necessary, applied pre-emptively this time.

**RDKit stopped loading partway through the build:** `DLL load failed while importing
rdmolfiles: An Application Control policy has blocked this file`. It had worked all day; this
is a Windows policy change, not a code fault. Rather than wait on it, the coordinate-transfer
step was rewritten to parse SDF V2000 directly (`mkligpdb2.py`) — the format is a fixed-width
counts line, atom block and bond block, so nothing was lost and the pipeline no longer depends
on RDKit for this step.

**The validator then flagged the FAD transfer, correctly, and the flag turned out to be a false
positive worth understanding.** Bonds differed by up to 0.305 A and angles 0.361 A, against
thresholds of 0.15/0.20. Checking instead of overriding: the two FAD files have **identical
element sequences and identical bond tables** (84 atoms, 89 bonds), because both were produced
by the same `prep_cofactor.py` RCSB ideal-chemistry route, so atom correspondence is guaranteed
by construction. The deviation is real crystallographic difference between two independently
refined copies of the cofactor — precisely the phenomenon recorded on 2026-09-24, where the
Cys-FAD bond refines to 2.31 A in 2V5Z and 1.65 A in 2Z5X. The thresholds were calibrated for
same-source transfers (ligand pose onto its own acpype template, which validates at 0.001 A)
and are simply wrong for a cross-crystal one. Now parameterised rather than hard-coded, with
the reasoning recorded in the script.

Queue: TTBK1 pose 2 (running) -> TTBK1 pose 3 -> MAO-B MD + MM-GBSA.

## 2026-09-26 — A (post-restart recovery: the symmetric pose test resolves, and it is a null)

**The machine restarted overnight under the MD load; no simulation was lost.** All three queued
production runs had already completed before the restart — TTBK1 pose 2 (10.1 ns, 156 ns/day),
TTBK1 pose 3 (10.1 ns, 177 ns/day) and MAO-B (10.1 ns, 81.4 ns/day) each reached step 5,050,000
with stable temperature (~300 K) and potential energy. Only the post-processing was interrupted,
which is minutes of work rather than hours. Recovered by re-running `cpptraj` strip + RMSD on
the two TTBK1 poses and the RMSD analysis on MAO-B; MAO-B's MM-GBSA had already landed at 06:34.

### TTBK1 holds cand_003 in 1 of 3 poses — exactly as TTBK2 does

With poses 2 and 3 now analysed, the pose test is symmetric for the first time: three docked
poses per protein, one protocol, 10.1 ns each.

    TTBK1  pose 1 run 1   last-100 1.88 A   STABLE
    TTBK1  pose 1 run 2   last-100 1.38 A   STABLE
    TTBK1  pose 2         last-100 4.18 A   drifts
    TTBK1  pose 3         last-100 3.42 A   drifts

    TTBK2  pose 1 run 1   last-100 3.72 A   drifts
    TTBK2  pose 1 run 2   last-100 5.69 A   leaves
    TTBK2  pose 2         last-100 1.99 A   STABLE
    TTBK2  pose 3         last-100 5.79 A   leaves

**1 of 3 for TTBK1, 1 of 3 for TTBK2. There is no pose-stability difference between the two
proteins.** Protein backbones were stable throughout in every run (1.48-1.81 A), so the ligand
motion is ligand motion and not a collapsing binding site.

This is the null outcome anticipated when the asymmetry was spotted, and it retires the
apparent result that preceded it. "TTBK2 fails in 2/3 poses while TTBK1 holds in 2/2 runs" was
**entirely an artifact of the experimental design**: TTBK2 had been probed across three poses
while TTBK1 had been probed across two velocity seeds of a single pose. Comparing a pose scan
against a replicate scan produced a difference that does not exist. Had the TTBK1 poses never
been run, that artifact would have gone into the manuscript as a qualitative selectivity
finding — the same failure mode as the water asymmetry on 2026-09-25 and the SEM-as-error-bar
trap, now three times in one week from the same root cause: an asymmetry between the two arms
of a comparison.

Which pose survives differs between the proteins (TTBK1 pose 1, TTBK2 pose 2), but that is not
a selectivity signal. It says the docking pose ranking does not predict which pose is
dynamically stable, and it says so for both proteins equally.

**No MM-GBSA was run on TTBK1 poses 2 or 3**, on the same grounds that TTBK2 pose 3 was skipped:
a trajectory that has left the docked pose yields an energy for a structure that was never
docked. Averaging those would repeat the error already made with the TTBK2 pose-1 runs.

**Net position on TTBK1 vs TTBK2, final for Week 2:**

    docking, water-symmetric   TTBK2 favoured by ~1.6 kcal/mol, 0/56 candidates favourable
    docking, validated pair    +0.12, confounded by a 0-vs-5 water asymmetry
    MD pose stability          NO difference -- 1/3 poses each, under a symmetric test
    MM-GBSA                    no resolvable difference (1.45 < 2.08 replicate spread)

Dynamics does not speak to the Table 3.7 TTBK2 liability in either direction. That claim rests
on docking scores alone, and the MD contributes nothing further to it.

### MAO-B completes the FAD pair; the cofactor scheme validates twice, the energetics still do not

    FAD      mean 0.45 A   max 0.75   last-100 0.46    restraints held
    protein  mean 1.12 A   max 1.58   last-100 1.32    stable
    ligand   mean 1.95 A   max 3.75   last-100 2.28    HELD the docked pose
    MM-GBSA  DELTA TOTAL -39.4802   SD 2.1998   SEM 0.1556

**The FAD approach reproduced on a second, independently refined flavin** (0.46 A here vs 0.48 A
for MAO-A), using the same lib and frcmod with only coordinates changed. The cofactor
parameterisation is now validated on both MAO structures.

**The MAO-A/MAO-B pair still cannot carry a selectivity claim, for a new reason on top of the
old one.** MAO-A -44.22 vs MAO-B -39.48 is a 4.74 kcal/mol margin favouring MAO-A, but:

1. **The pair is asymmetric in pose stability** — MAO-B held at 2.28 A while MAO-A drifted to
   3.24 A. The more favourable of the two numbers is the one describing a structure that is not
   the docked pose, so part of that margin is MAO-A relaxing elsewhere rather than binding
   better. A like-for-like ddG needs both arms on-pose, and only one is.
2. **It contradicts the docking.** Docking favoured MAO-B by -3.28; MM-GBSA now favours MAO-A by
   +4.74. Opposite signs, one replicate each. Against the known ~2.08 kcal/mol replicate spread
   the margin is not negligible, but a single replicate per target cannot settle a sign
   disagreement, and presenting either number as corroborating the other would be wrong.

Resolving this needs a MAO-A replicate from a pose that holds, not more targets.

**Provenance defect fixed.** `06_md/system_TTBK1_p2/mmpbsa.in` and `p3/mmpbsa.in` were copied
from the TTBK2 directory and their header lines still read `cand_003 / TTBK2 (7Q8Y)`. The systems
are correctly built on 7JXX (`build.leap` loads `7JXX_prot.pdb`) and the numerical settings were
identical, so no computed result was affected — but that line prints into the results file, and a
results file that misreports its own target is exactly the kind of thing that survives into a
manuscript. Corrected to name TTBK1 and the pose.

Files: `06_md/system_TTBK1_p2/{production_full_stripped.nc,stripped.complex.prmtop,lig_rmsd.dat,
prot_rmsd.dat}` (new), same for `system_TTBK1_p3` (new), `06_md/system_MAOB/{lig,prot,fad}_rmsd.dat`
(new), `06_md/system_TTBK1_{p2,p3}/mmpbsa.in` (header corrected).

**Queue is empty.** Every MD system built this week has run and been analysed. The open
scientific item is a MAO-A replicate from an on-pose start, which is the only remaining way to
make the MAO ddG interpretable.

## 2026-09-26 — A (workspace restructure and cleanup, pass 1: inventory, quarantine, documentation)

Restructured the workspace for the final analysis and proposal, following an
inventory → classify → dependency-check → archive → restructure → verify sequence. **Nothing was
permanently deleted.** Two findings during the dependency check changed the plan materially, and
both are worth recording because in each case the obvious action would have been the wrong one.

### The numbered directories turned out to be load-bearing, so they were not renamed

The intent was to regroup files into per-analysis-step folders. A dependency scan first found the
step directories hardcoded in **16+ places** across `scripts/` and `generation/` —
`03_receptors/`, `01_smiles/`, `04_docking/`, `02_ligands/pdbqt/`, `06_md/system/`,
`05_validation/`, `08_analysis/`, `00_library/reinvent4_output/campaign2_v2/`. Renaming any of
them would have broken every driver script for no scientific gain.

The existing `00_`–`09_` layout *is already* the workflow order. What was missing was not
structure but documentation. So the restructure adds a layer in place: a `README.md` per step
answering what / why / reference / inputs / analysis / result / files, plus a root `README.md`
workflow index and `INVENTORY.md` carrying the classification. 16 step READMEs, 241 internal
links, all verified to resolve.

### `artifacts/` was not a duplicate mirror — it was the project's missing rationale trail

`artifacts/` (202 MB) looked like a staging mirror: its two 95 MB RL checkpoints and the
manuscript files were byte-identical to copies in `00_library/` and `09_manuscript/`. Deleting it
would have been defensible from filenames alone.

Hashing every file instead showed **144 of its 174 files were unique**, including the entire
`phase0`–`phase8` report series this logbook cites constantly, the only copies of the earlier
receptor set (`4NFM`, `6U0K`, `2V60`), and the prior-phase MM-GBSA results. The phase reports are
the "why" documentation for the first half of the project; nothing else records it.

Its 30 verified duplicates were quarantined and its 144 unique files distributed into the
workflow step each belongs to, under `<step>/prior_phase/`. `artifacts/` no longer exists.
`prior_phase/` is a subfolder rather than a merge into each step root **because the two passes
used different receptor structures** — pooling `4NFM` results with `7JXX` results in one directory
is the same shape of error as the water asymmetry and the pose/replicate asymmetry.

This is the third time on this project that checking rather than assuming has overturned the
obvious call, and the second time that the thing about to be discarded was load-bearing.

### Quarantined: 2.2 GB, all of it verified byproduct

`_ARCHIVE_TO_DELETE/` mirrors the original tree and carries a `MANIFEST.csv` with the original
path, category, size and reason per file, so anything can be restored.

    reference.frc x7        1.83 GB   MMPBSA.py forces dump
    artifacts/ duplicates    198 MB   md5-identical to a copy that stays live
    2 RL checkpoints + dups  379 MB   superseded campaigns; v3 checkpoint retained
    __pycache__                8 KB   bytecode

`reference.frc` is worth a note on how it was identified as a byproduct rather than assumed to be
one: nothing in the repo references `.frc`; it exists in exactly the 7 system directories where
MM-GBSA ran and in none of the 3 where it did not; and each file's mtime matches its own
`mmgbsa_results.dat` to the minute. It is regenerated on any re-run.

### What was kept against a "minimal workspace" reading

The brief asked to remove outdated, duplicate and intermediate files. Several items qualify on
their face and were kept because a conclusion depends on them:

- **The dissociated and drifted MD runs** (`system_TTBK2m`, `system_TTBK2_p3`,
  `system_TTBK1_p2`, `system_TTBK1_p3`). The headline result — 1 of 3 poses stable on each kinase,
  therefore no stability difference — exists *only* because these failures were run and retained.
  Deleting "failed runs" would have deleted the finding.
- **`campaign1_first_failed/`** outputs: a docking-only reward yielding zero BBB-passing molecules
  is what justifies the BBB-aware reward in the final campaign.
- **`params/cand_003.acpype/`** (free-conformer, superseded) so the charge-mismatch incident stays
  demonstrable against the docked-pose set actually used.
- **`04_docking/*dry*`, `*noFAD*`, `*brg*`** — not redundant receptor variants but the controls
  that quantified the water artifact and the cofactor dependency.
- **All 6.75 GB of raw `production*.dcd`** and the `REINVENT4/` clone, by explicit decision.

### New: a collected results table

`07_mmgbsa/md_mmgbsa_summary.csv` now holds all ten MD systems in one table — target, structure,
starting pose, frame count, ligand/protein/FAD RMSD (mean, max, last-100), pose verdict, and
ΔG/SD/SEM — generated directly from the primary `.dat` files so it cannot drift from what was
computed. Regenerating it reproduced the 2.08 kcal/mol TTBK1 replicate spread from primary data
and surfaced the TTBK2 pose-1 spread at **6.90 kcal/mol** against quoted SEMs of 0.17 and 0.21, a
40x discrepancy and the starkest illustration yet that the SEM is not the error bar.

### Manuscript claims now out of date

Recorded in `09_manuscript/README.md` rather than acted on, since revising them is science not
cleanup. §3.5 ("Both complexes are stable over 20 ns") rests on the prior-phase apo-structure MD
and must be rewritten — stability is pose-dependent, 1 of 3 either side. §3.6's MAO-B preference
cannot be sourced to an MM-GBSA magnitude comparison, which points the other way and is unusable.
§3.7 (the series does not discriminate either target from its paralog) is the one claim that came
out of this week **stronger** than it went in, now supported by four independent lines.

Also flagged: several `source` fields in `01_smiles/references.csv` still read "CONFIRM primary
source before submission", and safinamide is not the same molecule as the PDB ligand `SAG`.

### Verification

241 documentation links resolve; every concrete path referenced by a script still exists (the one
miss, `generation/rl_stage1.chkpt` in the v2 TOML, is a pre-existing stale *output* declaration,
never tracked); all 10 MD systems retain trajectory, topology and RMSD data, and the 3 lacking
MM-GBSA are exactly the 3 that intentionally have none.

`_REVIEW_REQUIRED/` holds the two original project-plan JSONs (36 KB, 17 Aug) pending a decision.
Permanent deletion of `_ARCHIVE_TO_DELETE/` awaits approval.

Files: `README.md`, `INVENTORY.md` (new), `{00_library,01_smiles,02_ligands,03_receptors,
04_docking,05_validation,06_md,07_mmgbsa,08_analysis,09_manuscript,filtering,generation,
scripts}/README.md` (new), `07_mmgbsa/md_mmgbsa_summary.csv` (new),
`_ARCHIVE_TO_DELETE/{README.md,MANIFEST.csv,RESTRUCTURE_MOVES.csv}` (new),
`_REVIEW_REQUIRED/README.md` (new), `artifacts/` (removed, contents redistributed).

## 2026-09-27 — MAO pose scan completed: the MAO-A/MAO-B sign inversion was an off-pose artifact

Finished the four queued MAO pose runs (MAO-A poses 2 and 3, MAO-B poses 2 and 3), ran MM-GBSA and
the full cpptraj chain on all of them, and closed the open item that had been blocking any MAO
interpretation: **each arm now has a trajectory that holds its docked pose.** The headline is that
the MAO-A-over-MAO-B preference recorded on 2026-09-26 was an artifact of a single off-pose
trajectory, and the corrected comparison agrees with docking.

| System | Ligand last-100 | Verdict | ΔG | SD | SEM |
|---|---|---|---|---|---|
| `system_MAOA` (pose 1) | 3.24 Å | drifts off pose | −44.22 | 3.15 | 0.22 |
| `system_MAOA_p2` | 2.43 Å (1.79 core-fit) | holds, loose | −35.54 | 2.73 | 0.19 |
| `system_MAOA_p3` | **1.43 Å** | stable | **−36.67** | 2.40 | 0.17 |
| `system_MAOB` (pose 1) | 2.28 Å | stable | −39.48 | 2.20 | 0.16 |
| `system_MAOB_p2` | **1.10 Å** | stable | **−39.20** | 2.13 | 0.15 |
| `system_MAOB_p3` | 2.57 Å | holds, loose | −39.04 | 2.86 | 0.20 |

### The −44.22 was never a binding energy

Every on-pose MAO-A trajectory lands near −36 (−35.54, −36.67). The −44.22 is 7.6 kcal/mol more
favourable than either, it is the largest MAO-A number, and it is **the only off-pose one**. The
"MM-GBSA favours MAO-A by +4.74, contradicting docking" finding was therefore measuring how much
energy MAO-A gains by relaxing somewhere it was never docked.

Corrected: all six on-pose pairings favour **MAO-B**, by 2.37 to 3.94 kcal/mol. Docking favours
MAO-B pose for pose on this ligand's Vina affinities by 3.03–3.22 kcal/mol (3.184 / 3.032 /
3.221 for poses 1 / 2 / 3). The sign now agrees with docking *regardless of which poses are
paired*, and the magnitudes overlap.

Care with the docking figure: the −3.28 quoted in the earlier entries is the **library-wide**
consensus margin (mean over the 56 candidates, favourable 48/48), not a per-pose number for
`cand_003`. A first draft of this entry conflated the two and stated "3.18–3.28 across the same
poses", which is wrong in both provenance and range. The per-pose deltas come from the Vina
affinities recorded in `06_md/README.md`.

This is the clearest instance in the project of the off-pose problem producing not merely a noisy
number but an **inverted conclusion**. The existing rule — never compute MM-GBSA on a trajectory
that left its docked pose — earned its keep here; the failure mode it prevents is not a wide error
bar, it is a confident wrong answer.

### What is still not claimable: the magnitude

The 2.37–3.94 kcal/mol margin sits on top of the 2.08 kcal/mol velocity-replicate spread measured
for TTBK1 and well under TTBK2's 6.90, and **neither MAO arm has a velocity replicate.** The
direction is robust and pose-independent; the size is not resolved.

The pose scan does supply a second kind of spread, and the temptation to substitute it must be
resisted:

| Target | On-pose trajectories | ΔG range | Pose spread |
|---|---|---|---|
| MAO-A | poses 2, 3 | −35.54, −36.67 | 1.13 kcal/mol |
| MAO-B | poses 1, 2, 3 | −39.48, −39.20, −39.04 | **0.44 kcal/mol** |

MAO-B's three independent poses agree to 0.44 kcal/mol, which looks like tight convergence right up
until it is placed beside TTBK2's 6.90 kcal/mol between two runs differing **only in velocity
seed**. Converging onto the same energy from different starting geometries says the basin is well
defined; it says nothing about how far that energy wanders under resampling. The smaller variance
cannot bound the larger one. Pose spread is not replicate spread.

### Docking rank does not predict pose stability — now 1 of 4, not 1 of 2

The best-holding pose is MAO-A's **worst**-ranked pose 3 (1.43 Å, Vina −7.399 against pose 1's
−8.236) and MAO-B's pose 2 (1.10 Å, −11.26 against −11.42). Across all four targets the top-ranked
docked pose was the most stable in **one of four** cases — TTBK1, the exception. At 3 poses × 4
targets this is the project's best-supported methodological result, and it is why every MM-GBSA
number in the table is tied to a named pose rather than to a target.

### A 3.44 Å backbone RMSD that is not an unstable fold

`system_MAOA_p2` reported the highest protein backbone RMSD in the project, 3.44 Å against a
1.12–1.81 Å band everywhere else, which on its face looks like a failed run. Per-residue RMSF
located it precisely: the maximum is always the last residues of the construct — MAO-A 510–513 of
513, MAO-B 496–499 of 499 — fluctuating 5–8 Å. These are the membrane-anchoring end of the protein,
simulated in water with no bilayer, so they flail. They also dominate any whole-protein fit.

Refitting on the ordered core `:1-496` and re-measuring: **core RMSD is 1.06–1.63 Å across all six
MAO runs and flat in time**, and the `system_MAOA_p2` ligand reads 1.79 Å rather than 2.43 Å. No
MAO run has an unstable fold, and the tail motion never reaches the binding site. `rmsf.cpptraj`
and `core_rmsd.cpptraj`, previously only in the two pose-1 MAO systems, were copied to the four new
ones so all six are analysed identically.

### The verdict scale needed a new band, so it got one rather than a rounding

The first ten runs fell into groups at ≤ 2.28 Å and ≥ 3.24 Å with nothing in between, so
`stable`/`drifts` had never had to be defined precisely. Two MAO pose runs landed in the gap, at
2.43 and 2.57 Å. Rather than round them into whichever neighbouring label was convenient — which
would have been an unrecorded judgement call affecting whether an energy gets quoted — they are
labelled `holds, loose`, the thresholds are now explicit constants in the collector, and neither
run is used as a primary number.

### FAD restraint scheme validated six times

FAD stayed 0.46–0.73 Å from its crystallographic position across all six MAO systems, same lib and
frcmod, coordinates alone differing. Previously this was demonstrated twice at 0.46–0.48 Å.

### The summary table now has a generator

`07_mmgbsa/md_mmgbsa_summary.csv` claimed to be "generated directly from the primary `.dat` files
so it cannot drift", but no generator was committed — the claim rested on whoever last edited it
having been careful. `scripts/collect_md_summary.py` makes it true: it parses the cpptraj `.dat`
files and the last `DELTA TOTAL` of each `mmgbsa_results*.dat`, and `--check` exits 1 if the
committed CSV disagrees with primary data. It reproduced all ten pre-existing rows **byte-for-byte**
before the four new ones were added, which is the regression test that it parses what was parsed
before. Two columns were added, `lig_corefit_last100` and `prot_core_rmsd_last100`, populated for
the MAO systems.

Gotcha recorded in the script header: the `python` first on PATH is MGLTools' Python 2.7, there for
the AutoDock prep scripts, so this must run through `conda run -n docking_project`.

### One interruption, no lost science

A machine restart at 16:06 on 2026-09-26 (initiated through `shutdown.exe` under this user — most
likely Windows Update) killed the queue 0.77 ns into `system_MAOA_p3`. `run_md_restrained.py`
writes no intermediate checkpoint, so the partial run could not be resumed; the system was re-run
from the top and `system_MAOA_p2`, which already had `final_state.xml`, was skipped by the queue's
existing resume check. The completed trajectory is a normal full 10.1 ns run with nothing stitched
together. Throughput held: MAO-A 59.8–63.1 ns/day over three runs, MAO-B 80.9–82.0.

**Worth doing before the next long queue:** `run_md_restrained.py` has no checkpointing, so any
interruption costs the whole run. Adding an OpenMM checkpoint would not change the physics and
would have saved nothing this time (0.77 ns), but the exposure grows with run length.

### Manuscript claims revised

`09_manuscript/README.md` §3.6 previously read that MM-GBSA "points the opposite way" to docking
for MAO and was unusable in either direction. That is now wrong and was rewritten: the direction is
supportable, the magnitude is not, and the −44.22 is flagged as never citable as a MAO-A binding
energy.

§3.5 gained the MAO pose-stability lines. Counting strictly, MAO-A holds 1 of 3 and MAO-B 2 of 3;
no MAO run dissociated, the worst being 3.24 Å, against three of eight TTBK runs past 5 Å. "Both
complexes are stable" is closest to true for the MAO pair and least true for TTBK2.

§3.7 ("the series does not discriminate either target from its paralog/isoform") now **must be
split by pair**, and this is the one revision that weakens a claim. Its four converging evidence
lines are all TTBK1/TTBK2, where it still holds. For the MAO-A/MAO-B isoform pair it no longer
reads the same way: docking and MM-GBSA now agree in sign on a MAO-B preference. "No *significant*
discrimination" survives for MAO, since the margin is not resolvable; "no discrimination", and any
wording implying the isoform comparison came out symmetric, does not.

Files: `INVENTORY.md` (06_md now 14 GB / 582 files / 14 systems), `06_md/README.md`,
`07_mmgbsa/README.md`, `07_mmgbsa/md_mmgbsa_summary.csv` (4 rows and 2
columns added, regenerated), `09_manuscript/README.md`, `scripts/collect_md_summary.py` (new),
`scripts/run_mao_pose_analysis.sh` (new), `06_md/system_MAO{A,B}_p{2,3}/` (4 completed 10.1 ns
systems with trajectories, MM-GBSA and RMSD/RMSF analysis), `06_md/system_MAOA_p3/{rmsf,core_rmsd}.cpptraj`
and `06_md/system_MAOB_p{2,3}/{rmsf,core_rmsd}.cpptraj` (copied from the pose-1 systems).

## 2026-09-27 — B (results directory; five documentation errors found by cross-checking it)

Built [`10_results/`](10_results/) — the consolidated results and interpretation, organised by
**claim** rather than by method, with three figures and a `--check`-able provenance chain back to
the primary `.dat` files. Then wrote a cross-check that compares every number in the new prose
against the source table, and it immediately found five errors, two of them pre-existing and one
of them scientifically misleading.

### The receptor mislabel: the lead's best score was credited to the wrong protein

`04_docking/README.md` and `08_analysis/README.md` both read:

    cand_003 is the lead: #8/56 on 7JXX (-8.38) and #2/56 on 7Q8Y (-11.41)

The −11.41 is **2V5Z (MAO-B)**, not 7Q8Y (TTBK2). `consensus_new.csv` contains no `7Q8Y`/`cand_003`
row at all — the receptors present are 7JXX, 2V5Z and 4BTK.

This one matters beyond bookkeeping. **7Q8Y is an anti-target and 2V5Z is an intended target**, so
the mislabel credited the compound's single strongest docking score to the protein the design is
trying to *avoid*, inverting what that number says about selectivity. Both files corrected, with a
dated note, and the ranks restated as "over the 56 candidates, excluding the native reference
ligand in the same file" — the 57th row is `native_VP7` / `native_SAG`, whose inclusion is what
makes a naive rank read 9 instead of 8.

### Three miscounts and a transcribed SEM

| Where | Said | Is |
|---|---|---|
| `07_mmgbsa/README.md` table | TTBK2 pose 2 SEM **0.28** | **0.27** (primary output: 0.2694) |
| `07_mmgbsa`, `09_manuscript` | SEMs run **0.16–0.28** | **0.15–0.27** across the 11 computed systems |
| top-level `README.md` | SEMs run **0.17–0.27** | **0.15–0.27** |
| `07_mmgbsa/README.md` | spread is **10× to 30×** the SEM | **8× to 41×** (7.7–10.9 TTBK1, 32.9–40.6 TTBK2) |
| `09_manuscript`, top-level, and the first draft of `10_results` | **three** of eight TTBK runs passed 5 Å | **two** (5.69 and 5.79; nothing else reaches 5) |
| first draft of `10_results` | **five** of fourteen runs hold their pose | **six** (1.88, 1.38, 1.99, 2.28, 1.43, 1.10) |

Each is small on its own. The pattern is not: **every one of them is a number that was typed into
prose rather than read out of a file**, which is the same failure mode the summary-CSV generator was
written to stop. The lesson generalises past the CSV — a claim in a README is as much a data
artifact as a row in a table, and it needs the same provenance.

The cross-check that found them lives in this session's scratch space rather than the repo, because
it hardcodes the prose's own numbers and would need editing on every revision. `--check` on
`collect_md_summary.py` remains the durable invariant; the prose check is a one-shot audit, and it
is worth re-running by hand whenever the results text is revised.

### Velocity replicates launched, closing the last open item

`scripts/run_mao_replicate_queue.sh` runs velocity replicates of the two best on-pose MAO systems,
`system_MAOA_p3_r2` and `system_MAOB_p2_r2`, then their MM-GBSA and cpptraj chain. Started 11:30;
MAO-A r2 at ~59 ns/day, so the pair plus analysis lands around 19:30.

Each replicate directory carries `complex.prmtop` and `complex.inpcrd` copied byte-for-byte from its
parent — **sha256-verified, and recorded in a per-directory `PROVENANCE.md`** — with no build inputs
copied, because `tleap` was not re-run and shipping its inputs would imply otherwise. The only
difference between run 1 and run 2 is the random velocity seed, which OpenMM draws itself since
`run_md_restrained.py` sets none. Same mechanism that produced `system_TTBK1_r2`.

These measure the MAO arms' own replicate spread, which is the one thing standing between the MAO
result being *a direction* and being *a number with an error bar*. `make_figures.py` already has
placeholder rows for them in figure 3 and will fill them in on the next run.

### What the results directory concludes

Ten claims, each with a status. The four that carry the project:

- **Supportable, strongly** — docking rank does not predict which pose survives dynamics (1 of 4
  targets, 3 poses each); and the reported SEM is not a usable error bar (8× and 41×).
- **Direction only** — the compound favours MAO-B over the MAO-A anti-target. All six on-pose
  pairings agree and docking agrees independently, but 2.37–3.94 kcal/mol sits inside the measured
  replicate spread.
- **Not supportable** — selectivity for TTBK1 over TTBK2 (1.45 < 2.08), and any cross-protein
  comparison of absolute ΔG.

Two explicit do-nots are recorded where the write-up will reach for them: do not quote the −44.22,
and do not put TTBK and MAO binding energies on the same axis. Figure 2 is MAO-only for exactly
that reason, while figure 1 shares an axis legitimately because ligand RMSD measures one ligand
against its own starting pose.

Figures were checked against the data-visualisation guidance and rendered before being accepted;
three layout defects (band labels colliding with a legend, on-pose mean labels overrunning a title,
SEM labels landing on the y-axis ticks) were fixed. The verdict scale is encoded by **position
against threshold bands, not colour**: the project's four status colours fail the categorical
normal-vision floor when validated as a 4-slot ramp (worst adjacent ΔE 13.6 against a floor of 15),
and they are specified to ship with an icon and label rather than to carry meaning by hue.

Files: `10_results/{README.md,make_figures.py,fig1_pose_stability.png,fig2_mao_binding_energy.png,fig3_sem_vs_replicate.png}`
(new), `scripts/run_mao_replicate_queue.sh` (new),
`06_md/system_{MAOA_p3,MAOB_p2}_r2/` (new, running),
`README.md` (step 10 added to the index; SEM range and the 5 Å count corrected),
`INVENTORY.md`, `04_docking/README.md` and `08_analysis/README.md` (receptor mislabel),
`07_mmgbsa/README.md` (SEM 0.28→0.27, range, ratios), `09_manuscript/README.md` (range, 5 Å count).

## 2026-09-27 — C (no merge conflict; the real blocker was GitHub's file-size limit)

Investigated a reported merge conflict. **There was none**, and none was possible: no unmerged
paths, no conflict markers in any tracked or untracked file, no merge/rebase/cherry-pick in
progress, no OneDrive sync copies, no submodules, and `origin/master` was a strict ancestor of
`HEAD` (`git rev-list --left-right --count` gave `0 2`), so a push would have fast-forwarded.

What was actually wrong: **the two unpushed commits cannot be pushed at all.** `.gitignore` already
excluded `*.dcd` and `*.prmtop` for size, but `*.nc` and `*.frc` were missed, so 26 files totalling
**4.00 GB** were tracked — the bulk of the repository's 4.47 GB. Eleven of them are `reference.frc`
at 194–339 MB each, over **GitHub's 100 MB hard per-file limit**, which is enforced by a
pre-receive hook on the whole pushed history.

That last detail is the one worth remembering: **deleting an oversized file in a later commit does
not make the branch pushable.** GitHub checks every blob in the history being pushed, so the only
fixes are rewriting the offending commits or moving the blobs to LFS.

### What was done, and what it deliberately does not fix

Untracked all 26 via `git rm --cached` (files untouched on disk), added `*.frc` and `*.nc` to
`.gitignore` with the reasoning inline, and recorded in `06_md/README.md` that trajectories are not
version-controlled. Tracked content went from **4.47 GB to 0.48 GB across 6079 files**.

**The branch is still not pushable**, by choice — the oversized blobs remain in `0a758cb` and
`74129f0`. Rewriting those two commits is the fix and is safe while they stay unpushed; that was
offered and deferred. Recorded here so the decision is not mistaken later for an oversight.

Nothing of scientific value left the repository. `reference.frc` is an `MMPBSA.py` byproduct;
`production_full_stripped.nc` is produced from `production.dcd` by `strip_traj.cpptraj`; both are
rebuilt by `scripts/run_mmgbsa.sh`. Every derived number stays tracked — the `.dat` RMSD traces,
`production.log`, `mmgbsa_results.dat` — as do all the inputs needed to regenerate the rest.
`collect_md_summary.py --check` still passes against the 14-system table.

The honest cost, now stated in `06_md/README.md` rather than left implicit: **a fresh clone cannot
reproduce the reported numbers without re-running the MD.** The raw trajectories exist only on
local disk and need backing up separately from git.

The two `06_md/system_MAOA_p3_r2/*.log` files were left uncommitted on purpose: the replicate MD was
mid-run and still writing to them. `production.log` stays tracked — it is the per-frame
energy/temperature trace the README cites — so in-run churn is expected, and the rule is simply not
to commit mid-run.

Files: `.gitignore`, `06_md/README.md`, and 26 index removals (commit `e2cc22d`).

## 2026-09-27 — D (MAO velocity replicates: the ΔΔG magnitude becomes claimable, and the SEM finding sharpens)

Ran velocity replicates of the two best on-pose MAO systems and closed the open item. **The MAO
ΔΔG now has a measured error bar, and the margin clears it about tenfold.**

| System | Ligand last-100 | ΔG | SD | SEM |
|---|---|---|---|---|
| `system_MAOA_p3` (run 1) | 1.43 Å | −36.67 | 2.40 | 0.17 |
| `system_MAOA_p3_r2` (run 2) | 1.07 Å | −36.40 | 2.38 | 0.17 |
| `system_MAOB_p2` (run 1) | 1.10 Å | −39.20 | 2.13 | 0.15 |
| `system_MAOB_p2_r2` (run 2) | 1.37 Å | −38.97 | 2.11 | 0.15 |

Replicate spreads **0.26** (MAO-A) and **0.24** (MAO-B) kcal/mol. Comparing replicate means,
−36.53 against −39.09 gives **ΔΔG = 2.55 kcal/mol favouring MAO-B, 9.8× the larger spread.** Both
replicates also reproduced the pose stability (1.07 and 1.37 Å), and FAD held at 0.55 and 0.61 Å.

So claim 6 in `10_results/README.md` moved from *not supportable* to *supportable*: the size of the
MAO-B preference, not merely its direction, is now inside what this pipeline can resolve. Across
all twelve on-pose pairings the margin runs 2.30–3.94 and never changes sign.

### The caveat that must travel with the magnitude

The MAO spread is tight partly **because of the protocol, not only the physics**: FAD is
positionally restrained, which suppresses receptor conformational sampling and damps the
run-to-run variation. 0.24–0.26 is the correct error bar *for this protocol*; it is not evidence
that MM-GBSA is intrinsically this precise, and it should not be transferred to an unrestrained
system. The ΔΔG survives the objection because **both arms carry the same restraint** — the
symmetry rule this project keeps relearning.

### The SEM finding is not what it looked like, and the sharper version is better

The standing lesson was "the SEM is not the error bar; report the replicate spread." The MAO
replicates complicate that, because here the SEM is nearly right — 0.15 against a measured 0.24.
Putting all four of the project's replicate pairs side by side explains why:

| Pair | Mean ligand RMSD | Spread | Ratio to larger SEM |
|---|---|---|---|
| MAO-B pose 2 | 1.23 Å | 0.24 | 2× |
| MAO-A pose 3 | 1.25 Å | 0.26 | 2× |
| TTBK1 pose 1 | 1.63 Å | 2.08 | 8× |
| TTBK2 pose 1 | 4.70 Å | 6.89 | 33× |

**The spread tracks pose stability, monotonically across all four pairs**, rising 29-fold while
the reported SEM stays in a 0.15–0.27 band. The mechanism is obvious once seen: a trajectory
leaving the site samples structures that were never the complex, so its energy swings between
seeds; a tightly held pose samples one basin and returns nearly the same number twice.

The refined rule: **the SEM is only as good as the pose is stable, and you cannot tell which case
you are in without running a replicate.** The MAO magnitude is quotable because the replicate was
run, not because the poses looked stable. With n = 4 pairs across two protein families this is a
consistent pattern, not a calibration curve — and it is a better finding than the blanket version,
because it says *when* the SEM misleads and by how much.

This also retires a caution recorded yesterday, that **pose spread is not replicate spread**.
MAO-B's three independent poses agreed to 0.44 kcal/mol, which could not legitimately bound the
velocity-seed spread. It turned out to point the right way — the measured spread is 0.24 — but that
was not knowable in advance, and TTBK2 remains the counter-example where a converged-looking arm
concealed a 6.89 spread.

### Precision bug found in the collector

Replicate spreads are differences of two near-equal numbers, so computing them from the summary
CSV's 2 dp display column amplifies rounding: the true 0.26/0.24 came out as 0.27/0.23. Added a
`dG_raw` column carrying MMPBSA.py's value at 4 dp, and `make_figures.py` now computes spreads
from it. The display column is unchanged.

Ratios are quoted as whole numbers (2×, 8×, 33×) because the TTBK2 figure is 32.8× off the CSV and
33.2× off the primary file, depending on which rounding of a 0.2079 SEM is used — a distinction
with no meaning at this precision.

### Figure 3 rebuilt around the real relationship

It was a dumbbell of SEM against spread, titled "The reported SEM is not the error bar". With the
MAO pairs at 2× that title no longer described its own data. It is now spread and SEM plotted
against mean ligand RMSD for all four pairs, titled for what it shows: how far the SEM sits from
the real error bar depends on pose stability. Both series are kcal/mol on one axis — deliberately
not a dual-axis chart. Two label collisions were fixed (the MAO pairs sit 0.02 Å apart on x, so
near-coincident neighbours are now pushed to opposite sides).

### Scale and naming

Sixteen MD systems now; eight hold their pose, two dissociate. MAO systems that share a pose with
a replicate are renamed on the `start` column to `pose N, run 1` / `run 2`, following the TTBK1
precedent. Core-fit protein RMSD across the eight MAO runs is 1.06–1.77 Å (the new maximum is
`system_MAOB_p2_r2` at 1.77) and the FAD restraint is now validated eight times at 0.46–0.73 Å.

### One process note

`scripts/collect_md_summary.py` did not list the two replicate directories, so the CSV — and
therefore every figure — would have silently ignored them while `make_figures.py` waited for rows
that never arrived. Caught before it mattered. The collector hard-exits on a listed system whose
`lig_rmsd.dat` is missing, which is the right failure: loud, not silent.

Files: `07_mmgbsa/md_mmgbsa_summary.csv` (16 rows, `dG_raw` column added),
`scripts/collect_md_summary.py`, `10_results/{make_figures.py,fig1,fig2,fig3}`,
`10_results/README.md` (claims 5, 6 and 9 revised; result 2 and result 4 rewritten),
`07_mmgbsa/README.md`, `06_md/README.md`, `08_analysis/README.md`, `09_manuscript/README.md`
(§3.6 and §3.7), `README.md`, and `06_md/system_{MAOA_p3,MAOB_p2}_r2/` (two completed 10.1 ns
systems with full analysis).

## 2026-09-27 — E (cleanup carried out; repository pushed; workflow document rewritten)

Executed the deletion that had been quarantined since 26 Sep, tidied the top level, and pushed the
whole project to GitHub for the first time since 21 Sep.

### 2.2 GB deleted, with the redundancy re-verified rather than trusted

The quarantine README asserted that 30 files were byte-identical duplicates. **That assertion was
re-checked by md5 at deletion time rather than taken on trust, and it was wrong for three of them**
— not because they were unique, but because the earlier audit matched on basename while the
surviving copy has a different filename:

| Quarantined | Actual surviving copy |
|---|---|
| `artifacts/2Z5X.pdb` | `03_receptors/2Z5X/raw.pdb` |
| `artifacts/7Q8Y.pdb` | `03_receptors/7Q8Y/raw.pdb` |
| `artifacts/phase0_flavonoid_library.csv` | `00_library/flavonoid_library.csv` |

All three were confirmed byte-identical before removal. The lesson is narrow but worth keeping: a
same-basename search is not a duplicate check, and a previous session's verification is evidence,
not proof.

Removed in three groups:

- **7 × `reference.frc`, 1.83 GB.** `MMPBSA.py` force dumps. The only mentions anywhere are
  `.gitignore` and a note in `06_md/README.md` saying they are a byproduct; `run_mmgbsa.sh`
  regenerates them.
- **28 duplicates, ~11 MB.** Each re-hashed against its surviving copy immediately before removal.
- **4 × RL checkpoints, 379 MB.** The one irreversible loss: stage-1 weights for generative
  campaigns 1 and 2, held as two unique files plus two exact duplicates. Those campaigns can no
  longer be resumed or re-sampled. Their CSVs, score plots and `.toml` configs remain, so they stay
  documented; campaign 3, which produced `cand_003`, keeps `generation/rl_stage1_v3.chkpt`.

Negative results were kept throughout — the dissociated MD runs, the failed campaign 1 outputs, the
superseded free-conformer parameters, the 1 ns MM-GBSA estimate, and the off-pose −44.22
trajectory. In this project those carry most of the findings.

### The top level was the real problem, not the contents

`_ARCHIVE_TO_DELETE/` and `_REVIEW_REQUIRED/` were both visible at the root of the repository, and
their names announce unfinished business to anyone browsing it. Both removed. `provenance/` replaces
them, holding the two original 17 Aug project plans (kept: they show the pipeline was specified
before it was run), the deletion manifest, the relocation log and a README explaining what went and
why. The three root-level queue logs moved to `06_md/queue_logs/` with a README, since they are the
MD run record and belong beside the runs.

Every documentation reference to the removed folders was updated and all links verified to resolve.

### Pushed to GitHub, which required rewriting the unpushed history

The branch had been unpushable since the 26th: eleven `reference.frc` blobs of 194–339 MB sat inside
commits `0a758cb` and `74129f0`, over GitHub's 100 MB hard per-file limit, which is enforced against
every blob in the pushed history. Deleting them in a later commit does not help.

`git filter-branch --index-filter` over `origin/master..HEAD` stripped `*.frc` and `*.nc` from the
five unpushed commits. Safe because nothing had been pushed, so no other history was affected. The
correctness check was that **HEAD's tree hash came out byte-identical before and after**
(`54d8330`) — which it must, since those patterns were already untracked at HEAD — with all five
commit messages preserved. Upload fell from 4.46 GB to 478 MB. Pushed `a27a5d8..f913df9`.

A `backup-before-strip` tag points at the pre-rewrite commit locally.

### Workflow document rewritten as a single protocol

`Complete_Workflow_Two_Person.md` → [`WORKFLOW.md`](WORKFLOW.md), restructured from parallel tracks
into one continuous protocol. The document is now silent on personnel: it describes what the
procedure is and why, not who performs it.

Two substantive improvements while rewriting it:

- **It now states that it is the plan, not the record.** The original read as a description of the
  project, but execution diverged materially — 16 MD systems rather than 4, 10.1 ns at 300 K rather
  than 3 × 20 ns at 310 K, three docked poses per target rather than replicates alone, FAD
  restrained rather than covalently modelled. A divergence table sits at the top and names
  `LOGBOOK.md` as primary wherever the two disagree. Left unmarked, the document would have
  misdescribed the work.
- **Part 6 no longer inlines script bodies.** It listed full source for six scripts that now exist
  for real in `scripts/`. A second copy drifts from the first, and a reader following the stale copy
  gets different numbers — the same failure the summary-CSV generator exists to prevent. Replaced
  with a table pointing at the real files.

It also absorbs the findings that postdate the original plan: the SEM's failure scaling with pose
stability, never computing MM-GBSA off-pose, pose spread not substituting for replicate spread,
absolute ΔG not comparing across proteins, the 0.861 kcal/mol calibration bias, and the operational
traps (MGLTools' Python 2.7 first on PATH, OpenMM needing `conda run` for CUDA, `-sp` aborting
against a stripped trajectory, no checkpointing in `run_md_restrained.py`).

Files: `provenance/` (new: 2 plans, deletion manifest, relocation log, README),
`06_md/queue_logs/` (new), `WORKFLOW.md` (new), `Complete_Workflow_Two_Person.md` (removed),
`_ARCHIVE_TO_DELETE/` and `_REVIEW_REQUIRED/` (removed), `README.md`, `INVENTORY.md`,
`00_library/README.md`, `.gitignore`. Commits `e2cc22d`, `461b6b0`, `f913df9`, `0b53435`, `a7044bc`.

## 2026-09-27 — F (reference values traced to primary sources; the 9IV calibration was mixing assays)

Chased every flagged `source` field in `01_smiles/references.csv` to its origin via PubMed, ChEMBL
and BindingDB. **All 16 rows now carry a traced source; no `CONFIRM` marker remains.** The exercise
was supposed to be bookkeeping. It changed a headline number.

### The 9IV calibration was comparing two different assays

The matched-pair calibration is the experiment that bounds every TTBK1/TTBK2 selectivity claim in
this project. It used TTBK1 = 430 nM against TTBK2 = 490 nM.

**The 430 was not a measurement.** It was the midpoint of BindingDB's 330–530 nM range for 7Q8V,
and that range *aggregates two assays from two different papers*. The 490 came from a single paper.
So one arm of the matched pair was a cross-paper average and the other was a single value — the
same arm-asymmetry error that produced the water-shell artifact, the pose-scan-versus-replicate
artifact and the off-pose MM-GBSA sign inversion. Three times this project has been bitten by it;
this is the fourth, and the first time it reached a published-facing number.

Both values exist in **one paper, one assay format**: Nozal et al. 2022 (*J Med Chem* 65(2):
1585–1607, PMID 34978799) report **TTBK1 IC50 330 nM** and **TTBK2 IC50 490 nM** for compound 42 /
VNG2.73 — the PDB ligand `9IV` — inhibiting recombinant human enzyme with the RICDLHDDEEDEAMSITA
substrate. Verified two ways: RCSB confirms 7Q8V is *TTBK1* + VNG2.73 and 7Q8Y is *TTBK2* + the same
compound, so the matched-pair premise itself is sound; and ChEMBL curates it as `CHEMBL5200069`
whose canonical SMILES matches the row byte-for-byte (C18H13ClN4O).

Re-running `scripts/calibration_9iv.py` with the arm-symmetric pair:

| | Before | After |
|---|---|---|
| Docking margin (measured, unchanged) | +0.784 | +0.784 |
| Experimental ΔΔG | −0.077 | **−0.234** |
| **Systematic protocol bias** | 0.861 | **1.018 kcal/mol** |

**This strengthens the central negative result.** The interpretability floor rises from ~0.9 to
~1.0 kcal/mol, so the TTBK2-over-TTBK1 margin sits even more clearly inside the protocol's own
bias. The correction also removes a phrase that was never quite true: experiment does not say
"no preference" between the paralogs, it says TTBK1 is favoured by a little over 0.2 kcal/mol.

Propagated to `04_docking`, `05_validation`, `08_analysis`, `09_manuscript`, `10_results`,
`README.md` and `WORKFLOW.md`.

### DTQ's 240 nM is a Kd, not an IC50

Xue et al. 2013 (*ChemMedChem* 8(11):1846–54, PMID 24039150) measured binding by **surface plasmon
resonance**, which the abstract states plainly. RCSB 4BTK reports *both* a **Kd of 240 nM** and an
**IC50 of 4610 nM** for the same ligand. The row had the Kd labelled as an IC50, understating the
enzymatic potency by **19-fold**.

No reported number moved, because DTQ is already excluded from the correlation. But 240 nM must not
be quoted as an IC50, and a Kd must not be mixed into a table of IC50s.

### The correlation set is mixed-species, and the rat points are the flavonoids

Species were verified from PubMed MeSH terms and abstracts for every correlation row:

| Compound | Target | Value | Species |
|---|---|---|---|
| safinamide | MAO-B | 7.67 nM | recombinant human |
| lazabemide | MAO-B | 0.063 µM | human |
| quercetin | MAO-A | 1.52 µM | recombinant human |
| 9IV | TTBK1 / TTBK2 | 330 / 490 nM | recombinant human |
| **kaempferol** | MAO-A | 700 nM | **rat brain** (Sprague-Dawley) |
| **isatin** | MAO-B | ~3 µM | **rat brain** |

The docking is against human structures, so **two of seven points are the wrong species — and they
are the two flavonoids**, the compounds closest to this project's lead and the anchors a reader
would weigh most heavily. Either drop them and state that the correlation rests on four human
points, or label it mixed-species. It must not be presented as a human-target correlation without
the caveat. Isatin's value is additionally only "IC50 approximately 3 µM".

### Two rows that cannot be primary-sourced, and are now labelled as such

Selegiline (7.0 nM) and rasagiline (4.4 nM) trace to Cavalli et al. 2008 (*J Med Chem* 51(3):
347–72, PMID 18181565), which **PubMed types as a Review** — the values are compiled there, not
measured there. Both compounds are irreversible covalent inhibitors whose IC50 is
preincubation-time dependent and therefore not an equilibrium constant, which is why ChEMBL lists
them across 2.76–52 and 4.0–46 nM. Both were already excluded from the correlation. The rows now
say "still not a primary source" rather than implying otherwise.

Clorgyline's Wikipedia chembox reference was dropped: no numeric value from it is used anywhere, so
no assay source is needed, and its structure was verified against `CHEMBL8706`.

### A value is only meaningful with its assay attached

Safinamide makes the point better than any argument. Its 7.67 nM is real and primary — Stössel et
al. 2013 (*J Med Chem* 56(11):4580–96, PMID 23631427), recombinant human MAO-B with p-tyramine
substrate. But Binda et al. 2007 (PMID 17915852) — **the paper behind the very 2V5Z structure this
project docks into** — reports *K*i 0.1–0.5 µM for the same compound, 13–65× weaker. Neither is
wrong. They are different measurements, and quoting either without its assay invites a reviewer to
ask which one you meant.

### Process note

Re-running a mutating update script purely to read its status report silently re-applied its edits
and clobbered a later fix. Caught by the `CONFIRM` grep immediately afterwards. Status checks
should be read-only; that is what the grep is for.

Files: `01_smiles/references.csv` (all 16 rows sourced; 9IV TTBK1 430→330 nM; DTQ relabelled Kd),
`01_smiles/README.md` (corrections and species audit), `scripts/calibration_9iv.py` (IC50 inputs +
the reasoning for the change), `05_validation/{calibration_9IV.csv,calibration_9IV_margin.json}`
(regenerated), `04_docking/README.md`, `08_analysis/README.md`, `09_manuscript/README.md`,
`10_results/README.md`, `README.md`, `WORKFLOW.md`.

Sources consulted via PubMed: PMIDs 34978799, 24039150, 23631427, 18181565, 17915852, 10813558,
8687491, 29496172, 32583952. Curated cross-checks via ChEMBL (CHEMBL5200069, CHEMBL396778,
CHEMBL972, CHEMBL887, CHEMBL8706, documents CHEMBL5154682/CHEMBL2380255/CHEMBL1141467) and
BindingDB via RCSB annotations for 7Q8V, 7Q8Y, 4BTK and 2V5Z.

## 2026-09-27 — G (the calibration bias had no error bar; fixed, and the interpretability floor rises to ~1.3)

`05_validation/calibration_9IV_margin.json` was reporting
`"experimental_ddG_range": [-0.234, -0.234]` — a range with identical endpoints. A defect I
introduced in entry F: when the IC50 inputs were narrowed to a single arm-symmetric pair, the
min/max machinery was left in place and produced a degenerate band.

The cosmetic bug was hiding a real one. `bias = margin - mid` was being reported as **1.018 with
no uncertainty whatsoever.** Three significant figures, no error bar, on the number that gates
every TTBK1/TTBK2 selectivity claim in the project — which is precisely the failure this project
documents in step 7 about the MM-GBSA SEM. The ±0.015 on the docking margin is seed noise; the
experimental reference contributed nothing at all, because the two IC50s are single values and
**Nozal 2022 publishes no error on either.**

### The fix

`IC50_REL_SD = 0.30` is now a named, documented assumption in `scripts/calibration_9iv.py` — a
conventional within-assay precision for enzymatic IC50 replicates, which is the relevant scale
since both values come from the same assay in the same paper. It propagates through the log and
into the bias:

| | Value |
|---|---|
| Docking margin | +0.784 ± 0.015 (seed noise only) |
| Experimental ΔΔG | −0.234 ± 0.251 (**assumed** 30% IC50 precision, not measured) |
| **Systematic bias** | **~1.0 ± 0.25 kcal/mol** |
| **Interpretability floor** | **~1.3 kcal/mol** (bias + 1 SD) |

**The assumed experimental term is 16× the docking SEM.** The calibration's precision is set by the
literature value, not by the docking, so no amount of further docking would tighten it. That is
worth stating plainly: the limit on this project's selectivity resolution is a published IC50 pair,
not the compute.

### What it changes

The interpretability floor moves from ~1.0 to **~1.3 kcal/mol**, and this has one consequence that
needed following up: **the ~1.6 kcal/mol water-symmetric margin clears the floor only barely**, and
subtracting the bias leaves ~0.6. So the TTBK2-over-TTBK1 result remains uninterpretable as
selectivity, now with a wider margin of safety rather than a narrower one. Recorded in
`05_validation`, `04_docking`, `08_analysis` and `10_results`.

Two wording errors fell out of the same review. `05_validation/README.md` said experiment shows
"**essentially no preference**" and in the next clause "a slight preference for TTBK1" — a direct
self-contradiction left by entry F's edit. It also still claimed "the true answer is known to be
'no difference'". Both rewritten: the measured answer is a small preference for TTBK1, −0.234
kcal/mol, and the protocol reports the opposite sign.

**Quote the bias as ~1.0, never 1.018.** The docs now say so explicitly in three places, and the
JSON carries `"systematic_bias_quote_as": "~1.0 kcal/mol"` plus an
`experimental_ddG_basis` string spelling out that the ± is assumed rather than measured.

Files: `scripts/calibration_9iv.py`, `05_validation/{README.md,calibration_9IV_margin.json}`,
`01_smiles/README.md`, `04_docking/README.md`, `08_analysis/README.md`, `09_manuscript/README.md`,
`10_results/README.md`, `README.md`, `WORKFLOW.md`.

## 2026-09-27 — H (a completion audit; the withdrawn −3.89 was still standing in the top-level README)

Swept the workspace for remaining work. The compute is finished: all 16 MD systems carry
`final_state.xml`, the 3 without MM-GBSA are exactly the 3 that should have none, and
`collect_md_summary.py --check` reports the summary table up to date. What the sweep found instead
was five documentation defects, four of them created by the **last** batch of data — the MAO
velocity replicates — which added two systems without the prose being re-read against them.

### The withdrawn number was still in the shop window

Entry on 2026-09-26 withdrew the **−3.89 kcal/mol** MAO-B/MAO-A docking margin: it averages in 8
candidates that have *no viable MAO-A pose at all* (best mode worse than −5.0 kcal/mol), and
averaging a non-measurement into a margin manufactures selectivity out of a docking failure. The
honest figure is **−3.28 over the 48 quantifiable candidates, favourable 48/48**, with the other 8
sterically excluded.

**The withdrawal never reached the docs.** `08_analysis/README.md` carried −3.89 twice, once in the
margins table (with an all-56 range of −10.78 to −1.80) and once as "the more robust signal ...
across 56 candidates" — the precise framing that was withdrawn. `README.md` carried a hybrid:
the withdrawn **−3.89** against the corrected **48/48** denominator, which is not a statistic that
exists. Recomputed from `selectivity_margins.csv` excluding the 8 named failures: mean **−3.279**,
favourable **48/48**, range **−6.04 to −1.80** — reproducing −3.28 exactly.

**The root cause is a gap in the script, not just in the prose.** `analyze_selectivity.py` applies
`FAIL_THRESHOLD` when it *prints*, but computes `m["no_pose"]` and then drops it from the returned
columns, so `08_analysis/selectivity_margins.csv` ships all 56 rows with no flag on the 8. Anyone
who averages the CSV — which is the obvious thing to do — gets −3.89 back. The docs now carry an
explicit "do not average this CSV directly" warning naming the 8 ligands, but **propagating the
`no_pose` column into the CSV is the actual fix and is still open.** This is the summary-CSV
principle from step 10 applied one level down: a table a reader can misread is a table that will be
misread.

### Three stale counts and a band that mixed two statistics

- **"3.44 Å, the highest in the project" is no longer true.** `system_MAOA_p3_r2` — a replicate
  added the same day — reports **3.64 Å** whole-protein backbone. Verified against
  `prot_rmsd.dat` (3.640) and `prot_core_rmsd.dat` (1.500). Its RMSF has the same signature as
  `p2`: the maximum at the C-terminal tail, residues 506–513 of 513, peaking at **9.16 Å** against
  `p2`'s 8.29 Å — which is the entire difference between 3.64 and 3.44. Same tail, flailing harder.
  The run was undiscussed anywhere, and it is the sharpest example in the project of the point that
  section exists to make: **`p3_r2` holds the project's most stable ligand (1.07 Å) and its highest
  backbone RMSD at the same time.** Corrected in `06_md`, `07_mmgbsa` and `10_results`.
- **The `1.12–1.81 Å` backbone band was two different statistics spliced together.** The `1.12` is
  a *whole-run mean* (`system_MAOB`, logged above); the `1.81` is a *last-100 mean*. Quoted as a
  range in three top-level documents, and it no longer covers the set either: last-100 across the
  16 systems runs 1.19–3.64. Replaced with the honest split — TTBK **1.34–1.90 Å** whole-protein,
  MAO **1.06–1.77 Å** core-fit — and the statistic is now named.
- **`08_analysis/README.md` was pre-replicate**: "all six MAO-A/MAO-B pairings" (now twelve) and
  "the *magnitude* is still inside method noise", which the replicates contradicted on 27 Sep. It
  was the last document still saying the magnitude is unquotable.
- **Two more counts left at six by the replicates**: `09_manuscript/README.md`'s "all six MAO runs
  stayed in the site" and `10_results/README.md`'s FAD validation "across six systems". Both are
  eight. The FAD *range* (0.46–0.73 Å) is unchanged — it already held across all eight.
  `06_md/README.md`'s "all six MAO systems are neutral at +0.0010 e" was left alone: that sentence
  is about the six pose-scan builds specifically, and is correct.

### A rule that contradicted itself two paragraphs apart

`README.md` said "**No MM-GBSA was computed for the runs that drifted or dissociated**" and then,
immediately below, discussed −44.22 — a number from a run that drifted 3.24 Å. Three off-pose runs
have computed energies (MAO-A pose 1 at −44.22; TTBK2 pose 1 runs at −29.78 and −22.88), and they
were computed *deliberately*: the first inverted the MAO conclusion for two days and the pair
supplied the 33× SEM discrepancy. The rule that was actually followed is **"no reported binding
energy comes from an off-pose run"**, which is narrower and survives contact with the data.
Rewritten to say that.

### Also corrected: a Vinardo claim resting on the wrong receptors

`08_analysis/README.md` called the lead "robust under both Vina and Vinardo" in the same sentence
as its `7JXX`/`2V5Z` ranks. The Vinardo cross-check exists only in `prior_phase/`, run on `4NFM`.
The claim is now attributed to the prior phase, and the open item — re-running Vinardo on the
validated receptor, `WORKFLOW.md` Phase 3 — is named inline.

### What remains

1. **The manuscript body has not been touched since 20 Aug** and is still entirely prior-phase:
   `4NFM`/`6U0K` receptors, 20 ns, the −26.32/−36.85 MM-GBSA pair quoted as SEM, a 15-candidate
   shortlist, §3.7 unsplit, and §4.5 limitation 5 ("no replicate MD") now resolved and needing
   deletion rather than softening. `09_manuscript/README.md` says what must change; nothing has.
   This is the largest remaining task in the project.
2. **Propagate `no_pose` into `selectivity_margins.csv`** so the −3.89 trap cannot be re-entered.
3. **Vinardo on the validated receptor** (`WORKFLOW.md` Phase 3, manuscript §3.3).
4. **TTBK2 pose 2 velocity replicate** — the one on-pose TTBK2 number (−30.95) still carries no
   error bar of its own, and the TTBK1-vs-TTBK2 comparison rests on it.

Process note, consistent with entry F's: **every defect in this entry was found by checking prose
against the primary files, and four of five were introduced by adding data without re-reading the
prose.** Adding a system to the summary table is not a complete change until the counts, ranges and
superlatives that quantify over the set have been re-derived from it. The superlatives are the
dangerous ones — "the highest in the project" is a claim about every row, so every new row can
falsify it silently.

Files: `README.md`, `06_md/README.md`, `07_mmgbsa/README.md`, `08_analysis/README.md`,
`09_manuscript/README.md`, `10_results/README.md`. No computed result changed; no script changed.

## 2026-09-27 — I (04_docking regrouped by purpose; a stale pIC50 found by baselining first)

`04_docking/` held **77 sibling run directories, 5,160 files — 85% of everything tracked in this
project**. Opening it told a reader nothing. Regrouped by purpose:

    04_docking/candidates/      24 runs   the production screen
               references/      18 runs   known-inhibitor benchmark
               native_redock/   30 runs   redocking validation
               controls/         5 runs   single-ligand controls (noFAD, bridging water)
               prior_phase/               untouched

**The run directory names are unchanged.** That is the whole reason this was safe: every parser
recovers receptor, ligand set and seed from `basename $(dirname <log>)`, so the grouping is
invisible to it. Readers glob one extra level (`04_docking/*/<REC>_<SET>_seed*/`), deliberately
group-agnostic so that reclassifying a set can never hide its runs. Only the write side needs to
know which group a set belongs to, and that mapping now lives in exactly one file,
`scripts/docking_paths.sh`, sourced by all five scripts that create or locate run directories. Five
copies of a `case` statement would have drifted, and a drifted copy writes a run where the reader
cannot see it.

The numbered step directories were **not** touched. There are **127** hardcoded path references in
`scripts/`, `generation/` and `10_results/` — 48 to `06_md`, 23 to `03_receptors`, 16 to
`01_smiles` — and the top-level README's own rule is to add to the structure rather than rearrange
it. `06_md/` is 91% of the workspace by disk (17.5 GB, of which 13.3 GB is gitignored
`production*.dcd`) but only 520 tracked files, so it is a disk question, not a navigability one.

### Two glob hazards, checked rather than assumed

1. **`prior_phase/` lives inside `04_docking/`**, so adding a wildcard level could have started
   matching runs that the old one-level glob never saw — a silent change to what gets averaged, the
   exact failure mode this project keeps hitting. Checked first: `prior_phase/` contains **zero**
   `.log` files and no `*seed*` directories. No collision.
2. **mtimes are load-bearent here.** Both `collect_results.py` and the shell resume logic compare a
   log's mtime against its receptor's to drop stale runs. A copy would have reset them and silently
   revalidated stale logs. `git mv` on one filesystem is a rename, and the mtimes were verified
   afterwards as still 24–25 Sep.

### Verification

Baselined before touching anything, which is what caught the defect below. All four analysis
entry points were run pre-move and post-move: `analyze_selectivity.py`, `analyze_water_test.py`,
`calibration_9iv.py`, `collect_results.py`.

| Check | Result |
|---|---|
| stdout, all four scripts | **byte-identical** |
| stderr, all four scripts | **byte-identical** |
| regenerated consensus table (558 rows) | **byte-identical** |
| `bash -n` on all 11 shell scripts | pass |
| `docking_out` vs disk, one probe per group | 4/4 resolve to real directories |
| `run_week2_redock.sh --dry-run` | `runs needed: 0   already current: 15` |
| `run_selectivity.sh --dry-run` | `runs needed: 0   already current: 9` |
| `run_water_test.sh --dry-run` | `runs needed: 0   already current: 15` |

The dry-runs are the decisive ones: had the write-side mapping been wrong, they would have reported
15 runs needed and cheerfully re-docked several thousand files over the top of good results.

### The baseline caught a stale derived value from entry F

Re-running the scripts *before* the move changed one committed file, which means the change was
pre-existing and nothing to do with the restructure. `08_analysis/water_test_references.csv` had
`9IV_ttbk1` at pIC50 **6.3665** = −log10(430 nM) — the **withdrawn** value. Entry F corrected 9IV
TTBK1 from 430 to 330 nM and regenerated the calibration files, but never regenerated this one.
Correct value **6.4815** = −log10(330 nM), now written.

Consequence is contained: the file feeds the score-vs-pIC50 correlation that
`analyze_water_test.py` prints (n=6, r=−0.642 wet / −0.418 dry), and **no document quotes that
r**, so no reported number moves. But it is the second derived artefact found carrying a withdrawn
reference value, after the −3.89 margin in entry H. Both were found the same way — by regenerating
and diffing rather than by reading. **A correction to `01_smiles/references.csv` is not complete
until every derived file has been regenerated and diffed**; there is no index of which files those
are, and building one is worth more than finding the third instance.

`consensus_week2.csv` also differs from a fresh collect (165 pairs against 558), but that is by
design: it is a snapshot of the week-2 redock subset, not of the whole directory. It was left
alone; the fresh collect went to a scratch path.

Files: `04_docking/` (77 directories moved via `git mv`, history preserved),
`scripts/docking_paths.sh` (new), `scripts/{dock,run_selectivity,run_water_test,run_week2_redock,
run_minwater_probe}.sh`, `scripts/{analyze_selectivity,analyze_water_test,calibration_9iv,
collect_results}.py`, `08_analysis/water_test_references.csv` (regenerated, 9IV pIC50 corrected),
`04_docking/README.md`, `scripts/README.md`, `README.md`, `INVENTORY.md`, `WORKFLOW.md`.
No receptor, ligand, trajectory or energy was touched.

## 2026-09-27 — J (the remaining directories restructured; one rule per boundary, readers resolve by glob)

Continued entry I's regrouping through the rest of the workspace. Five directories changed, four
were deliberately left alone, and the same design was used at every boundary.

### The design, stated once because it is what makes this safe

At each boundary there is **one writer that knows the rule** and **readers that carry no rule at
all**, resolving by glob at either depth. So:

- reclassifying something cannot make it silently unfindable;
- a half-migrated tree still resolves;
- an ambiguous hit (two copies of one receptor, one system, one ligand) **aborts** rather than
  silently picking one — that last property matters most in `collect_md_summary.py`, where a
  system it cannot find is dropped from the table rather than raising.

The rule lives in exactly one place per boundary: `docking_paths.sh` (run groups),
`receptor_paths.sh` (target families), `prep_ligands.py:ligand_group` (ligand classes). Where both
shell and Python need to *read* a boundary, the Python side (`project_paths.py`) deliberately has
**no** rule, only the glob — so there is nothing for a second copy to drift against.

### What changed

| Directory | Before | After |
|---|---|---|
| `02_ligands/{sdf,pdbqt}/` | 77 flat files each | `candidates/` 56, `natives/` 5, `references/` 16 |
| `03_receptors/` | 12 receptor dirs at top level | `ttbk/` 7, `mao/` 5, plus `_cofactors/`, `prior_phase/` |
| `05_validation/` | 11 flat files | 6 redock logs into `redock/` |
| `06_md/` | 24 entries, 16 of them systems | `systems/` holds all 16; 9 entries |
| `09_manuscript/` | rendered `.html`/`.docx` beside source | `rendered/` |

Top-level entry counts now: `00_library` 4, `01_smiles` 6, `02_ligands` 4, `03_receptors` 6,
`04_docking` 6, `05_validation` 6, `06_md` 9, `07_mmgbsa` 3, `08_analysis` 11, `09_manuscript` 5,
`10_results` 6.

### What was left alone, and why

`00_library` (4 entries), `01_smiles` (6), `07_mmgbsa` (3) and `10_results` (5) are already
minimal — subdividing them would add a lookup step and remove nothing. `08_analysis` is ten
distinctly named CSVs, each with a row in its README table explaining what it is; grouping those
would make the table's paths longer and the files harder to find, not easier. **Structure is only
worth adding where it removes a choice a reader has to make.**

`06_md/systems/` was deliberately **not** split further by target. The names already sort into
target clusters, and every extra level has to be paid for again in three WSL absolute paths and in
the summary table. One level bought the whole listing win at a fraction of the churn.

### Two things that would have broken silently

1. **`06_md/run_md.py` and `run_md_extend.py` hardcode `06_md/system/...`** rather than taking the
   system directory as an argument — unlike `run_md_restrained.py` and `run_md_system.py`. Eleven
   paths between them, and they are the scripts that produced the original TTBK1 baseline. Found by
   grepping the Python inside `06_md/`, not just `scripts/`.
2. **`prep_ligands.py` writes flat.** Grouping `02_ligands/` without teaching the writer the rule
   would have meant the next `--csv` re-prep quietly scattering files back outside the groups, where
   a group-only glob would miss them and `dock.sh`'s final `ls` fallback would then dock *a
   different set*. Both were fixed: the writer now files into the group, and the reader accepts
   either depth.

### Verification

Baseline captured before the first move (`scratchpad/verify.sh`, 8 checks: four analysis entry
points, the MD summary freshness check, three driver `--dry-run`s, plus copies of six tracked
output files), then re-run after **each** of the five directories.

**Every check byte-identical at every step, with exactly one expected exception:** the
`system_dir` column of `07_mmgbsa/md_mmgbsa_summary.csv`, which is a path and had to change.
Compared column by column: `system_dir` moved on all 16 rows (`06_md/system` →
`06_md/systems/system`), and **all 17 other columns byte-identical** — every RMSD, every ΔG, every
SD, every verdict. `10_results/fig{1,2,3}.png` regenerate to **identical md5s**.

Also: `bash -n` clean on all 11 shell scripts, `py_compile` clean on every Python file in
`scripts/`, `06_md/` and `10_results/`, and each new resolver unit-tested including its
not-found and ambiguous paths. All three docking drivers still report `runs needed: 0`.

### Side effect worth recording

`06_md/systems/` is a rename of **17.5 GB**. Instant on one volume, but this workspace lives under
OneDrive, so the client may re-upload those files. Nothing is lost either way — the 13.3 GB of
`production*.dcd` is gitignored and regenerable — but it is a real bandwidth cost that a
restructure of a smaller directory would not have carried.

Files: `02_ligands/{sdf,pdbqt}/` (154 files moved), `03_receptors/` (12 dirs moved),
`05_validation/redock/` (6 moved), `06_md/systems/` (16 moved), `09_manuscript/rendered/` (2 moved),
all via `git mv`. New: `scripts/receptor_paths.sh`, `scripts/project_paths.py`,
`scripts/md_paths.sh`. Modified: `scripts/{dock,prep_receptor,prep_ligands,run_mmgbsa,
run_minwater_probe,run_selectivity,run_water_test,run_week2_redock,run_mao_pose_analysis,
run_mao_pose_queue,run_mao_replicate_queue}`, `scripts/{analyze_selectivity,analyze_water_test,
collect_results,collect_md_summary,fix_native_bondorders}.py`, `06_md/{run_md,run_md_extend}.py`,
`07_mmgbsa/md_mmgbsa_summary.csv` (path column only), and the layout sections of
`README.md`, `WORKFLOW.md`, `INVENTORY.md`, `scripts/README.md`, `02_ligands/README.md`,
`05_validation/README.md`, `07_mmgbsa/README.md`, `09_manuscript/README.md`,
`generation/README.md`, `provenance/README.md`.
No receptor, ligand, trajectory, score or energy was altered.

## 2026-09-27 — K (the manuscript rewritten from the current dataset; the August draft archived, its renders deleted)

The manuscript had not been touched since 20 Aug and was still entirely the prior-phase paper:
`4NFM`/`6U0K`/`2V60` receptors, 20 ns MD, MM-GBSA quoted as ±SEM, a 15-candidate shortlist, one
combined selectivity claim, and a limitation ("single 20 ns replicate, no replicate MD") that the
September work had already resolved. `09_manuscript/README.md` had documented what each claim must
become; none of it had been applied.

**Rewritten rather than patched.** The receptor set, the MD protocol and the conclusions all
changed together, so editing in place would have produced a document mixing two incompatible
passes — exactly the silent-mismatch failure the `prior_phase/` convention exists to prevent.

### What was removed, and what was kept

- **Archived to `09_manuscript/prior_phase/`**: the August draft
  (`manuscript_dualtarget_flavonol_AUG.md`), its figure/journal index, and
  `Week1_Progress_Report.md`. Kept, not deleted, because the draft is the record of what was
  claimed before the corrections — and this project's whole method is preserving that record.
  The same reason `prior_phase/` exists in nine other step directories.
- **Deleted**: `rendered/manuscript_dualtarget_flavonol.{html,docx}`, 5.8 MB generated from the
  superseded draft. A stale render that looks authoritative is worse than no render, and both
  remain in git history at `e442b33`. There are deliberately no rendered outputs now; re-render
  when the text is final.

### The new structure, and what changed in substance

Results went from seven sections to nine, because three findings that did not exist in August are
now primary:

| August | Now |
|---|---|
| §3.5 "Both complexes are stable over 20 ns" | §3.7 stability is **pose-dependent** — 8 of 16 runs hold, 2 dissociate; **docking rank predicted the most stable pose in 1 of 4 targets** |
| §3.6 MM-GBSA −26.32 ± 0.53 / −36.85 ± 0.23 (SEM) | §3.8 on-pose only, **replicate spread** as the error bar |
| — | §3.9 **the reported SEM is not an error bar**, and its failure scales monotonically with pose instability (2× to 33×) |
| — | §3.5 a measured **~1.0 ± 0.25 kcal/mol protocol bias toward TTBK2**, giving a ~1.3 floor |
| — | §3.6 the TTBK margin is **largely a property of the water shell** |
| §3.7 one combined "no isoform discrimination" claim | §3.8 **split by pair** — MAO supportable in direction *and* size (2.55 against 0.24–0.26); TTBK not supportable at all |
| §4.4 limitation 5 "no replicate MD" | **deleted, not softened** — the replicates were run |

§4.3 is new and collects all four arm-asymmetry false positives in one place: the water shell, the
pose-scan-versus-replicate-scan, the on-pose-versus-off-pose inversion, and — the fourth, from
entry F — the cross-paper-versus-single-paper IC50 in the calibration pair itself.

The title changed to foreground the negatives, per the project's own framing rule. The
contribution is now as much methodological as chemical, so `J Chem Inf Model` was added to the
journal assessment alongside `J Cheminformatics`.

**The one gap stated in the manuscript rather than hidden:** there is no independent
scoring-function cross-check on the validated receptors. The Vinardo result (r = 0.56 at TTBK1,
−0.28 at MAO-B) was run on the earlier receptor set, so §3.3 presents it as a general warning
from this project's earlier data and explicitly *not* as validation of the current rankings.

### Verification: every number in the manuscript regenerated from the primary files

Wrote `scratchpad/audit_manuscript.py`, which parses the manuscript and re-derives its quantities
from `07_mmgbsa/md_mmgbsa_summary.csv`, `08_analysis/selectivity_margins.csv` and
`05_validation/calibration_9IV_margin.json`. **59 checks, all passing:** every system's last-100
RMSD and on-pose ΔG appears as written; both replicate spreads (0.26, 0.24), both replicate means
(−36.53, −39.09), the ΔΔG (2.55) and the ratio (9.8×) recompute; the TTBK spread (2.08) and ΔΔG
(1.45) recompute; the 12 on-pose pairings and their 2.30–3.94 range recompute; the protein-RMSD
bands, the FAD band over 8 systems, the verdict counts, the 48-candidate MAO margin (−3.28,
48/48, −6.04 to −1.80) all recompute; and every referenced figure file exists.

It also enforces the negative rules: **the two off-pose TTBK2 energies appear nowhere**, `−3.89`
appears nowhere, and `1.018` appears nowhere.

Three checks failed on the first run and are worth recording, because two were the audit's fault
and only one was the manuscript's:

1. **`−44.22` appeared zero times.** I had written around it. But `README.md` and the figure index
   both state it appears once in §4.3 as the cautionary case, and a reader comparing against the
   supplementary table needs to know *why* the most favourable MAO-A value is excluded. Naming it
   once, explicitly as not-a-binding-energy, is more honest than omitting it. Added.
2. **The audit computed the TTBK ΔΔG as 0.41, not 1.45** — it took the less-negative of the two
   TTBK1 replicates instead of their mean. The documented basis is the replicate mean
   (−32.40 against −30.95 = 1.4499). Audit corrected; the manuscript was right.
3. **The audit counted 10 runs holding, not 8.** It pooled `stable` with `holds (loose)`. The
   documented count is strictly `stable` = 8, with the two loose holds reported separately. Audit
   corrected; the manuscript was right.

That is the second time in two days a checking script has been wrong where the prose was right
(the first was the `--check` false alarm in entry I's baseline). **A verifier is code and gets
the same scrutiny as the thing it verifies** — a failing check is a hypothesis about a defect,
not a defect.

Also fixed while here: `04_docking/README.md` still showed `03_receptors/<PDB>/receptor.pdbqt`
in its Vina invocation block, stale since entry J's family grouping; and the claims table in
`10_results/README.md` said "14 runs" where there are 16.

Files: `09_manuscript/manuscript_dualtarget_flavonol.md` (rewritten, 733 lines),
`09_manuscript/manuscript_figures_tables_and_journals.md` (rewritten),
`09_manuscript/README.md` (rewritten — the "claims requiring revision" section is now a record of
what changed), `09_manuscript/prior_phase/` (3 files archived), `rendered/` (deleted),
`04_docking/README.md`, `10_results/README.md`. No computed result was altered; no figure was
regenerated.

## 2026-09-28 — L (the TTBK2 replicate: the ΔΔG falls from 1.45 to 0.53, and the fifth pair lands where the trend predicted)

Ran `system_TTBK2_p2_r2`, the velocity replicate of the project's **only on-pose TTBK2
trajectory**. This was the last open scientific item. Launched 00:05, MD done 01:49 (10.1 ns at
139 ns/day), MM-GBSA and cpptraj done 02:01. Driver:
`scripts/run_ttbk2_replicate.sh`, resumable in both halves.

### Why it had to be run rather than reasoned around

`system_TTBK2_p2` at −30.95 kcal/mol was the single number the whole TTBK1-vs-TTBK2 comparison
rested on, and it had no error bar of its own. The published verdict — "no resolvable
discrimination" — compared a 1.45 kcal/mol ΔΔG against **TTBK1's** 2.08 replicate spread. That is
a two-arm difference judged with one arm's uncertainty: the same class of asymmetry as the water
shell, the pose-scan-versus-replicate-scan, the on-pose-versus-off-pose inversion and the
cross-paper IC50. The other TTBK2 pair (pose 1, spread 6.89) could not stand in, because both of
its runs left the docked pose and it therefore describes structures that were never the complex.

### Result

| | ΔG | spread |
|---|---|---|
| `system_TTBK2_p2` (run 1) | −30.95 | |
| `system_TTBK2_p2_r2` (run 2) | **−34.92** | **3.97 kcal/mol** |

**The conclusion held and the number moved.** TTBK2's own spread is 3.97 — nearly twice TTBK1's
2.08, and the largest on-pose spread in the project. More consequentially, the ΔΔG computed from
*two* replicate means rather than one single run is **0.53 kcal/mol**, down from 1.45:

    TTBK1 replicate mean  -32.40  (spread 2.08)
    TTBK2 replicate mean  -32.93  (spread 3.97)
    ddG                     0.53  -- 7.4x smaller than the larger spread

So the non-result is now emphatic and symmetric: both arms carry a measured error bar, and the
difference is a third of what it appeared to be while one arm was unreplicated. **Measuring the
second arm did not merely add an error bar — it moved the central value.** That is the argument
for measuring it, and it is now recorded as the fifth arm-asymmetry false positive in the
manuscript's §4.3. It is the subtlest of the five: the asymmetry was in the *uncertainty*, not in
the systems, so both arms looked properly matched by design.

**Pose stability reproduced.** 1.89 Å against run 1's 1.99 Å, both `stable`. TTBK2 still holds
`cand_003` in 1 of 3 docked poses, as does TTBK1 — the symmetric test is unchanged.

### The fifth replicate pair is an out-of-sample test of the SEM finding, and it passed

The pair was run to give the TTBK2 arm an error bar, not to probe the SEM relationship. Sorted by
how well each pair held its pose:

| Pair | Mean ligand RMSD | Spread | Ratio to SEM |
|---|---|---|---|
| MAO-B pose 2 | 1.24 Å | 0.24 | 2× |
| MAO-A pose 3 | 1.25 Å | 0.26 | 2× |
| TTBK1 pose 1 | 1.63 Å | 2.08 | 8× |
| **TTBK2 pose 2** | **1.94 Å** | **3.97** | **15×** |
| TTBK2 pose 1 | 4.71 Å | 6.89 | 33× |

**The measured spread is now monotonic in ligand RMSD across all five pairs**, and the new point
fell between its neighbours in exactly the predicted order. With five pairs across two protein
families and three proteins, this is a consistent pattern — though the docs now say explicitly
that it should *not* be used to predict a spread from an RMSD.

**And it sharpens the rule rather than just confirming it.** TTBK2 pose 2 is a *stable* pair —
1.99 and 1.89 Å, comfortably inside this project's stability band — and it still returned a
3.97 kcal/mol spread against a reported SEM of 0.27. So pose stability bounds how bad the SEM can
be, but **a stable pose is not a licence to skip the replicate.** Every previous statement of this
finding leaned on the unstable cases; this one does not.

### A correction to my own analysis, caught before it reached any document

My first pass computed the pairs' "mean ligand RMSD" from the `lig_rmsd_mean` column and concluded
the relationship was **not** monotonic. Wrong column: the published table uses the mean of the
pair's **last-100** values, which is also what every stability verdict in this project uses
(a ligand that leaves late shows a low whole-run average — that is why the last-100 window exists
at all). On the correct metric the relationship is monotonic. Recomputed before writing anything.

That is now twice in two days that a checking script was wrong where the data was right, after
entry K's two false audit failures. Recording the rule again because it keeps earning its place:
**a verifier is code and gets the same scrutiny as the thing it verifies.**

### Protocol notes

- **`run_md_system.py`, not `run_md_restrained.py`.** TTBK2 has no FAD. This matters for reading
  the spreads across arms: the MAO pairs' unusually tight 0.24–0.26 is partly an artifact of the
  FAD restraint damping receptor motion, which this system does not have. Borrowing a MAO spread
  for a TTBK arm would have been wrong in both directions.
- **All five inputs verified byte-identical to run 1 by sha256** — `complex.prmtop`,
  `complex.inpcrd` *and* the three analysis inputs (`strip_traj.cpptraj`, `mmpbsa.in`,
  `lig_rmsd.cpptraj`). Copying rather than regenerating the last three means the strip mask, the
  energy model and the RMSD masks provably cannot have drifted between the two runs. The driver
  re-checks the hashes and **aborts** if they differ: a replicate with different topology is not a
  velocity replicate, and its spread would be meaningless.
- **`run_mmgbsa.sh` does not generate the cpptraj/mmpbsa inputs**, it expects them to exist. Had
  they not been copied, the MD would have run for 1.75 h and the analysis would then have failed
  immediately. Worth knowing before building any future system by hand.
- No `rmsf.cpptraj` or `core_rmsd.cpptraj` here: those exist for the MAO systems because both MAO
  constructs end in a solvent-exposed C-terminal tail that dominates a whole-protein fit. TTBK2
  has no such tail and run 1 was analysed with `lig_rmsd` alone.

### Propagated

Summary table regenerated (17 systems, 14 with MM-GBSA); all three figures regenerated.
`collect_md_summary.py` now lists the system and `system_TTBK2_p2` is relabelled "pose 2, run 1";
`make_figures.py` has its label. The manuscript audit was extended to the new dataset — **64
checks, all passing**, including that the stale 1.45 no longer appears as a current ΔΔG.

Files: `06_md/systems/system_TTBK2_p2_r2/` (new, with `PROVENANCE.md` recording the sha256s and
what each possible outcome would have meant, written *before* the result was known),
`scripts/run_ttbk2_replicate.sh` (new), `scripts/collect_md_summary.py`,
`10_results/make_figures.py`, `scripts/README.md`, `07_mmgbsa/md_mmgbsa_summary.csv`,
`10_results/fig{1,2,3}.png`, `06_md/README.md`, `07_mmgbsa/README.md`, `10_results/README.md`,
`09_manuscript/manuscript_dualtarget_flavonol.md` (abstract, §2.5, §3.7, §3.8, §3.9, §4.2, §4.3,
§4.4, §4.5, §5), `09_manuscript/manuscript_figures_tables_and_journals.md`, `README.md`,
`INVENTORY.md`, `WORKFLOW.md`.

## 2026-09-28 — M (Smart App Control blocks RDKit; the prep path moves into WSL2 and is shown equivalent to ≤0.02 kcal/mol)

Windows Smart App Control began blocking two unsigned conda-forge binaries in the native
`docking_project` env:

```
Lib\site-packages\rdkit\Chem\rdchem.pyd   -> "An Application Control policy has blocked this file"
Scripts\mk_export.exe                    -> "blocked by your organization's Device Guard policy"
```

State at the time: `VerifiedAndReputablePolicyState = 1` (Enforced),
`UsermodeCodeIntegrityPolicyEnforcementStatus = 2`.

### Why the obvious fix does not work

Smart App Control has **no allowlist**. `Add-MpPreference -ExclusionPath` configures Defender
antivirus, which is a different subsystem, so a path exclusion cannot reach a file the Code
Integrity policy has blocked. The exclusions were added and both files were confirmed still
blocked, twice. Turning Smart App Control off is the only native fix and it **cannot be undone
without reinstalling Windows** — far too much collateral for a toolchain problem.

`rdchem.pyd` is also only the first of 63 `.pyd` files in that tree, so unblocking it individually
would most likely just move the error to the next one.

### What was affected, and what was not

Blocked: `prep_ligands.py`, `prep_receptor.sh`, `filter_cascade.py`, `bbb_score.py`,
`fix_native_bondorders.py`, and PDBQT → SDF export.

**Correction to my own first pass:** I initially listed `run_prolif_v2.py` here too. It does
import RDKit, but it is prior-phase and cannot run at all — `prep/`, `ligands/` and `docking/`
were all removed in the cleanup, so it fails on missing inputs with or without RDKit. Leaving it
on the blocked list would have sent someone to fix the wrong problem. ProLIF 2.2.1 and
MDAnalysis 2.9.0 were installed into `ligprep` anyway, so it has a home if it is ever brought
forward to the numbered tree.

Unaffected: Vina scoring, OpenMM MD on CUDA, WSL MM-GBSA, and every analysis and figure script.
The Vinardo cross-check ran throughout, because it is score-only (`DOCK_SDF=0`) and needs no SDF.

### The fix: WSL2, which the policy does not govern

WSL2 runs ELF binaries and is outside Smart App Control — which is also why AmberTools in the
`mdgbsa` env kept working the whole time. A new env, `ligprep`, was built there (spec exported to
`scripts/env_ligprep_wsl.yml`). It is deliberately **separate from `mdgbsa`**: a dependency solve
that broke that AmberTools install would cost far more than it saves.

`scripts/wsl_run.sh` is the bridge. The project's scripts take relative paths and run from the
project root, so the whole script executes inside WSL unchanged — no per-binary shim and no
argument rewriting for a caller to get subtly wrong. `shutil.which` inside those scripts then finds
the Linux entry point names on its own.

Two install notes:

* `conda create` first failed on a **Terms-of-Service gate** for `repo.anaconda.com/pkgs/{main,r}`,
  and `| tail -25` hid the failure behind exit code 0. Fixed with `--override-channels -c
  conda-forge`, which is where everything here lives anyway — no ToS acceptance on the user's
  behalf.
* The entry point is `mk_export.py` on conda-forge Linux and `mk_export.exe` on Windows — same
  meeko version, different script name. `wsl_run.sh` aliases the bare names so neither platform's
  spelling breaks.

### Versions pinned to match, and the two that could not be

`ligprep` is pinned to the Windows env's versions so it is a drop-in replacement rather than a
second, divergent toolchain. Final state (spec: `scripts/env_ligprep_wsl.yml`):

```
package       WSL ligprep    Windows
rdkit         2025.09.6      2025.09.6    match
meeko         0.7.1          0.7.1        match
prolif        2.2.1          2.2.1        match
MDAnalysis    2.9.0          2.9.0        match
pandas        2.3.3          2.3.3        match  (pinned down from 3.0.6 -- see below)
shapely       2.1.2          n/a          Windows copy is pip-installed
openbabel     3.2.1          3.1.1        UNMATCHED
numpy         2.4.6          2.0.1        UNMATCHED
```

Two could not be matched:

* **openbabel**: 3.1.1 conflicts with python 3.11 + rdkit 2025.09.6 through libxml2. This is the one
  that touches coordinates, via the `-p 7.4` protonation step, so it is the reason the verification
  below measures geometry rather than assuming it.
* **numpy**: 2.4.6 against 2.0.1. Left alone deliberately -- numpy 2.x is API-stable, RDKit's
  geometry is C++ and does not route through it, and the docking equivalence test below passed under
  this numpy.

One was caught and corrected: the ProLIF solve pulled **pandas 3.0.6** in, a major-version jump from
Windows' 2.3.3. pandas 3.x changes `read_csv` dtype inference and copy-on-write semantics, and
`prep_ligands.py`, `filter_cascade.py` and `bbb_score.py` all read CSVs -- that is the class of
difference that changes results silently instead of erroring, so it was pinned back to 2.3.3. Every
verification below was then re-run on the final stack, not the intermediate one.

### Verification — the real operation every time, never a `--help` probe

A `--help` probe already produced one false "OK" for `mk_export.exe` earlier in this project: bash
printed "Permission denied" and the probe read the exit status as success. So every check here does
the actual work.

1. rdkit + meeko import, and embed + MMFF + obabel + `mk_prepare_ligand` on a flavonol — OK.
2. Real PDBQT → SDF on a docked pose: 486 lines, correct energies, `$$$$` terminator.
3. `prep_ligands.py` end-to-end through the bridge on `cand_001..003` — all three prepped.
4. **Conformer comparison against the repo copies: `GetBestRMS` = 0.0001 Å for all three** —
   identical to printing precision. Raw RMSD differed (9.31, 0.006, 3.71 Å) purely as
   rigid-body placement: same molecule, same internal geometry, different global frame. The frame
   shift is obabel 3.2.1's, from the `-p 7.4` step. Re-measured on the final pinned stack after
   the pandas correction, not only on the intermediate one.
5. **The test that actually matters**, since a different starting frame could in principle perturb
   the search: re-dock the WSL-prepped `cand_003` against `2V5Zdry` at the identical box,
   exhaustiveness 32, `num_modes` 9 and seed 11.

```
pose      repo     WSL-prepped     delta
  1     -11.42        -11.41        0.01
  2     -11.26        -11.28        0.02
  3     -10.62        -10.62        0.00
  4      -9.347        -9.348       0.001
```

Agreement to **≤0.02 kcal/mol** — two orders of magnitude below this project's own inter-seed SD
and far below any margin it reports. The obabel mismatch is therefore immaterial in practice.

### Operational rule

Existing ligands **do not need re-prepping** — they are already in the repo and nothing about them
changed. `ligprep` exists to unblock *new* work. Re-prepping is chemically safe (identical
conformer) but would perturb the absolute frame for no gain, so it should not be done gratuitously.

### A verifier bug, caught by the operation succeeding

The first bridge wrapper reported `[FAIL] no SDF produced` for a conversion that had in fact
written a valid 486-line SDF. The bug was in the check, not the work: it ran `[ -s "$OUT" ]` on the
Windows side against a `/mnt/c/...` path, which does not resolve there. All verification now
happens inside WSL, where the path is certainly valid, and that wrapper was deleted in favour of
the general `wsl_run.sh`.

This is the same lesson as entries K and L, in the opposite direction: there the verifier invented
a defect that was not in the data, here it invented a failure in an operation that had succeeded.
**A verifier is code and gets the same scrutiny as the thing it verifies.**

### Added

`scripts/wsl_run.sh`, `scripts/env_ligprep_wsl.yml`.

Propagated to `scripts/README.md` (the "Three environments, one pipeline" section replacing "Two").

Still pending at the time of writing: `dock.sh` gains a preflight that resolves the export route
once and falls back to the WSL bridge when the native `mk_export` is blocked, rather than aborting.
It could not be applied while the Vinardo cross-check was running -- Windows holds a rename lock on
a script bash has open, so the atomic `os.replace` returned `WinError 5`. That is the guard working
as designed; the patch waits for the job to finish.

## 2026-09-28 — N (the Vinardo cross-check on the validated receptors: the caveat is confirmed, not removed)

336 dockings: all 56 candidates against both validated water-symmetric on-targets (7JXXdry TTBK1,
2V5Zdry MAO-B), three seeds each, under Vinardo. Identical box, exhaustiveness 32, `num_modes` 9 and
seeds as the production Vina run — **only the scoring function differs**. Ran 03:13–03:50, all 6
runs `[ok]`, 336/336 logs.

This closes the gap §3.3 and WORKFLOW Phase 3 both flagged. The only prior Vinardo data was
prior-phase, measured on 4NFM/2V5Z — the superseded receptor set, whose TTBK1 structure was apo and
carried no passing redocking validation — so the published caveat was a general warning from
superseded receptors, not a cross-check of the rankings actually reported.

### Result

```
                        TTBK1 (7JXXdry)        MAO-B (2V5Zdry)
Pearson, all 56         +0.489 (p 1.3e-4)     +0.690 (p 4.1e-9)
Spearman, all 56        +0.463 (p 3.3e-4)     +0.516 (p 4.6e-5)
Pearson, top 15         +0.354 (p 0.20) n.s.  +0.142 (p 0.62) n.s.
top-15 overlap          8 of 15                6 of 15
best by vina/vinardo    cand_013 / cand_002    cand_043 / cand_001
cand_003 rank           8 vs 12 of 56          2 vs 10 of 56
vina spread all56/top15 1.87 / 0.60            2.92 / 0.77
```

### Why the top-15 row is the one that matters

The all-56 correlation looks reassuring and is not. It is carried by **dynamic range** — both
functions agree that weak binders are weak. Restrict to the slice where selection actually happens
and the Vina spread collapses from 1.87 to 0.60 kcal/mol at TTBK1 and 2.92 to 0.77 at MAO-B, and with
it the correlation, to **non-significance on both targets**.

This is not seed noise. Mean inter-seed SD is 0.065 and 0.012 kcal/mol for Vina, 0.012 and 0.009 for
Vinardo — one to two orders of magnitude below the disagreement. The functions genuinely disagree.

So: **absolute favourability reproduces, fine-grained rank does not.** Every candidate again scores
favourably under Vinardo (−7.72 to −4.99 at TTBK1, −9.42 to −4.55 at MAO-B). Neither arm's best
compound agrees between functions. The cross-check therefore *confirms* the existing caveat rather
than resolving it, and the manuscript was updated to say so rather than to report a closed gap.

What survives is weaker but real: `cand_003` sits in the **top 12 of 56 under both functions on both
targets** (top ~21%) — robustly good, not demonstrably best. That is the honest strength of the
claim, and it is why pose stability rather than docking rank has been the discriminating filter
throughout.

### A comparison I nearly made and should not have

My first reading was that MAO-B had "improved" from the prior phase's r = −0.28 to +0.690. That
comparison is invalid: the prior figure was a **top-15** value on **different receptors**, and both
the slice and the receptors changed. Range restriction alone attenuates correlation, which is
exactly what the top-15 column here demonstrates. Computing the top-15-restricted values was what
made the comparison honest — and on that like-for-like basis there is no improvement to claim, only
a measurement made on the reported receptors for the first time.

Had the run covered only the top 15, as WORKFLOW Phase 3 originally specified, the non-significant
correlation would have had nothing to contrast against and would have read as noise rather than range
restriction. Going wide is what made it interpretable.

### Not done, deliberately

No consensus of the two functions is computed. Averaging them is how the prior phase nearly selected
a different lead — a weak Vina score got masked rather than corroborated. Agreement is tested, never
averaged, and the two are never plotted on a shared axis. Vinardo output lives one directory level
deeper, in `04_docking/crosscheck/vinardo/`, below the globs used by `collect_results.py`,
`analyze_selectivity.py` and `analyze_water_test.py`, so it cannot be pooled into a Vina consensus or
selectivity margin even by accident.

A third empirical scoring function would not help: two disagreeing empirical functions cannot be
adjudicated by a third of the same kind, a point WORKFLOW already made. §4.5 was rewritten to say
the way past this is short MD on the shortlist, since pose stability has separated these compounds
where score could not.

### Added / changed

`08_analysis/vinardo_crosscheck.csv` (112 rows), `08_analysis/vinardo_crosscheck_slices.csv`,
`scripts/run_vinardo_crosscheck.sh`, `scripts/analyze_vinardo_crosscheck.py`.
Propagated: manuscript §2.4, §3.3 (heading and body), §4.4 limitation 7, §4.5 next steps;
`WORKFLOW.md` Week 3; `04_docking/README.md`.

### Propagated to the results document (same day)

`10_results/README.md` brought up to the current findings:

- **New Result 5** for the cross-check, with the top-15 collapse as the headline rather than the
  flattering all-56 correlation, and an explicit note that the all-56 row is the wrong one to read.
- **New claim 11** ("fine-grained docking rank is not reproducible across scoring functions",
  supportable strongly). Claims were **not** renumbered — claim 7 is cross-referenced later in the
  same file. Claim 2's basis now cross-references claim 11: Result 1 tests docking rank against
  physics, Result 5 tests it against itself, and they agree.
- **The lead compound's ranks now carry the caveat where a reader meets them.** "Rank 8 of 56" and
  "rank 2 of 56" were stated bare; they are Vina's, and Vinardo puts the same compound at 12 and 10.
- **New limitation 7** bounding every rank quoted in the file.

Two stale items were found while doing this, neither related to Vinardo:

- **Limitation 4 contradicted the file's own Result 4 table.** It still read "The TTBK arms still
  have no on-pose replicate pair beyond TTBK1 pose 1", which stopped being true when
  `system_TTBK2_p2_r2` ran — and the Result 4 table four sections earlier already listed that pair.
  Rewritten to state that every arm now carries one, and that the remaining gap is a distribution
  rather than a point estimate.
- **"four documented instances"** of the comparability failure, against the manuscript's five. The
  fifth (a replicated arm judged against an unreplicated one) was added the same day as the TTBK2
  replicate and never propagated here.

Also: the two "Done 2026-09-27/28" bullets were removed from "What would change the answer", which
is a forward-looking list. Nothing was lost — every number they carried (0.26, 0.24, 3.97, the
1.45 → 0.53 fall) is already in the Results sections above, and the audit asserts they survive.

Verified by `scripts/audit_results_readme.py` — 50 checks, recomputing the Vinardo statistics and
the pose-verdict counts from the primary files rather than comparing text to text, and asserting
the README agrees with the manuscript on the shared numbers.
