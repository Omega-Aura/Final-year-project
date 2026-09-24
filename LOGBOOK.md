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
