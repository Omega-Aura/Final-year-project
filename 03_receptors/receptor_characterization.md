# Receptor Characterization — Six Docking Targets

All values below were extracted directly from the downloaded coordinate files in
`03_receptors/<PDB>/raw.pdb` (RCSB originals) and from the prepared
`clean_noH.pdb` / `receptor.pdbqt` files, not copied from literature. Generated 2026-08-28.

---

## 1. Protein IDs and identity

| PDB ID | Protein | UniProt | Entry name | EC / class | Title in file |
|---|---|---|---|---|---|
| **7JXX** | TTBK1 (tau-tubulin kinase 1), kinase domain | Q5TCY1 | TTBK1_HUMAN | Transferase (kinase) | TTBK1 kinase domain in complex with compound 3 |
| **4BTK** | TTBK1, kinase domain | Q5TCY1 | TTBK1_HUMAN | Transferase (kinase) | TTBK1 in complex with inhibitor |
| **7Q8V** | TTBK1, kinase domain | Q5TCY1 | TTBK1_HUMAN | Transferase (kinase) | TTBK1 in complex with VNG2.73 (compound 42) |
| **7Q8Y** | TTBK2, kinase domain | Q6IQ55 | TTBK2_HUMAN | Transferase (kinase) | TTBK2 in complex with VNG2.73 (compound 42) |
| **2V5Z** | MAO-B (monoamine oxidase B) | P27338 | AOFB_HUMAN | Oxidoreductase, EC 1.4.3.4 | Human MAO-B with the selective inhibitor safinamide |
| **2Z5X** | MAO-A (monoamine oxidase A) | P21397 | AOFA_HUMAN | Oxidoreductase, EC 1.4.3.4 | Human monoamine oxidase A with harmine |

All six are **X-ray diffraction** structures of **human** proteins.

### Residue-numbering caveat (matters when comparing pocket residues)

| Structure | PDB numbering vs UniProt (from DBREF) | Offset |
|---|---|---|
| 7JXX | A 15–343 → UNP 15–343 | **0** (native UniProt numbering) |
| 7Q8V | A 13–320 → UNP 13–320 | **0** |
| **4BTK** | A 25–337 → UNP **1–313** | **+24** — 4BTK numbers are 24 higher than 7JXX/7Q8V |
| 7Q8Y | A/B 1–299 → UNP 1–299 | 0 (but TTBK2 sequence; equivalent residues sit ~13 lower than TTBK1) |
| 2V5Z | A/B 1–520 → UNP 1–520 | 0 |
| 2Z5X | A 12–524 → UNP 12–524 | 0 |

Worked example: the hinge glutamine is **GLN110 in 7JXX/7Q8V**, **GLN134 in 4BTK**
(110 + 24), and **GLN97 in 7Q8Y** (TTBK2 equivalent). Any cross-structure pocket
table must state which numbering it uses.

---

## 2. Number of amino acids

"SEQRES" = residues in the deposited construct. "Observed" = residues with actual
coordinates (what docking/MD actually sees). "Missing" = REMARK 465 disordered residues.

| PDB | Chains in file | SEQRES per chain | Observed (chain A) | Missing (chain A) | Full-length protein (UniProt) | Construct coverage |
|---|---|---|---|---|---|---|
| 7JXX | A | 332 | **288** | 44 | 1321 aa | kinase domain only (15–343) |
| 4BTK | A | 337 | **287** | 50 | 1321 aa | kinase domain only (UNP 1–313) |
| 7Q8V | A | 309 | **292** | 17 | 1321 aa | kinase domain only (13–320) |
| 7Q8Y | A, B | 300 each | **294** (B: 293) | 6 (B: 7) | 1244 aa | kinase domain only (1–299) |
| 2V5Z | A, B | 520 each | **499** (B: 494) | 21 (B: 26) | 520 aa | full-length |
| 2Z5X | A | 513 | **513** | **0** | 527 aa | 12–524, no disorder |

Notes:
- The four kinase structures are **isolated catalytic domains**, roughly 20–25% of the
  full-length protein. Fine for ATP-site docking, but the N-terminal extension and the
  long C-terminal tail of TTBK1/2 are absent.
- **2Z5X is the only structure with zero missing residues** — every SEQRES residue is
  modeled. It is the cleanest input of the six.
- 2V5Z chain A is used throughout the project; chain B is a crystallographic copy and
  is discarded by `prep_receptor.sh`.
- Prepared receptor sizes (`receptor.pdbqt`, chain A + retained cofactor, polar H added):
  7JXX 2803 atoms, 4BTK 2902, 7Q8V 2943, 7Q8Y 2950, 2V5Z 4898, 2Z5X 5007.
- Construct lengths and modeled ranges are taken from each file's DBREF record; the
  full-length column is the UniProt reference length for that entry.

### Incomplete side chains (REMARK 470)

| PDB | Residues with missing side-chain atoms |
|---|---|
| 7JXX | 9 — ASP22, PHE45, LYS75, ARG98, GLU100, LYS101, ASN195, ARG200, LYS212 |
| 7Q8Y | 2 — GLU6 (chain A), PHE32 (chain B) |
| 2V5Z | 1 — ILE501 (chain A) |
| 4BTK, 7Q8V, 2Z5X | none |

None of these sit in a docking pocket.

---

## 3. Resolution and refinement quality

| PDB | Method | Resolution | R-work | R-free | Assessment |
|---|---|---|---|---|---|
| 7JXX | X-ray | **1.56 Å** | 0.180 | 0.213 | Best of the set — high resolution, tight R-free gap |
| 7Q8Y | X-ray | **1.60 Å** | 0.178 | 0.205 | Excellent |
| 2V5Z | X-ray | **1.60 Å** | 0.208 | 0.227 | Excellent |
| 4BTK | X-ray | **2.00 Å** | 0.190 | 0.249 | Good; largest R-work/R-free gap (0.059) |
| 7Q8V | X-ray | **2.13 Å** | 0.187 | 0.244 | Acceptable |
| 2Z5X | X-ray | **2.20 Å** | 0.201 | 0.255 | Lowest resolution of the set |

All six are at or below 2.2 Å, comfortably within normal practice for structure-based
docking. The set spans a 0.64 Å resolution range, so side-chain placement in 2Z5X and
7Q8V is inherently less certain than in 7JXX/7Q8Y/2V5Z — worth one sentence in
Limitations.

---

## 4. Types of bonds present

### 4a. Covalent bonds — verified from coordinates, not just from header records

| Bond type | Present in | Detail |
|---|---|---|
| **Peptide (amide) bonds** | all six | Backbone; count = observed residues − number of contiguous chain segments. Chain breaks occur at every REMARK 465 gap (see §5c). |
| **Disulfide bonds** | **none in any structure** | SSBOND records are absent in all six, and a direct geometric scan of every CYS SG–SG pair (< 2.5 Å) found **zero** genuine disulfides. Free cysteines: 7JXX 3, 4BTK 3, 7Q8V 3, 7Q8Y 10, 2V5Z 18, 2Z5X 9. The scan initially flagged a 1.45 Å "pair" in 7JXX — that is CYS172 modeled in two alternate conformations (altloc A and B), not a disulfide. `mk_prepare_receptor` was run with `--default_altloc A`, so only conformer A reaches the pdbqt. |
| **Cis-peptide bonds** | 2V5Z, 2Z5X only | 2V5Z: ASN275–PRO276 and **CYS397–TYR398** (both chains). 2Z5X: ARG284–PRO285 and **CYS406–TYR407**. The Cys–Tyr cis bonds are the conserved strained geometry adjacent to the flavin attachment site — a real MAO-family feature, not a modeling error. |
| **Covalent flavinylation (8α-S-cysteinyl thioether)** | 2V5Z, 2Z5X | Measured SG···C8M distance: **2Z5X CYS406–FAD = 1.65 Å** (a proper C–S covalent bond, and flagged in that file's own REMARK 500) and **2V5Z CYS397–FAD = 2.31 Å in both chains** — long for a thioether; the link is real biologically but was refined with loose restraints here. Neither file carries a LINK record for it, so **no docking or MD tool will treat FAD as covalently attached unless the bond is added manually.** |
| **Phosphoanhydride / phosphodiester** | 2V5Z, 2Z5X (FAD); 7Q8V, 7Q8Y (PO4) | FAD contains the adenosine–pyrophosphate–ribityl backbone (2 P atoms per FAD). PO4 ions in 7Q8V (1) and 7Q8Y (4) are crystallization-buffer ions. |
| **Glycosidic C–N bond** | 2V5Z, 2Z5X | FAD adenine N9–ribose C1′. |

### 4b. Coordination / ionic bonds

Only **7JXX** contains a metal: 3 Na⁺ ions (7 LINK records). The one nearest the pocket,
NA 404, is coordinated by **GLN110 OE1 (2.72 Å), THR168 OG1 (2.72 Å), ARG164 N (2.90 Å)**
plus waters. It is almost certainly a crystallization artifact. `prep_receptor.sh` was
called with no cofactor list for 7JXX, so **all three Na⁺ were stripped from
`receptor.pdbqt`** — the right call, but note that GLN110 is the hinge residue, so its
side-chain rotamer in this structure was partly shaped by an ion that is no longer there.

### 4c. Non-covalent interaction types in each pocket (from the native ligand)

Measured directly from `raw.pdb`: polar contacts ≤ 3.5 Å between N/O/S atom pairs,
hydrophobic C···C ≤ 4.5 Å, halogen/F contacts ≤ 4.0 Å.

| PDB | Native ligand | Key polar / H-bonds (Å) | Halogen | Hydrophobic contact residues |
|---|---|---|---|---|
| 7JXX | VP7 | GLU77 OE1···O1 2.61; **GLN110 N···N1 2.88** (hinge); PHE177 N···O1 2.91; GLN108 O···N2 2.95 | — | 13: ILE40, ILE48, ALA61, LYS63, GLU77, LEU81, VAL105, MET107, LEU109, GLN110, LEU175, ASP176, PHE177 |
| 4BTK | DTQ | GLU101 OE1···O21 2.46; LYS87 NZ···O21 2.82; **GLN134 N···N3 2.90** (hinge) | — | 14: ILE64, GLY65, ILE72, ALA85, LYS87, GLU101, MET131, GLN132, LEU133, GLN134, GLY135, LEU199, ASP200, PHE201 |
| 7Q8V | 9IV | **GLN110 N···N2 2.86** (hinge); GLN108 O···N3 3.02 | none within 4.0 Å | 13: ILE40, GLY41, ILE48, ALA61, LYS63, MET107, LEU109, GLN110, ASN113, ALA115, ASP116, SER158, LEU175 |
| 7Q8Y | 9IV | **GLN97 N···N2 2.84** (hinge); GLN95 O···N3 2.86 | **ASP103 OD2···Cl 2.90** (halogen bond); ASN100 ND2···Cl 3.94 | 13: ILE27, GLY28, ILE35, ALA48, LYS50, MET94, GLN95, GLN97, ASN100, ALA102, ASP103, SER145, LEU162 |
| 2V5Z | SAG (safinamide) | **GLN206 OE1···N21 2.81**; GLN206 OE1···N16 3.30 | LEU164 O···F3 3.92 (weak F contact) | 14: PRO104, TRP119, LEU164, PHE168, LEU171, CYS172, ILE198, ILE199, GLN206, ILE316, TYR326, PHE343, TYR398, TYR435 |
| 2Z5X | HRM (harmine) | TYR407 OH···NAH 3.31; GLN215 NE2···NAI 3.44 | — | 12: TYR69, ILE180, ASN181, PHE208, GLN215, ILE325, ILE335, LEU337, MET350, PHE352, TYR407, TYR444 |

Pattern: **all four kinase sites are driven by a single backbone hinge H-bond to a
glutamine** (GLN110 / GLN134 / GLN97) plus a deep hydrophobic sleeve — classic type-I
ATP-competitive recognition. Both MAO sites are dominated by an **aromatic cage**
(TYR398/TYR435 in MAO-B, TYR407/TYR444 in MAO-A — the flavin-flanking tyrosine pair)
with only weak polar anchoring, so MAO scoring depends far more on shape complementarity
and π-stacking than on hydrogen bonding.

---

## 5. Missing content — inhibitors, cofactors, and chain gaps

### 5a. Co-crystallized ligands present in each file

| PDB | Native inhibitor | Chemical name (HETNAM) | Cofactors / other heteroatoms |
|---|---|---|---|
| 7JXX | **VP7** ×1 | 4-(2-amino-5,6,7,8-tetrahydropyrimido[4′,5′:3,4]cyclohepta[1,2-b]indol-11-yl)-2-methylbut-3-yn-2-ol | 3 × Na⁺, waters |
| 4BTK | **DTQ** ×1 | 4-[3-hydroxyanilino]-6,7-dimethoxyquinazoline | DMSO ×1 |
| 7Q8V | **9IV** ×1 | N-[4-(2-chlorophenoxy)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine | 1 × PO4, waters |
| 7Q8Y | **9IV** ×2 (one per chain) | same as 7Q8V | 4 × PO4, waters |
| 2V5Z | **SAG** ×2 (safinamide) | (S)-(+)-2-[4-(fluorobenzyloxy-benzylamino)propionamide] | **FAD ×2**, waters |
| 2Z5X | **HRM** ×1 (harmine) | 7-methoxy-1-methyl-9H-β-carboline | **FAD ×1**, DCX ×2 (decyl-dimethyl-phosphine oxide, a detergent), GOL ×3 (glycerol), waters |

**No inhibitor is missing from any deposited file** — all six are inhibitor-bound holo
complexes, which is exactly what is wanted for docking. Each native ligand was correctly
extracted to `native_<LIG>.pdb` and used to define the docking box.

### 5b. FAD is missing from the prepared MAO-A receptor

This is the one real defect found, and it should be fixed before any MAO-A result is
reported.

| File | 2V5Z (MAO-B) | 2Z5X (MAO-A) |
|---|---|---|
| `raw.pdb` FAD | present (53 atoms) | present (53 atoms) |
| `clean_noH.pdb` FAD | present (53 atoms) | **present (53 atoms)** |
| `receptor.pdbqt` FAD | **present** | **ABSENT — 0 atoms** |

The `awk` cofactor filter in `scripts/prep_receptor.sh` worked correctly for both files;
the FAD was lost at the **Meeko `mk_prepare_receptor` step**, for 2Z5X only.

Why it matters — the flavin isoalloxazine ring forms the back wall of the MAO substrate
cavity. Checking the actual `box.json` grids:

- **2Z5X box contains 13 FAD atoms** (N1, C2, O2, N3, C4, O4, C4X, N5, C5X, C6, C9A, N10,
  C10 — the full isoalloxazine face); nearest FAD atom 6.98 Å from box centre.
- 2V5Z box contains the 13 equivalent atoms, nearest 9.22 Å from box centre — and those
  atoms **are** in its pdbqt.

So every MAO-A docking run so far searched a pocket with a ~13-atom hole where the flavin
should be, while MAO-B searched the correct closed pocket. Ligands can drift into the
vacant flavin space, and MAO-A scores are therefore **not directly comparable to MAO-B
scores** — including the reference-set selectivity comparison already logged (harmine
2Z5X −8.66 vs 2V5Z −8.62; a 0.04 kcal/mol margin is well inside the error this
introduces).

**Action:** re-run receptor prep for 2Z5X so FAD survives into the pdbqt, then re-dock
everything scored against MAO-A. The call is `scripts/prep_receptor.sh 2Z5X A HRM "FAD"`;
if Meeko still drops the flavin, supply FAD as an explicit residue template or merge the
FAD block into the pdbqt manually, as was evidently done for 2V5Z.

### 5c. Chain gaps near the binding sites

Disordered stretches close to a docking box can leave the pocket partly open-walled.
Distances are from the nearest flanking observed Cα to the box centre in `box.json`.

| PDB | Gap | Length | Nearest flanking Cα → box centre | Concern |
|---|---|---|---|---|
| 7JXX | 12–21 | 10 | 16.3 Å | no |
| 7JXX | **44** | 1 | **11.8 Å** | glycine-rich loop edge, adjacent to pocket residues ILE40/ILE48 — minor |
| 7JXX | 187–188 | 2 | 28.2 Å | no |
| 7JXX | 313–343 | 31 | 25.9 Å | no (C-terminal tail) |
| 4BTK | 1–45 | 45 | 16.7 Å | no (expression tag + N-terminus) |
| 4BTK | **67–69** | 3 | **10.0 Å** | **in the glycine-rich loop, immediately after pocket residues ILE64/GLY65 — the P-loop roof over the ATP site is incomplete** |
| 4BTK | 336–337 | 2 | 21.6 Å | no |
| 7Q8V | 12–20, 313–320 | 9, 8 | 19.4, 23.9 Å | no |
| 7Q8Y | 0–5 | 6 | 22.1 Å | no |
| 2V5Z | 1–2, 502–520 | 2, 19 | 43.9, 44.1 Å | no |
| 2Z5X | — | — | — | **none — fully modeled** |

**4BTK is the weakest structure of the set**: 2.00 Å, the largest R-work/R-free gap,
50 missing residues, and a 3-residue break in the glycine-rich loop directly over the ATP
pocket. Treat 4BTK scores as the least reliable of the three TTBK1 structures and prefer
7JXX (1.56 Å, intact P-loop) as the primary TTBK1 target.

---

## 6. Summary of issues found

| # | Issue | Severity | Affected | Fix |
|---|---|---|---|---|
| 1 | FAD absent from `2Z5X/receptor.pdbqt` although present in `clean_noH.pdb`; 13 flavin atoms fall inside the MAO-A docking box | **High** | all MAO-A docking and MAO-A/MAO-B selectivity results | re-prep 2Z5X with FAD retained, re-dock MAO-A |
| 2 | Cys→FAD covalent link has no LINK record in either MAO file; 2V5Z SG···C8M is 2.31 Å (long) | Medium | MD parameterization of MAO systems | add the bond explicitly if MAO MD is run |
| 3 | 4BTK P-loop gap 67–69 sits directly above the ATP site | Medium | 4BTK docking scores | de-prioritize 4BTK; use 7JXX as primary TTBK1 |
| 4 | 4BTK residue numbering is offset +24 from UniProt / 7JXX / 7Q8V | Medium | any cross-structure pocket comparison | state the numbering convention explicitly in Methods |
| 5 | 7JXX Na⁺ ions stripped; one of them coordinated the hinge residue GLN110 | Low | 7JXX pocket geometry | note in Limitations; no action needed |
| 6 | 7JXX CYS172 has A/B alternate conformers | Low | already resolved by `--default_altloc A` | none |
| 7 | Resolution spread 1.56–2.20 Å across the set | Low | cross-target score comparison | one sentence in Limitations |

**No structure is missing its inhibitor.** The only genuinely missing molecule anywhere
in the pipeline is the **MAO-A FAD cofactor**, and only in the prepared pdbqt — the
downloaded structure has it.
