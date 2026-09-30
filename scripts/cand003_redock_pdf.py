#!/usr/bin/env python
"""Build the cand_003 clean-receptor redock report (PDF). Runs in WSL `ligprep` (reportlab).
Run from project root after cand003_redock_export.py and the PyMOL renders:
  bash scripts/wsl_run.sh python scripts/cand003_redock_pdf.py
"""
import json, re
import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle)

D = "10_results/cand_003_redock_clean"
S = {s["receptor"]: s for s in json.load(open(f"{D}/summary.json"))}
OV = json.load(open(f"{D}/overlap_with_removed_atoms.json"))
ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontSize=17, spaceAfter=6)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontSize=12.5, spaceBefore=8, spaceAfter=4)
B = ParagraphStyle("B", parent=ss["BodyText"], fontSize=9.2, leading=12.4)
SM = ParagraphStyle("SM", parent=B, fontSize=8, leading=10.4, textColor=colors.HexColor("#444444"))
CELL = ParagraphStyle("CELL", parent=B, fontSize=8.3, leading=10.4)
P = lambda t, st=B: Paragraph(t, st)


def table(rows, widths, head=True):
    t = Table([[P(str(c), CELL) for c in r] for r in rows], colWidths=widths, repeatRows=1 if head else 0)
    st = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#bbbbbb")), ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
    if head:
        st.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e6eef7")))
    t.setStyle(TableStyle(st))
    return t


def contacts_pdb(complex_pdb):
    lig = [(l[76:78].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
           for l in open(complex_pdb) if l.startswith("HETATM") and l[76:78].strip() != "H"]
    res, hb = {}, []
    for l in open(complex_pdb):
        if not l.startswith("ATOM") or l[76:78].strip() == "H":
            continue
        p = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]); el = l[76:78].strip()
        for i, (s, q) in enumerate(lig):
            d = float(np.linalg.norm(p - q)); k = f"{l[17:20]}{int(l[22:26])}"
            if d <= 4.0: res[k] = 1
            if d <= 3.3 and s in ("N", "O") and el in ("N", "O"):
                hb.append(f"{k}:{l[12:16].strip()} ... lig {s}{i + 1} ({d:.2f} A)")
    return sorted(res, key=lambda k: int(re.sub(r"\D", "", k))), hb


def cvals(R):
    return "/".join(f"{S[R]['seeds'][k][0][1]:.2f}" for k in sorted(S[R]["seeds"]))


doc = SimpleDocTemplate(f"{D}/cand_003_redock_report.pdf", pagesize=A4, leftMargin=1.8 * cm, rightMargin=1.8 * cm,
                        topMargin=1.6 * cm, bottomMargin=1.5 * cm, title="cand_003 redocking, clean receptors",
                        author="Final year project")
E = []
E += [P("cand_003 redocked against MAO-A, MAO-B, TTBK1, TTBK2", H1),
      P("Clean-receptor protocol (no water, no heteroatoms, missing atoms rebuilt, polar H, Gasteiger charges) "
        "&middot; AutoDock Vina 1.2.7 &middot; known active site, exhaustiveness 32 &middot; 3 seeds", SM), Spacer(1, 6)]

E.append(P("Result", H2))
mA = "-8.77 (FAD-clash-free) / -11.44 raw"
rows = [["Receptor", "Target", "Best score (kcal/mol)", "Seeds 11/22/33", "Pose vs crystal site", "Verdict"],
        ["2V5Z", "MAO-B", "<b>-11.15</b>", cvals("2V5Z"), "0.8 A from native SAG centroid", "Reliable: in the substrate cavity, no clash with FAD or waters"],
        ["2Z5X", "MAO-A", "<b>-8.77</b> (raw top mode -11.44 is an artefact)", cvals("2Z5X") + " raw; clash-free mode 6 = -8.77 in all 3 seeds",
         "raw pose 3.5 A off native and inside FAD", "Use -8.77. Raw modes 1-5 overlap the deleted FAD"],
        ["7JXX", "TTBK1", "<b>-8.49</b>", cvals("7JXX"), "1.5 A from native VP7 centroid", "Reliable, in the ATP site"],
        ["7Q8Y", "TTBK2", "<b>-9.87</b>", cvals("7Q8Y"), "3.5 A from native 9IV centroid", "Treat with caution: pose is off the crystal binding mode (see notes)"]]
E.append(table(rows, [1.7 * cm, 1.5 * cm, 3.2 * cm, 3.2 * cm, 3.4 * cm, 4.3 * cm]))
E.append(Spacer(1, 4))
E.append(P("Scores are Vina mode-1 affinities of the best of three seeds; the seeds agree within 0.03 kcal/mol for MAO-B, "
           "TTBK1 and TTBK2 and within 0.06 for the raw MAO-A top mode, so the search is converged and the differences "
           "between receptors are not sampling noise.", SM))

E.append(P("Two things you must know before using these numbers", H2))
E.append(P("<b>1. MAO-A: the top score is an artefact of removing FAD.</b> Modes 1-5 of every seed place 4-12 ligand heavy atoms "
           "within 2.5 A of where FAD sits in the crystal (closest approach 0.3-0.55 A); the CF3 end of the ligand lies inside the "
           "flavin ring system (see figure). With the cofactor deleted, nothing penalised that. The best mode that avoids FAD is "
           "<b>-8.77 kcal/mol</b>, and the project's earlier MAO-A run that kept FAD (no waters) gave -8.19 to -8.24, so -8.77 is "
           "the number to quote. Both poses are supplied.", B))
E.append(P("<b>2. Stripping every heteroatom is not neutral for these targets.</b> In the crystals FAD lines the MAO cavity and TTBK2's "
           "9IV inhibitor binds through a bridging water. In the same frame, 17 of the ligand's atoms in the TTBK2 pose overlap "
           "crystal waters (0 for MAO-B, 3 within 3.5 A for TTBK1), which is consistent with the pose settling 3.5 A away from the "
           "crystal ligand. The project's own dry-vs-water test found the same TTBK2 sensitivity. I followed your instruction as "
           "given; if you want a like-for-like comparison, the receptors with FAD and the pocket-water shell are already in "
           "03_receptors/ and can be docked the same way.", B))

E.append(P("Comparison with the project's earlier cand_003 runs (Vina, best of 3 seeds)", H2))
rows = [["Receptor", "This run (all heteroatoms removed)", "Earlier: no water, cofactor kept", "Earlier: pocket waters + cofactor kept"],
        ["MAO-B 2V5Z", "-11.15", "-11.42", "-12.10"],
        ["MAO-A 2Z5X", "-11.44 raw / -8.77 clash-free", "-8.24", "not run for candidates"],
        ["TTBK1 7JXX", "-8.49", "-8.41", "-9.22"],
        ["TTBK2 7Q8Y", "-9.87", "-9.89", "not run (needs bridging-water variant)"]]
E.append(table(rows, [3 * cm, 5 * cm, 4.6 * cm, 4.9 * cm]))
E.append(P("Removing waters costs roughly 0.7-0.8 kcal/mol on MAO-B and TTBK1, as expected. Only MAO-A moves in the other direction, "
           "and that is the FAD artefact above.", SM))

E.append(PageBreak())
E.append(P("Preparation and docking protocol (your steps 1-8)", H1))
steps = [["Step", "What was done"],
         ["1. Remove water", "All HOH removed (2V5Z 717, 2Z5X 182, 7JXX 291, 7Q8Y 739 in the raw file). Chain A only, first alternate location."],
         ["2. Remove heteroatoms", "2V5Z: FAD, SAG (native). 2Z5X: FAD, HRM (native), DCX, GOL. 7JXX: VP7 (native), Na. 7Q8Y: 9IV (native), PO4."],
         ["3. Missing atoms", "PDBFixer 1.x: 3 heavy atoms rebuilt in 2V5Z, 0 in 2Z5X, 38 in 7JXX, 4 in 7Q8Y (+ C-terminal OXT). Unresolved residues were "
          "<b>not</b> built (invented loop coordinates): 2V5Z N-term MET-SER and 19-residue C-terminal tail; 7JXX N-term 10 residues, 1 residue after 22, "
          "2 after 164, and a 31-residue C-terminal segment; 7Q8Y N-term 6 residues. 2Z5X is complete. The 7JXX gaps at 22 and 164 are chain breaks."],
         ["4. Polar hydrogens", "Hydrogens added at pH 7.4 (PDBFixer), then written to PDBQT by Meeko with non-polar H merged and polar H typed HD. "
          "Meeko was used rather than MGLTools, because MGLTools failed to bond PDBFixer's hydroxyl H (Ser/Thr/Tyr) and gave several Asn/Arg residues "
          "wrong charges."],
         ["5. Charges", "Gasteiger. Verified: every residue's charge is a whole number, and totals match the protonation states exactly: "
          "2V5Z -3, 2Z5X +2, 7JXX +12, 7Q8Y +11. Only chain termini are fractional. Ligand cand_003 sums to -0.001 (neutral, 26 atoms in the PDBQT, 2 polar H)."],
         ["6. Site and exhaustiveness", "Active site is known, so this is <b>not</b> blind docking. Box = the project's native-ligand box (18 A cube; 18 x 18.9 x 18 for 2V5Z), "
          "centred on the crystal ligand. Exhaustiveness 32 (not 100), 9 modes, seeds 11/22/33."],
         ["7. Results", "This report, all_modes_scores.csv (108 rows), and per-receptor complex PDBs."],
         ["8. Software", "AutoDock Vina 1.2.7, Vina scoring function, on AutoDock-format PDBQT files. AutoDock4 (ad4) was <b>not</b> used."]]
E.append(table(steps, [3.6 * cm, 13.7 * cm]))
E.append(Spacer(1, 6))
E.append(P("Files for Discovery Studio", H2))
E.append(P("Discovery Studio does not open PDF files as structures, so the structures are supplied as <b>PDB</b>: "
           "*_best_pose_complex.pdb (protein with hydrogens + docked ligand, residue C03, chain L, with full CONECT bonding) opens directly; "
           "*_all9_poses.pdb steps through all nine modes as MODELs; *_best_pose_ligand.pdb / .mol are the ligand alone. This PDF is the "
           "summary and figures. Ligand hydrogens are explicit; protein hydrogens are those added at pH 7.4.", B))

for R, title in (("2V5Z", "MAO-B (2V5Z)"), ("2Z5X", "MAO-A (2Z5X), raw top mode, overlaps FAD"),
                 ("2Z5X_clashfree", "MAO-A (2Z5X), best FAD-clash-free pose"), ("7JXX", "TTBK1 (7JXX)"), ("7Q8Y", "TTBK2 (7Q8Y)")):
    E.append(PageBreak())
    E.append(P(f"cand_003 in {title}", H1))
    base = R.split("_")[0]
    if R == "2Z5X_clashfree":
        cpx = f"{D}/2Z5X_MAO-A_cand_003_FADclashfree_pose_complex.pdb"; sc = "-8.77 (seed 11, mode 6)"
    else:
        cpx = f"{D}/{base}_{S[base]['target']}_cand_003_best_pose_complex.pdb"; sc = f"{S[base]['mode1']:.2f} (seed {S[base]['best_seed']}, mode 1)"
    res, hb = contacts_pdb(cpx)
    E.append(P(f"Vina score <b>{sc} kcal/mol</b>. Orange: cand_003. Blue: receptor residues within 4 A. Green: crystal ligand (removed before docking). "
               "Purple: crystal FAD (removed before docking). Black dashes: N/O...N/O pairs within 3.3 A.", SM))
    E.append(Image(f"{D}/img/{R}.png", width=17.2 * cm, height=17.2 * cm * 1000 / 1400))
    E.append(Spacer(1, 4))
    E.append(P("<b>Residues within 4 A:</b> " + " ".join(res), B))
    E.append(P("<b>Hydrogen-bond candidates (heavy-atom N/O...N/O &le; 3.3 A):</b> " + ("; ".join(hb) if hb else "none"), B))
    if base in OV:
        o = OV[base]
        bits = [f"{g}: closest {v['min_dist']} A" for g, v in o.items() if g not in ("water",) and g != {"2V5Z": "SAG", "2Z5X": "HRM", "7JXX": "VP7", "7Q8Y": "9IV"}[base]]
        E.append(P("<b>Distance to removed crystal heteroatoms (raw-pose):</b> " + "; ".join(bits) +
                   f"; waters: {o['water']['lig_atoms_within_3p5']} ligand atoms within 3.5 A.", SM))

E.append(PageBreak())
E.append(P("All modes, best seed per receptor", H1))
rows = [["Receptor", "Mode", "Affinity"]]
import csv
allrows = list(csv.DictReader(open(f"{D}/all_modes_scores.csv")))
cols = []
for R in ("2V5Z", "2Z5X", "7JXX", "7Q8Y"):
    bs = S[R]["best_seed"]
    m = [(int(r["mode"]), float(r["affinity_kcal_mol"])) for r in allrows if r["receptor"] == R and int(r["seed"]) == bs]
    cols.append((R, bs, m))
hdr = ["Mode"] + [f"{R} (seed {bs})" for R, bs, _ in cols]
rows = [hdr]
for i in range(9):
    rows.append([str(i + 1)] + [f"{m[i][1]:.2f}" if i < len(m) else "-" for _, _, m in cols])
E.append(table(rows, [2 * cm, 3.6 * cm, 3.6 * cm, 3.6 * cm, 3.6 * cm]))
E.append(Spacer(1, 6))
E.append(P("Full 108-row table (4 receptors x 3 seeds x 9 modes): all_modes_scores.csv. Overlap analysis: overlap_with_removed_atoms.json and "
           "mao_a_fad_clash_scan.py in the same folder.", SM))
doc.build(E)
print("wrote", f"{D}/cand_003_redock_report.pdf")
