# cand_003 redocked against MAO-A / MAO-B / TTBK1 / TTBK2 (all heteroatoms removed)

Open `*_best_pose_complex.pdb` in Discovery Studio (protein + hydrogens + ligand C03, chain L, with CONECT).
`*_all9_poses.pdb` = the nine modes as MODELs. `cand_003_redock_report.pdf` = summary and figures
(a PDF cannot be opened as a structure). `all_modes_scores.csv` = 4 receptors x 3 seeds x 9 modes.

| receptor | best (kcal/mol) |
|---|---|
| MAO-B 2V5Z | -11.15 |
| MAO-A 2Z5X | -8.77 (FAD-clash-free; raw -11.44 sits inside the deleted FAD) |
| TTBK1 7JXX | -8.49 |
| TTBK2 7Q8Y | -9.87 (pose 3.5 A off crystal mode) |

Vina 1.2.7, known-site box, exhaustiveness 32, seeds 11/22/33. Read the PDF before quoting the MAO-A or TTBK2 numbers.
Provenance: LOGBOOK 2026-09-29 entry O. Prepared receptors and logs: `04_docking/cand003_redock/`.
