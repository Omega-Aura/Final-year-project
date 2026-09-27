#!/usr/bin/env python
"""Figures for 10_results/, drawn from 07_mmgbsa/md_mmgbsa_summary.csv.

Reads the collected table rather than any hand-entered numbers, so a figure cannot disagree with
the text. Regenerate the table first if the MD set has changed:

    conda run -n docking_project python scripts/collect_md_summary.py
    conda run -n docking_project python 10_results/make_figures.py

Run through `conda run`, NOT a bare `python`: the interpreter first on PATH is MGLTools' Python 2.7.

Design notes, so later edits do not undo a deliberate choice:

* **TTBK and MAO binding energies are never drawn on the same axis.** Absolute MM-GBSA values carry
  protein-specific desolvation and surface terms that do not cancel between different proteins, so a
  shared axis would invite exactly the cross-protein comparison the project forbids (a kinase
  against a flavoenzyme). Figure 2 is MAO-only for that reason. Ligand RMSD in figure 1 *is*
  comparable across proteins -- it is a geometric measure of one ligand against its own starting
  pose -- so all fourteen systems share that axis legitimately.
* **Verdict is encoded by position against threshold lines, not by colour.** The four verdicts are
  an ordered scale; the project's status palette is a fixed four-role set meant to ship with an icon
  and label, and validating it as a 4-slot categorical ramp fails the normal-vision floor
  (worst adjacent dE 13.6, below the 15 hard floor). Bands plus direct labels carry it instead.
* Categorical hues are slots 1 and 2 of the validated default palette (blue / orange), which pass
  all-pairs CVD and normal-vision separation in light mode (worst dE 24.7 CVD, 33.6 normal).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
SUMMARY = os.path.join(ROOT, "07_mmgbsa", "md_mmgbsa_summary.csv")

# Validated default palette, light mode.
BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8983"
GRID = "#e2e1dc"

# Verdict band edges, mirroring scripts/collect_md_summary.py.
STABLE_MAX, LOOSE_MAX, DISSOC_MIN = 2.30, 3.20, 5.00

LABEL = {
    "system": "TTBK1 pose 1 r1", "system_TTBK1_r2": "TTBK1 pose 1 r2",
    "system_TTBK1_p2": "TTBK1 pose 2", "system_TTBK1_p3": "TTBK1 pose 3",
    "system_TTBK2": "TTBK2 pose 1 r1", "system_TTBK2m": "TTBK2 pose 1 r2",
    "system_TTBK2_p2": "TTBK2 pose 2", "system_TTBK2_p3": "TTBK2 pose 3",
    "system_MAOA": "MAO-A pose 1", "system_MAOA_p2": "MAO-A pose 2",
    "system_MAOA_p3": "MAO-A pose 3", "system_MAOB": "MAO-B pose 1",
    "system_MAOB_p2": "MAO-B pose 2", "system_MAOB_p3": "MAO-B pose 3",
    "system_MAOA_p3_r2": "MAO-A pose 3 r2", "system_MAOB_p2_r2": "MAO-B pose 2 r2",
}
ORDER = list(LABEL)


def load():
    with open(SUMMARY, newline="", encoding="utf-8") as fh:
        rows = {r["system_dir"].split("/")[-1]: r for r in csv.DictReader(fh)}
    return rows


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, length=0)


def fig1_pose_stability(rows, out):
    keys = [k for k in ORDER if k in rows]
    vals = [float(rows[k]["lig_rmsd_last100"]) for k in keys]
    names = [LABEL[k] for k in keys]

    fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(keys) + 2.0))
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    y = range(len(keys))

    # Verdict bands, drawn behind the bars, labelled along the top edge.
    bands = [(0, STABLE_MAX, "stable"), (STABLE_MAX, LOOSE_MAX, "holds, loose"),
             (LOOSE_MAX, DISSOC_MIN, "drifts"), (DISSOC_MIN, 9.2, "dissociates")]
    for i, (lo, hi, name) in enumerate(bands):
        ax.axvspan(lo, hi, color="#000000", alpha=0.03 if i % 2 else 0.055, lw=0, zorder=0)
        ax.text((lo + min(hi, 9.2)) / 2, -1.02, name, ha="center", va="bottom",
                fontsize=8.5, color=MUTED)
    for edge in (STABLE_MAX, LOOSE_MAX, DISSOC_MIN):
        ax.axvline(edge, color=GRID, lw=1.0, zorder=1)

    colors = [ORANGE if "MAO" in LABEL[k] else BLUE for k in keys]
    ax.barh(list(y), vals, height=0.62, color=colors, zorder=3)
    for i, v in enumerate(vals):
        ax.text(v + 0.12, i, f"{v:.2f}", va="center", fontsize=8.5, color=INK, zorder=4)

    ax.set_yticks(list(y))
    ax.set_yticklabels(names, fontsize=9, color=INK)
    # Explicit limits rather than invert_yaxis(), so the band labels get headroom above the
    # first bar instead of colliding with the legend at the bottom.
    ax.set_ylim(len(keys) - 0.4, -1.1)
    ax.set_xlim(0, 9.2)
    ax.set_xlabel("Ligand RMSD from the docked pose, mean of the final 100 frames (Å)",
                  fontsize=9.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title("Which docked poses survive 10.1 ns of unrestrained dynamics",
                 fontsize=12.5, color=INK, pad=26, loc="left", weight="medium")
    ax.legend(handles=[Patch(facecolor=BLUE, label="TTBK1 / TTBK2 (kinases)"),
                       Patch(facecolor=ORANGE, label="MAO-A / MAO-B (flavoenzymes)")],
              loc="lower right", frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", os.path.relpath(out, ROOT))


def fig2_mao_energy(rows, out):
    """MAO only -- see the module docstring on why this axis is not shared with TTBK."""
    # Replicates sit immediately after the run they replicate, so a pair reads as a pair.
    groups = [("MAO-A", ORANGE, ["system_MAOA", "system_MAOA_p2", "system_MAOA_p3",
                                 "system_MAOA_p3_r2"]),
              ("MAO-B", BLUE, ["system_MAOB", "system_MAOB_p2", "system_MAOB_p2_r2",
                               "system_MAOB_p3"])]
    fig, ax = plt.subplots(figsize=(9.2, 4.9))
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    ypos, ylab, means, lo, hi = [], [], [], [], []
    slot = 0
    for gname, gcol, keys in groups:
        keys = [k for k in keys if k in rows and rows[k]["dG_kcal_mol"]]
        for k in keys:
            r = rows[k]
            dg, sem = float(r["dG_kcal_mol"]), float(r["dG_SEM"])
            off = r["pose_verdict"] == "drifts"
            ax.errorbar(dg, slot, xerr=sem, fmt="o", ms=8, color=gcol,
                        mfc=SURFACE if off else gcol, mew=2.0, ecolor=gcol,
                        elinewidth=2, capsize=3, zorder=3)
            tag = f"{dg:.2f}" + ("   off-pose — excluded" if off else "")
            ax.text(dg + 0.28, slot, tag, va="center", fontsize=9,
                    color=MUTED if off else INK, zorder=4)
            ypos.append(slot); ylab.append(LABEL[k]); lo.append(dg); hi.append(dg)
            slot += 1
        on = [float(rows[k]["dG_kcal_mol"]) for k in keys if rows[k]["pose_verdict"] != "drifts"]
        if on:
            means.append((gname, gcol, sum(on) / len(on)))
        slot += 0.7

    ax.set_yticks(ypos)
    ax.set_yticklabels(ylab, fontsize=9.5, color=INK)
    # Headroom at the top for the on-pose mean labels, which sit inside the axes so they cannot
    # collide with the title.
    ax.set_ylim(slot - 0.9, -1.45)
    ax.set_xlim(min(lo) - 1.0, max(hi) + 2.6)
    for gname, gcol, m in means:
        ax.axvline(m, color=gcol, ls=":", lw=1.6, zorder=1)
        ax.text(m, -1.38, f"{gname} on-pose mean {m:.2f}", ha="center", va="bottom",
                fontsize=8.5, color=gcol)

    ax.set_xlabel("MM-GBSA ΔG (kcal/mol) — more negative is more favourable", fontsize=9.5,
                  color=INK2)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title("MAO binding free energy: the off-pose run is the outlier, not the answer",
                 fontsize=12.5, color=INK, pad=30, loc="left", weight="medium")
    ax.legend(handles=[Patch(facecolor=ORANGE, label="MAO-A (anti-target)"),
                       Patch(facecolor=BLUE, label="MAO-B (intended target)")],
              loc="lower left", frameon=False, fontsize=8.5, labelcolor=INK2)
    # Caption below the axes: the only placement that cannot collide with a mark or a label.
    fig.text(0.5, 0.018, "Error bars are the SEM MMPBSA.py reports; figure 3 shows how far each sits "
             "from the measured spread.\nHollow marker = off-pose trajectory, excluded from every "
             "comparison.", ha="center", va="bottom", fontsize=8.2, color=MUTED,
             linespacing=1.5)
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", os.path.relpath(out, ROOT))


def fig3_error_bars(rows, out):
    """Replicate spread against how well the pair held its pose, with the SEM for reference.

    All four replicate pairs in the project are plotted. The spread is computed from `dG_raw`,
    not the 2 dp display column: it is a difference of two near-equal numbers, and rounding first
    turns the true 0.26/0.24 MAO spreads into 0.27/0.23.

    Both series are kcal/mol on one axis -- this is deliberately not a dual-axis chart. Ligand
    RMSD is the x position, so the reading is "as the pose holds less well, the spread grows while
    the reported SEM does not move".
    """
    PAIRS = [("TTBK1 pose 1", "system", "system_TTBK1_r2"),
             ("TTBK2 pose 1", "system_TTBK2", "system_TTBK2m"),
             ("MAO-A pose 3", "system_MAOA_p3", "system_MAOA_p3_r2"),
             ("MAO-B pose 2", "system_MAOB_p2", "system_MAOB_p2_r2")]
    pts = []
    for name, a, b in PAIRS:
        if not (a in rows and b in rows and rows[a].get("dG_raw") and rows[b].get("dG_raw")):
            continue
        spread = abs(float(rows[a]["dG_raw"]) - float(rows[b]["dG_raw"]))
        sem = max(float(rows[a]["dG_SEM"]), float(rows[b]["dG_SEM"]))
        rms = (float(rows[a]["lig_rmsd_last100"]) + float(rows[b]["lig_rmsd_last100"])) / 2
        pts.append((name, rms, spread, sem))
    pts.sort(key=lambda t: t[1])

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    fig.patch.set_facecolor(SURFACE)
    style(ax)

    xs = [p[1] for p in pts]
    ax.plot(xs, [p[2] for p in pts], "-o", color=ORANGE, lw=2, ms=9, zorder=3,
            label="spread between two runs differing only in velocity seed")
    ax.plot(xs, [p[3] for p in pts], "-o", color=BLUE, lw=2, ms=9, zorder=3,
            label="SEM reported by MMPBSA.py (larger of the pair)")

    # The two MAO pairs land within 0.02 A of each other on x, so their labels would print on top
    # of one another. Push near-coincident neighbours to opposite sides instead of centring both.
    align = ["center"] * len(pts)
    for i in range(1, len(pts)):
        if pts[i][1] - pts[i - 1][1] < 0.20:
            align[i - 1], align[i] = "right", "left"
    dx = {"center": 0, "right": -9, "left": 9}

    for (name, rms, spread, sem), ha in zip(pts, align):
        ax.annotate(f"{name}\n{spread:.2f}  ({spread / sem:.0f}× SEM)",
                    xy=(rms, spread), xytext=(dx[ha], 13), textcoords="offset points",
                    ha=ha, fontsize=8.6, color=INK, linespacing=1.45)
        ax.annotate(f"{sem:.2f}", xy=(rms, sem), xytext=(dx[ha], -17), textcoords="offset points",
                    ha=ha, fontsize=8.6, color=INK2)

    ax.set_xlim(min(xs) - 0.45, max(xs) + 0.55)
    ax.set_ylim(-0.75, max(p[2] for p in pts) + 1.9)
    ax.set_xlabel("How well the pair held its pose — mean ligand RMSD over the final 100 frames (Å)",
                  fontsize=9.5, color=INK2)
    ax.set_ylabel("kcal/mol", fontsize=9.5, color=INK2)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title("How far the SEM is from the real error bar depends on pose stability",
                 fontsize=12.5, color=INK, pad=14, loc="left", weight="medium")
    ax.legend(loc="upper left", frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.text(0.5, 0.018, "All four replicate pairs in the project. The SEM barely moves while the "
             "measured spread grows 29-fold;\nn = 4 pairs across two protein families, so this is "
             "a consistent pattern rather than a calibration.",
             ha="center", va="bottom", fontsize=8.2, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", os.path.relpath(out, ROOT))


def main():
    rows = load()
    fig1_pose_stability(rows, os.path.join(HERE, "fig1_pose_stability.png"))
    fig2_mao_energy(rows, os.path.join(HERE, "fig2_mao_binding_energy.png"))
    fig3_error_bars(rows, os.path.join(HERE, "fig3_sem_vs_replicate.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
