#!/usr/bin/env python3
"""Regenerate Fig 3: Actual temporal evolution curves for 15 targets on D224.

Improved version:
  - 5×3 grid (all 15 slots filled, no blank cell)
  - Li6/Li7 shown as % depletion (more intuitive than raw log-ratio)
  - Irradiation/cooling phase boundary marked with vertical dashed line
  - Legend placed outside axes to avoid overlap
  - Consistent formatting across panels

Output:
  figures/Fig_P3_actual_evolution_curves.pdf
  figures/Fig_P3_actual_evolution_curves.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Computer Modern Roman"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.formatter.use_mathtext": True,
    "legend.fontsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.4,
})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
DATASET_DIR = REPO_ROOT / "data" / "dataset"
LABELS_CSV = DATASET_DIR / "labels_long_bzsum_elapsed_total_act9_with_yref_logratio.csv"
SPLITS_JSON = DATASET_DIR / "splits_d224_seed42.json"

# PAPER_ROOT already set as REPO_ROOT above
FIGURES_DIR = REPO_ROOT / "figures"

# ---------------------------------------------------------------------------
# Irradiation duration (years) — phase boundary
# ---------------------------------------------------------------------------
T_IRRAD = 10.0

# ---------------------------------------------------------------------------
# Target definitions — ordered by category, 5 cols × 3 rows
# Row 1: Breeding (3) + Aggregate (2)
# Row 2: Aggregate (1) + PbLi-origin (4)   [dose completes aggregate]
# Wait — let me do a cleaner grouping:
# Row 1: H3, Li6, Li7, Decay_heat, Activity         (breeding + aggregate)
# Row 2: Contact_dose, Po210, Bi210, Hg203, Tl204  (aggregate + PbLi)
# Row 3: Bi207, Pb203, Mn54, Co60, Fe55              (PbLi + steel)
# ---------------------------------------------------------------------------
TARGETS = [
    # (column_name, display_label, unit, category, color)
    # Row 1
    ("atoms_H3",           r"$^{3}$H",         "atoms", "Breeding",        "#1f77b4"),
    ("atoms_Li6_logratio", r"$^{6}$Li",         None,    "Breeding",        "#1f77b4"),
    ("atoms_Li7_logratio", r"$^{7}$Li",         None,    "Breeding",        "#1f77b4"),
    ("decay_heat_W",       "Decay heat",        "W",     "Aggregate",       "#2ca02c"),
    ("activity_Bq",        "Activity",          "Bq",    "Aggregate",       "#2ca02c"),
    # Row 2
    ("contact_dose_Svph",  "Contact dose",      "Sv/h",  "Aggregate",       "#2ca02c"),
    ("atoms_Po210",        r"$^{210}$Po",       "atoms", "PbLi activation", "#ff7f0e"),
    ("atoms_Bi210",        r"$^{210}$Bi",       "atoms", "PbLi activation", "#ff7f0e"),
    ("atoms_Hg203",        r"$^{203}$Hg",       "atoms", "PbLi activation", "#ff7f0e"),
    ("atoms_Tl204",        r"$^{204}$Tl",       "atoms", "PbLi activation", "#ff7f0e"),
    # Row 3
    ("atoms_Bi207",        r"$^{207}$Bi",       "atoms", "PbLi activation", "#ff7f0e"),
    ("atoms_Pb203",        r"$^{203}$Pb",       "atoms", "PbLi activation", "#ff7f0e"),
    ("atoms_Mn54",         r"$^{54}$Mn",        "atoms", "Steel activation","#d62728"),
    ("atoms_Co60",         r"$^{60}$Co",        "atoms", "Steel activation","#d62728"),
    ("atoms_Fe55",         r"$^{55}$Fe",        "atoms", "Steel activation","#d62728"),
]

CATEGORY_COLORS = {
    "Breeding":        "#1f77b4",
    "Aggregate":       "#2ca02c",
    "PbLi activation": "#ff7f0e",
    "Steel activation":"#d62728",
}


def main():
    # Load splits
    with open(SPLITS_JSON) as f:
        splits = json.load(f)
    test_ids = splits["test"]
    sample_id = test_ids[0]
    print(f"Using test sample: {sample_id}")

    # Load data
    cols_needed = ["sample_id", "time_years"] + [t[0] for t in TARGETS]
    df = pd.read_csv(LABELS_CSV, usecols=cols_needed)
    sub = df[df["sample_id"] == sample_id].sort_values("time_years").copy()
    print(f"  {len(sub)} time points, t=[{sub['time_years'].min():.3f}, {sub['time_years'].max():.1f}] y")

    t = sub["time_years"].values

    # ---------------------------------------------------------------------------
    # Figure: 3 rows × 5 cols = 15 panels, all filled
    # ---------------------------------------------------------------------------
    nrows, ncols = 3, 5
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 7.5))

    for i, (col, label, unit, category, color) in enumerate(TARGETS):
        row, c = divmod(i, ncols)
        ax = axes[row, c]

        y = sub[col].values

        # Transform Li6/Li7 log-ratio → percentage depletion for display
        if "logratio" in col:
            # log-ratio = log10(y/y_ref), convert to % change: (10^lr - 1) * 100
            y_display = (np.power(10, y) - 1) * 100
            ax.plot(t, y_display, color=color, linewidth=1.4)
            ax.set_yscale("linear")
            title_str = f"{label} depletion (%)"
        else:
            ax.plot(t, y, color=color, linewidth=1.4)
            if unit:
                title_str = f"{label} ({unit})"
            else:
                title_str = label
            # Log scale for positive quantities
            y_pos = y[y > 0]
            if len(y_pos) > 0:
                ax.set_yscale("log")
            else:
                ax.set_yscale("linear")

        ax.set_title(title_str, fontsize=9, pad=3)
        ax.set_xscale("log")
        ax.grid(True, alpha=0.2, linewidth=0.5)

        # Irradiation/cooling phase boundary
        ax.axvline(T_IRRAD, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

        # Light shading for irradiation phase
        ax.axvspan(t.min(), T_IRRAD, color="#e8e8e8", alpha=0.3, zorder=0)

        # Panel label
        letter = chr(ord('a') + i)
        ax.text(0.03, 0.93, f"({letter})", transform=ax.transAxes,
                fontweight="bold", fontsize=9, va="top")

    # X-axis labels only on bottom row
    for row in range(nrows):
        for c in range(ncols):
            if row < nrows - 1:
                axes[row, c].tick_params(labelbottom=False)
            else:
                axes[row, c].set_xlabel("Time (years)", fontsize=9)

    # Shared x-limits
    for ax in axes.flat:
        ax.set_xlim(t.min() * 0.9, t.max() * 1.1)

    # Category legend — horizontal bar at the bottom, outside plot area
    legend_handles = []
    for name, color in CATEGORY_COLORS.items():
        legend_handles.append(Line2D([0], [0], color=color, linewidth=2.5, label=name))
    # Add phase boundary indicator
    legend_handles.append(Line2D([0], [0], color="gray", linestyle="--",
                                 linewidth=0.8, alpha=0.6, label="Shutdown"))

    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01), ncol=5, frameon=True,
               fontsize=9, edgecolor="lightgray", fancybox=False,
               columnspacing=1.5, handletextpad=0.5)

    fig.subplots_adjust(hspace=0.35, wspace=0.35, bottom=0.1, top=0.96,
                        left=0.05, right=0.98)

    # Save
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    out_pdf = FIGURES_DIR / "Fig_P3_actual_evolution_curves.pdf"
    out_png = FIGURES_DIR / "Fig_P3_actual_evolution_curves.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {out_pdf}")
    print(f"WROTE {out_png}")


if __name__ == "__main__":
    main()
