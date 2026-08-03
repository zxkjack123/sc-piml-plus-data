#!/usr/bin/env python3
"""
Generate Figure S5: Per-bin width comparison between no-aggregation and max-aggregation modes.
Shows median interval width (log10 scale) for each temporal bin across all targets.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

# Data from canonical split (seed 42) analysis
# Bin order: <=1d, 1d-1y, 1y-10y, 10y-100y, 100y-1ky, 1ky-10ky, >10ky
bins = ["<=1d", "1d-1y", "1y-10y", "10y-100y", "100y-1ky", "1ky-10ky", ">10ky"]

# Median width (log10 scale) for each target and aggregation mode
# No aggregation (pooled)
data_none = {
    "atoms_H3": [np.nan, 0.1630, 0.0848, 0.0911, 1.0915, 5.6025, 5.6025],
    "atoms_Li6": [np.nan, 0.0000, 0.0001, 0.0002, 0.0002, 0.0002, 0.0002],
    "atoms_Li7": [np.nan, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    "decay_heat_W": [np.nan, 0.0364, 0.1549, 0.0020, 0.0018, 0.0018, 0.0018],
}

# Max aggregation (geometry-level)
data_max = {
    "atoms_H3": [np.nan, 0.3498, 0.1728, 0.1897, 0.5730, 4.5761, 4.5761],
    "atoms_Li6": [np.nan, 0.0000, 0.0003, 0.0002, 0.0002, 0.0002, 0.0002],
    "atoms_Li7": [np.nan, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    "decay_heat_W": [np.nan, 0.0285, 0.2700, 0.0040, 0.0016, 0.0016, 0.0016],
}

# Compute width ratios (max / none)
ratios = {}
for target in data_none.keys():
    ratios[target] = []
    for w_none, w_max in zip(data_none[target], data_max[target]):
        if np.isnan(w_none) or np.isnan(w_max) or w_none == 0:
            ratios[target].append(np.nan)
        else:
            ratios[target].append(w_max / w_none)

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

targets = list(data_none.keys())
colors_none = "#1f77b4"  # blue
colors_max = "#d62728"  # red

for idx, target in enumerate(targets):
    ax = axes[idx]

    # Plot bars
    x = np.arange(len(bins))
    width = 0.35

    bars_none = ax.bar(
        x - width / 2,
        data_none[target],
        width,
        label="No aggregation",
        color=colors_none,
        alpha=0.8,
    )
    bars_max = ax.bar(
        x + width / 2,
        data_max[target],
        width,
        label="Max aggregation",
        color=colors_max,
        alpha=0.8,
    )

    # Formatting
    ax.set_xlabel("Temporal bin", fontsize=11)
    ax.set_ylabel("Median width (log$_{10}$ scale)", fontsize=11)
    ax.set_title(target.replace("_", " ").title(), fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(bins, rotation=45, ha="right", fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    # Add ratio annotations for bins where both values are non-zero
    for i, (w_none, w_max, ratio) in enumerate(
        zip(data_none[target], data_max[target], ratios[target])
    ):
        if not np.isnan(ratio) and ratio > 0:
            max_height = max(w_none, w_max)
            ax.text(
                i,
                max_height * 1.1,
                f"{ratio:.1f}×",
                ha="center",
                va="bottom",
                fontsize=8,
                color="darkgreen",
            )

plt.tight_layout()

# Save figure
from pathlib import Path

output_path = (
    Path(__file__).resolve().parents[1]
    / "figures"
    / "Fig_S13_binwise_width_comparison.pdf"
)
fig.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Figure saved to: {output_path}")

# Also save PNG for quick preview
png_path = output_path.replace(".pdf", ".png")
fig.savefig(png_path, dpi=150, bbox_inches="tight")
print(f"PNG preview saved to: {png_path}")

# Print summary statistics
print("\n=== Width Ratio Summary (Max / No Aggregation) ===")
for target in targets:
    valid_ratios = [r for r in ratios[target] if not np.isnan(r)]
    if valid_ratios:
        print(f"{target}:")
        print(f"  Mean ratio: {np.mean(valid_ratios):.2f}×")
        print(f"  Range: {np.min(valid_ratios):.2f}× - {np.max(valid_ratios):.2f}×")
