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

# Data from canonical split (seed 42) analysis.
# The per-bin widths are loaded from the versioned CSV
# (tables/fig_s13_binwise_width_data.csv) rather than hardcoded inline, so the
# figure is reproducible from a tracked, diffable artifact. See PROVENANCE.md.
# Bin order: <=1d, 1d-1y, 1y-10y, 10y-100y, 100y-1ky, 1ky-10ky, >10ky
from pathlib import Path

HERE = Path(__file__).resolve().parent
_DATA_CSV = HERE.parent / "tables" / "fig_s13_binwise_width_data.csv"
_df = pd.read_csv(_DATA_CSV)

bins = list(dict.fromkeys(_df["bin"]))


def _series(column):
    """{target: [width_log10 per bin]}, with blank cells mapped to NaN."""
    out = {}
    for target, grp in _df.groupby("target", sort=False):
        vals = []
        for v in grp[column]:
            s = str(v).strip()
            vals.append(float(s) if s not in ("", "nan") else np.nan)
        out[target] = vals
    return out


# Median width (log10 scale) for each target and aggregation mode
data_none = _series("median_width_log10_no_aggregation")
data_max = _series("median_width_log10_max_aggregation")

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
output_path = HERE.parent / "figures" / "Fig_S13_binwise_width_comparison.pdf"
fig.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Figure saved to: {output_path}")

# Also save PNG for quick preview
png_path = output_path.with_suffix(".png")
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
