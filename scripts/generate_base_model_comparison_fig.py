#!/usr/bin/env python3
"""Generate base-model comparison figure for SC-PIML+ (HGBR-anchored).

Data source: main.tex Table base-model (SC-PIML+ base-learner comparison on
D224+act9, alpha=0.05). HGBR is the default base learner and is highlighted
across all three panels.

SCI-style choices:
  - Okabe-Ito colorblind-safe categorical palette for bars.
  - Cividis (perceptually uniform, colorblind-safe, B&W-printable) for heatmap.
  - Unified model ordering across all three panels (sorted by pass count,
    ties broken by mean WIS).
  - HGBR highlighted by bold black edge in (a)/(b) and by dashed boundary in (c).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

# ── Data (from main.tex Table base-model, D224+act9, alpha=0.05) ────────────
MODELS_RAW = ["HGBR", "XGBoost", "RF", "LightGBM", "MLP", "GP", "LSTM"]
PASS_RAW = [15, 14, 13, 11, 13, 10, 8]
MEAN_WIS_RAW = [1.08, 1.09, 1.10, 1.10, 1.63, 2.38, 2.88]

TARGETS = [
    r"$^3$H", r"$^6$Li", r"$^7$Li", "Decay heat", "Activity",
    "Contact dose", r"$^{210}$Po", r"$^{210}$Bi", r"$^{203}$Hg",
    r"$^{204}$Tl", r"$^{207}$Bi", r"$^{203}$Pb", r"$^{54}$Mn",
    r"$^{60}$Co", r"$^{55}$Fe",
]
PICP_RAW = np.array([
    [97.7, 96.7, 96.6, 97.4, 97.8, 95.6, 95.3],
    [99.7, 98.0, 95.8, 97.1, 97.1, 98.6, 97.1],
    [100.0, 100.0, 90.0, 96.1, 98.5, 99.8, 96.8],
    [98.8, 98.9, 96.7, 97.1, 98.2, 97.9, 95.9],
    [99.1, 97.9, 96.3, 98.8, 95.0, 96.4, 80.3],
    [98.8, 97.6, 96.4, 96.7, 96.3, 96.7, 91.9],
    [96.4, 96.2, 95.8, 94.9, 97.9, 97.4, 81.9],
    [97.8, 95.0, 97.1, 94.9, 99.5, 97.6, 99.3],
    [97.7, 97.6, 95.6, 97.2, 99.0, 93.9, 96.4],
    [95.9, 96.0, 96.0, 94.4, 94.0, 94.2, 94.8],
    [97.6, 95.8, 98.0, 96.2, 96.0, 96.6, 87.8],
    [97.6, 97.7, 95.0, 96.6, 96.9, 94.2, 96.3],
    [95.9, 97.6, 97.2, 93.8, 99.4, 90.0, 90.0],
    [98.8, 96.8, 97.8, 96.1, 97.5, 97.2, 95.2],
    [96.4, 97.3, 94.5, 98.5, 97.2, 94.8, 94.3],
])

# ── Unified ordering: sort by pass desc, then mean WIS asc ──────────────────
order = sorted(range(len(MODELS_RAW)),
               key=lambda i: (-PASS_RAW[i], MEAN_WIS_RAW[i]))
MODELS = [MODELS_RAW[i] for i in order]
PASS = [PASS_RAW[i] for i in order]
MEAN_WIS = [MEAN_WIS_RAW[i] for i in order]
PICP = PICP_RAW[:, order]
ANCHOR = "HGBR"
anchor_idx = MODELS.index(ANCHOR)

# ── Colorblind-safe categorical palette (Okabe-Ito) ─────────────────────────
OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
             "#56B4E9", "#D55E00", "#F0E442"]
COLORS = {m: OKABE_ITO[i] for i, m in enumerate(MODELS)}

# ── Figure ──────────────────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5),
                         gridspec_kw={"width_ratios": [1, 1, 1.3]})
n = len(MODELS)
x = np.arange(n)


def _highlight(bar):
    bar.set_edgecolor("black")
    bar.set_linewidth(2.2)


# --- Panel (a): Pass count ---
ax = axes[0]
bars = ax.bar(x, PASS, color=[COLORS[m] for m in MODELS],
              edgecolor="white", linewidth=0.6)
_highlight(bars[anchor_idx])
for i, (v, m) in enumerate(zip(PASS, MODELS)):
    ax.text(i, v + 0.35, str(v), ha="center", va="bottom", fontsize=9,
            fontweight="bold" if m == ANCHOR else "normal")
ax.axhline(y=15, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
ax.set_xticks(x)
ax.set_xticklabels(MODELS, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Pass count (≥95% PICP)", fontsize=10)
ax.set_title("(a) Coverage pass count", fontsize=11, fontweight="bold")
ax.set_ylim(0, 18)
ax.yaxis.set_major_locator(mticker.MultipleLocator(3))
ax.grid(axis="y", alpha=0.3, linewidth=0.5)
ax.get_xticklabels()[anchor_idx].set_fontweight("bold")

# --- Panel (b): Mean WIS ---
ax = axes[1]
bars2 = ax.bar(x, MEAN_WIS, color=[COLORS[m] for m in MODELS],
               edgecolor="white", linewidth=0.6)
_highlight(bars2[anchor_idx])
for i, (v, m) in enumerate(zip(MEAN_WIS, MODELS)):
    ax.text(i, v + 0.06, f"{v:.2f}", ha="center", va="bottom", fontsize=9,
            fontweight="bold" if m == ANCHOR else "normal")
ax.set_xticks(x)
ax.set_xticklabels(MODELS, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Mean WIS (↓ better)", fontsize=10)
ax.set_title("(b) Interval efficiency", fontsize=11, fontweight="bold")
ax.set_ylim(0, 3.3)
ax.grid(axis="y", alpha=0.3, linewidth=0.5)
ax.get_xticklabels()[anchor_idx].set_fontweight("bold")

# --- Panel (c): Per-target PICP heatmap (cividis) ---
ax = axes[2]
im = ax.imshow(PICP, aspect="auto", cmap="cividis",
               vmin=80, vmax=100, origin="upper")
for i in range(len(TARGETS)):
    for j in range(len(MODELS)):
        val = PICP[i, j]
        text_color = "white" if val < 92 else "black"
        ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                fontsize=7.5, color=text_color,
                fontweight="bold" if val < 95 else "normal")
ax.set_xticks(x)
ax.set_xticklabels(MODELS, rotation=30, ha="right", fontsize=9)
ax.set_yticks(np.arange(len(TARGETS)))
ax.set_yticklabels(TARGETS, fontsize=8.5)
ax.set_title("(c) Per-target PICP heatmap", fontsize=11, fontweight="bold")
ax.get_xticklabels()[anchor_idx].set_fontweight("bold")
ax.axvline(x=anchor_idx - 0.5, color="black", linewidth=1.6, linestyle="--")
ax.axvline(x=anchor_idx + 0.5, color="black", linewidth=1.6, linestyle="--")
ax.text(anchor_idx, -1.0, "default", ha="center", va="bottom",
        fontsize=8.5, fontweight="bold", color="black",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#F0E442",
                  edgecolor="black", linewidth=0.8))
cbar = plt.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
cbar.set_label("PICP (%)", fontsize=9)
cbar.ax.axhline(y=95, color="red", linewidth=1.4, linestyle="--")

plt.tight_layout(pad=1.5, w_pad=2.0)

# ── Save ────────────────────────────────────────────────────────────────────
png_path = os.path.join(OUT_DIR, "Fig_P3_base_model_comparison.png")
pdf_path = os.path.join(OUT_DIR, "Fig_P3_base_model_comparison.pdf")
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"✅ Saved: {png_path}")
print(f"✅ Saved: {pdf_path}")
plt.close(fig)
