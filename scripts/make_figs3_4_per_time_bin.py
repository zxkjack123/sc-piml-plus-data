#!/usr/bin/env python3
"""Regenerate Fig 4 & 5: Per-time and per-bin coverage/width figures for D224.

Improvements over original:
  - Per-time figure: step-line instead of dense scatter; 95% nominal line;
    clean axis labels; legend outside plot; irradiation/cooling shading
  - Per-bin figure: grouped bar chart; legend outside plot area
  - Width panels use "excess width factor (w − 1)" for Li targets to avoid
    confusing "1+2.4e-4" labels

Output:
  figures/Fig_P3_timeaware_per_time_d224.pdf/.png
  figures/Fig_P3_timeaware_per_bin_d224.pdf/.png
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Computer Modern Roman"],
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.3,
    "lines.markersize": 3.5,
})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

# ---------------------------------------------------------------------------
# Paths (same as original)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
REPORT_ROOT = REPO_ROOT / "data" / "d224_runs"
# PAPER_ROOT = REPO_ROOT (set above)
FIGURES_DIR = REPO_ROOT / "figures"

EPS = 1.0
T_SHUTDOWN = 10.0  # years

D224_GROUPING_RUNS = {
    "atoms_H3": {
        "Global": REPORT_ROOT / "d224_cvplus_global_H3",
        "Time-binned": REPORT_ROOT / "d224_cvplus_smoke",
        "Phase×time": REPORT_ROOT / "d224_cvplus_phase_timebin_H3",
    },
    "atoms_Li6": {
        "Global": REPORT_ROOT / "d224_cvplus_global_Li6",
        "Time-binned": REPORT_ROOT / "d224_cvplus_Li6",
        "Phase×time": REPORT_ROOT / "d224_cvplus_phase_timebin_Li6",
    },
    "atoms_Li7": {
        "Global": REPORT_ROOT / "d224_cvplus_global_Li7",
        "Time-binned": REPORT_ROOT / "d224_cvplus_Li7",
        "Phase×time": REPORT_ROOT / "d224_cvplus_phase_timebin_Li7",
    },
}

BIN_EDGES_YEARS = [0.0, 1.0 / 365.0, 1.0, 10.0, 100.0, 1000.0, 10000.0, float("inf")]
BIN_LABELS = ["<=1d", "1d-1y", "1y-10y", "10y-100y", "100y-1ky", "1ky-10ky", ">10ky"]

TARGETS = ["atoms_H3", "atoms_Li6", "atoms_Li7"]
TARGET_PRETTY = {
    "atoms_H3": r"$^{3}$H",
    "atoms_Li6": r"$^{6}$Li",
    "atoms_Li7": r"$^{7}$Li",
}

METHOD_COLORS = {"Global": "#1f77b4", "Time-binned": "#ff7f0e", "Phase×time": "#2ca02c"}
METHOD_MARKERS = {"Global": "o", "Time-binned": "s", "Phase×time": "D"}


# ---------------------------------------------------------------------------
# Data loading (reused from original)
# ---------------------------------------------------------------------------
def _read_eps(run_dir: Path) -> float:
    cfg_path = run_dir / "timeaware_uq_config.json"
    if cfg_path.exists():
        with open(cfg_path) as f:
            return float(json.load(f).get("eps", 1.0))
    return 1.0


def _read_shutdown_years(run_dir: Path) -> float:
    report_path = run_dir / "timeaware_uq_report.md"
    if report_path.exists():
        txt = report_path.read_text(encoding="utf-8")
        m = re.search(r"^\- shutdown_years \(detected\):\s*([0-9eE+\-\.]+)\s*$", txt, flags=re.MULTILINE)
        if m:
            return float(m.group(1))
    cfg_path = run_dir / "timeaware_uq_config.json"
    if cfg_path.exists():
        with open(cfg_path) as f:
            cfg = json.load(f)
        val = cfg.get("midterm_shutdown_years") or cfg.get("midterm_max_irr_years")
        if val is not None:
            return float(val)
    return 10.0


def assign_time_bin(t_yr: float) -> str:
    for (lo, hi), label in zip(zip(BIN_EDGES_YEARS[:-1], BIN_EDGES_YEARS[1:]), BIN_LABELS):
        if t_yr == 0.0 and lo == 0.0:
            return label
        if lo < t_yr <= hi:
            return label
    return BIN_LABELS[-1]


def per_time_width_from_predictions(run_dir: Path, target: str) -> pd.DataFrame:
    eps = _read_eps(run_dir)
    p = pd.read_csv(run_dir / f"predictions_test_{target}.csv")
    yhat = np.maximum(p["y_pred"].astype(float).to_numpy(), 0.0) + eps
    ylo = np.maximum(p["y_pred_lo_cal"].astype(float).to_numpy(), 0.0) + eps
    yhi = np.maximum(p["y_pred_hi_cal"].astype(float).to_numpy(), 0.0) + eps
    width_factor = np.maximum(yhi / yhat, yhat / ylo)
    tmp = pd.DataFrame({
        "time_years": p["time_years"].astype(float).to_numpy(),
        "uq_half_width_factor": width_factor,
    })
    return (
        tmp.groupby("time_years", as_index=False)
        .agg(uq_half_width_factor=("uq_half_width_factor", "median"))
        .sort_values("time_years")
    )


def load_per_time_df(run_dir: Path, target: str, method: str) -> pd.DataFrame:
    df = pd.read_csv(run_dir / f"metrics_per_time_{target}.csv").copy()
    df["time_bin"] = df["time_years"].map(assign_time_bin)
    w = per_time_width_from_predictions(run_dir, target)
    df = df.merge(w, on="time_years", how="left")
    df = df.drop_duplicates(subset=["time_years"], keep="first")
    df["method"] = method
    df["target"] = target
    return df


def load_all_data():
    """Load all per-time data for the three targets and three methods."""
    frames = []
    methods = list(METHOD_COLORS.keys())
    for target in TARGETS:
        for method in methods:
            rd = D224_GROUPING_RUNS[target][method]
            if not rd.exists():
                print(f"  [WARN] missing: {rd}")
                continue
            frames.append(load_per_time_df(rd, target, method))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ---------------------------------------------------------------------------
# Figure 4: Per-time coverage + width
# ---------------------------------------------------------------------------
def make_per_time_figure(all_df: pd.DataFrame) -> None:
    print("\n=== Fig_P3_timeaware_per_time_d224 ===")

    methods = list(METHOD_COLORS.keys())
    nrows = len(TARGETS)

    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12, 7.5),
                             gridspec_kw={"width_ratios": [1, 1]})

    for r, target in enumerate(TARGETS):
        ax_cov = axes[r, 0]
        ax_wid = axes[r, 1]

        for method in methods:
            sub = all_df[(all_df["target"] == target) & (all_df["method"] == method)].sort_values("time_years")
            if sub.empty:
                continue
            t = sub["time_years"].astype(float).values
            cov = sub["coverage_90_cal"].astype(float).values
            wf = sub["uq_half_width_factor"].astype(float).values
            color = METHOD_COLORS[method]

            # Step-line for coverage (better than scatter for time-series)
            ax_cov.step(t, cov, where="mid", color=color, linewidth=1.2, alpha=0.85)
            # Sparse markers at sampled positions (every ~8th point)
            idx = np.linspace(0, len(t) - 1, min(15, len(t)), dtype=int)
            ax_cov.plot(t[idx], cov[idx], marker=METHOD_MARKERS[method],
                        color=color, linestyle="none", markersize=4, alpha=0.7)

            # Width: line plot
            ax_wid.plot(t, wf, color=color, linewidth=1.2, alpha=0.85)
            ax_wid.plot(t[idx], wf[idx], marker=METHOD_MARKERS[method],
                        color=color, linestyle="none", markersize=4, alpha=0.7)

        # Coverage panel configuration
        ax_cov.axhline(0.95, color="gray", ls="--", lw=0.8, alpha=0.6, zorder=0)
        ax_cov.set_ylim(0.45, 1.05)
        ax_cov.set_ylabel(f"{TARGET_PRETTY[target]} PICP")
        ax_cov.set_xscale("log")
        ax_cov.grid(True, alpha=0.2, linewidth=0.5)

        # Irradiation/cooling shading
        xlims = ax_cov.get_xlim()
        ax_cov.axvspan(xlims[0], T_SHUTDOWN, color="#e8e8e8", alpha=0.3, zorder=0)
        ax_cov.axvline(T_SHUTDOWN, color="gray", ls="--", lw=0.7, alpha=0.5)

        # Width panel configuration
        ax_wid.set_xscale("log")
        ax_wid.set_ylabel(f"{TARGET_PRETTY[target]} width factor")

        # Determine if we need "excess width" formatting
        lines_data = all_df[(all_df["target"] == target)]
        wf_all = lines_data["uq_half_width_factor"].dropna()
        w_range = wf_all.max() - wf_all.min()

        if w_range < 0.01:
            # For Li targets: show as "excess width factor (w − 1)"
            ax_wid.set_ylabel(f"{TARGET_PRETTY[target]}\nexcess width ($w - 1$)")
            # Replot in excess space
            ax_wid.cla()
            ax_wid.set_xscale("log")
            for method in methods:
                sub = all_df[(all_df["target"] == target) & (all_df["method"] == method)].sort_values("time_years")
                if sub.empty:
                    continue
                t = sub["time_years"].astype(float).values
                wf = sub["uq_half_width_factor"].astype(float).values - 1.0
                color = METHOD_COLORS[method]
                ax_wid.plot(t, wf, color=color, linewidth=1.2, alpha=0.85)
                idx = np.linspace(0, len(t) - 1, min(15, len(t)), dtype=int)
                ax_wid.plot(t[idx], wf[idx], marker=METHOD_MARKERS[method],
                            color=color, linestyle="none", markersize=4, alpha=0.7)
            ax_wid.set_ylabel(f"{TARGET_PRETTY[target]}\nexcess width ($w - 1$)")
            # Scientific notation for small values
            ax_wid.ticklabel_format(axis='y', style='scientific', scilimits=(-3, -3))
        else:
            # For H3: normal scale
            pass

        ax_wid.grid(True, alpha=0.2, linewidth=0.5)

        # Irradiation shading on width panel too
        xlims_w = ax_wid.get_xlim()
        ax_wid.axvspan(xlims_w[0], T_SHUTDOWN, color="#e8e8e8", alpha=0.3, zorder=0)
        ax_wid.axvline(T_SHUTDOWN, color="gray", ls="--", lw=0.7, alpha=0.5)

        # Panel labels
        ax_cov.text(0.03, 0.93, f"({chr(ord('a') + r * 2)})", transform=ax_cov.transAxes,
                    fontweight="bold", fontsize=10, va="top")
        ax_wid.text(0.03, 0.93, f"({chr(ord('a') + r * 2 + 1)})", transform=ax_wid.transAxes,
                    fontweight="bold", fontsize=10, va="top")

    # Column titles
    axes[0, 0].set_title("Coverage (PICP per time point)", fontsize=10)
    axes[0, 1].set_title("Interval width factor per time point", fontsize=10)

    # X-labels only on bottom row
    for c in range(2):
        for r in range(nrows - 1):
            axes[r, c].tick_params(labelbottom=False)
        axes[nrows - 1, c].set_xlabel("Time (years)")

    # Legend at bottom, outside plot area
    legend_handles = []
    for method in methods:
        legend_handles.append(Line2D([0], [0], color=METHOD_COLORS[method],
                                     marker=METHOD_MARKERS[method], markersize=5,
                                     linewidth=1.5, label=method))
    legend_handles.append(Line2D([0], [0], color="gray", ls="--", lw=0.8, label="95% nominal"))
    legend_handles.append(Line2D([0], [0], color="gray", ls="--", lw=0.7, alpha=0.5, label="Shutdown"))

    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=5, frameon=True,
               fontsize=9, edgecolor="lightgray", fancybox=False,
               columnspacing=1.5, handletextpad=0.5)

    fig.subplots_adjust(hspace=0.2, wspace=0.3, bottom=0.1, top=0.94,
                        left=0.08, right=0.97)

    # Save
    for suffix in ["d224"]:
        out_pdf = FIGURES_DIR / f"Fig_P3_timeaware_per_time_{suffix}.pdf"
        out_png = FIGURES_DIR / f"Fig_P3_timeaware_per_time_{suffix}.png"
        fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        print(f"  WROTE {out_pdf}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: Per-bin aggregated coverage + width
# ---------------------------------------------------------------------------
def make_per_bin_figure(all_df: pd.DataFrame) -> None:
    print("\n=== Fig_P3_timeaware_per_bin_d224 ===")

    methods = list(METHOD_COLORS.keys())

    # Aggregate per bin, per method, PER TARGET first, then average across targets
    all_df = all_df.copy()
    all_df["time_bin"] = all_df["time_years"].map(assign_time_bin)

    # Per-target, per-bin, per-method aggregation
    target_bin_agg = (
        all_df.groupby(["target", "method", "time_bin"], as_index=False)
        .agg(
            coverage=("coverage_90_cal", "mean"),
            width_median=("uq_half_width_factor", "median"),
        )
    )
    # Then average across 3 targets
    bin_agg = (
        target_bin_agg.groupby(["method", "time_bin"], as_index=False)
        .agg(
            coverage=("coverage", "mean"),
            width_median=("width_median", "mean"),
        )
    )

    bins = [b for b in BIN_LABELS if b in set(bin_agg["time_bin"])]
    x = np.arange(len(bins))
    bar_width = 0.25

    fig, (ax_cov, ax_wid) = plt.subplots(1, 2, figsize=(12, 3.8))

    for i, method in enumerate(methods):
        sub = bin_agg[bin_agg["method"] == method].set_index("time_bin")
        y_cov = [float(sub.loc[b, "coverage"]) if b in sub.index else np.nan for b in bins]
        y_w = [float(sub.loc[b, "width_median"]) if b in sub.index else np.nan for b in bins]
        offset = (i - 1) * bar_width

        ax_cov.bar(x + offset, y_cov, bar_width, color=METHOD_COLORS[method],
                   alpha=0.85, label=method, edgecolor="white", linewidth=0.5)
        ax_wid.bar(x + offset, y_w, bar_width, color=METHOD_COLORS[method],
                   alpha=0.85, edgecolor="white", linewidth=0.5)

    # Coverage panel
    ax_cov.axhline(0.95, color="gray", ls="--", lw=0.8, alpha=0.6)
    ax_cov.set_ylim(0.55, 1.05)
    ax_cov.set_ylabel("PICP (bin-averaged)")
    ax_cov.set_xticks(x)
    ax_cov.set_xticklabels(bins, rotation=30, ha="right", fontsize=8)
    ax_cov.set_title("Coverage by time bin", fontsize=10)
    ax_cov.grid(axis="y", alpha=0.2, linewidth=0.5)
    ax_cov.text(0.02, 0.95, "(a)", transform=ax_cov.transAxes,
                fontweight="bold", fontsize=10, va="top")

    # Width panel — use log scale to reveal structure near 1.0
    ax_wid.set_yscale("log")
    from matplotlib.ticker import ScalarFormatter
    sf = ScalarFormatter()
    sf.set_scientific(False)
    ax_wid.yaxis.set_major_formatter(sf)
    ax_wid.yaxis.set_minor_formatter(sf)
    ax_wid.set_ylabel("Width factor (bin-median)")
    ax_wid.set_xticks(x)
    ax_wid.set_xticklabels(bins, rotation=30, ha="right", fontsize=8)
    ax_wid.set_title("Interval width by time bin", fontsize=10)
    ax_wid.grid(axis="y", alpha=0.2, linewidth=0.5)
    ax_wid.text(0.02, 0.95, "(b)", transform=ax_wid.transAxes,
                fontweight="bold", fontsize=10, va="top")

    # Legend at bottom center, outside — only from coverage axis (avoid duplicates)
    handles, labels = ax_cov.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.05),
               ncol=3, frameon=True, fontsize=9, edgecolor="lightgray",
               fancybox=False, columnspacing=1.5)

    fig.subplots_adjust(bottom=0.18, wspace=0.3)

    for suffix in ["d224"]:
        out_pdf = FIGURES_DIR / f"Fig_P3_timeaware_per_bin_{suffix}.pdf"
        out_png = FIGURES_DIR / f"Fig_P3_timeaware_per_bin_{suffix}.png"
        fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        print(f"  WROTE {out_pdf}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    all_df = load_all_data()
    if all_df.empty:
        print("No data loaded!")
        sys.exit(1)
    print(f"Loaded {len(all_df)} rows across {all_df['target'].nunique()} targets, "
          f"{all_df['method'].nunique()} methods")

    make_per_time_figure(all_df)
    make_per_bin_figure(all_df)


if __name__ == "__main__":
    main()
