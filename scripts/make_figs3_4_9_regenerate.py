#!/usr/bin/env python3
"""Regenerate all D224-backed figures for Paper 3.

Generates:
  1. Fig_P3_method_comparison.pdf      — bar chart of Global PICP@90 by method × target
  2. Fig_P3_timeaware_per_time_n80.pdf — per-time coverage + width for D224 (reuses n80 filename)
  3. Fig_P3_timeaware_per_bin_n80.pdf  — per-bin aggregated coverage + width for D224
  4. Fig_P3_trajectory_uq_bands_atoms_H3.pdf  — trajectory UQ bands (D224)
  5. Fig_P3_trajectory_uq_bands_atoms_Li6.pdf
  6. Fig_P3_trajectory_uq_bands_atoms_Li7.pdf
  + corresponding .png previews and summary CSVs

Created: 2026-04-02
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.constrained_layout.use": True,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
})
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator, ScalarFormatter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
REPORT_ROOT = REPO_ROOT / "data" / "d224_runs"
# PAPER_ROOT = REPO_ROOT (set above)
FIGURES_DIR = REPO_ROOT / "figures"
TABLES_DIR = REPO_ROOT / "tables"

EPS = 1.0

# ---------------------------------------------------------------------------
# D224 run directory map
# ---------------------------------------------------------------------------
# Per-time/per-bin grouping runs: (grouping_strategy -> run_dir) for each target
D224_GROUPING_RUNS: Dict[str, Dict[str, Path]] = {
    "atoms_H3": {
        "global": REPORT_ROOT / "d224_cvplus_global_H3",
        "time_bin": REPORT_ROOT / "d224_cvplus_smoke",
        "phase_time_bin": REPORT_ROOT / "d224_cvplus_phase_timebin_H3",
    },
    "atoms_Li6": {
        "global": REPORT_ROOT / "d224_cvplus_global_Li6",
        "time_bin": REPORT_ROOT / "d224_cvplus_Li6",
        "phase_time_bin": REPORT_ROOT / "d224_cvplus_phase_timebin_Li6",
    },
    "atoms_Li7": {
        "global": REPORT_ROOT / "d224_cvplus_global_Li7",
        "time_bin": REPORT_ROOT / "d224_cvplus_Li7",
        "phase_time_bin": REPORT_ROOT / "d224_cvplus_phase_timebin_Li7",
    },
}

# Method comparison: target -> {method_label -> (run_dir, target_col)}
# We read global_picp from metrics.json or from d224_baseline_comparison_results.csv
D224_METHOD_COMP_PROPOSED = {
    "atoms_H3": {
        "Proposed CV+\n(time-bin)": REPORT_ROOT / "d224_cvplus_smoke",
    },
    "atoms_Li6": {
        "Proposed CV+\n(time-bin)": REPORT_ROOT / "d224_cvplus_Li6",
    },
    "atoms_Li7": {
        "Proposed CV+\n(time-bin)": REPORT_ROOT / "d224_cvplus_Li7",
    },
    "decay_heat_W": {
        "Proposed CV+\n(time-bin)": REPORT_ROOT / "d224_cvplus_decay_heat",
    },
}

# Canonical time bins
BIN_EDGES_YEARS = [0.0, 1.0 / 365.0, 1.0, 10.0, 100.0, 1000.0, 10000.0, float("inf")]
BIN_LABELS = ["<=1d", "1d-1y", "1y-10y", "10y-100y", "100y-1ky", "1ky-10ky", ">10ky"]
# LaTeX-safe display labels (usetex=True renders > as ¿ in text mode)
BIN_DISPLAY = {b: b.replace("<=", r"$\leq$").replace(">", r"$>$") for b in BIN_LABELS}

TARGETS = ["atoms_H3", "atoms_Li6", "atoms_Li7"]

TARGET_PRETTY = {
    "atoms_H3": r"$^{3}$H",
    "atoms_Li6": r"$^{6}$Li",
    "atoms_Li7": r"$^{7}$Li",
    "decay_heat_W": "Decay heat",
}

TARGET_DISPLAY = {
    "atoms_H3": r"Atoms of $^{3}$H",
    "atoms_Li6": r"Atoms of $^{6}$Li",
    "atoms_Li7": r"Atoms of $^{7}$Li",
    "decay_heat_W": "Decay heat (W)",
}


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
def _load_metrics_uq(run_dir: Path, target: str) -> dict:
    with open(run_dir / "metrics.json") as f:
        return json.load(f)["targets"][target]["uq"]


def _load_predictions(run_dir: Path, target: str) -> Optional[pd.DataFrame]:
    p = run_dir / f"predictions_test_{target}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    df["sample_id"] = df["sample_id"].astype(str)
    df["time_years"] = pd.to_numeric(df["time_years"], errors="coerce")
    for c in ["y_true", "y_pred", "y_pred_lo_cal", "y_pred_hi_cal"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _load_per_time(run_dir: Path, target: str) -> pd.DataFrame:
    return pd.read_csv(run_dir / f"metrics_per_time_{target}.csv")


def assign_time_bin(t_yr: float) -> str:
    for (lo, hi), label in zip(zip(BIN_EDGES_YEARS[:-1], BIN_EDGES_YEARS[1:]), BIN_LABELS):
        if t_yr == 0.0 and lo == 0.0:
            return label
        if lo < t_yr <= hi:
            return label
    return BIN_LABELS[-1]


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
    return 5.0  # CFETR default


def per_time_width_from_predictions(run_dir: Path, target: str) -> pd.DataFrame:
    """Compute per-time interval half-width from saved calibrated bounds."""
    eps = _read_eps(run_dir)
    p = pd.read_csv(run_dir / f"predictions_test_{target}.csv")

    yhat = np.maximum(p["y_pred"].astype(float).to_numpy(), 0.0) + eps
    ylo = np.maximum(p["y_pred_lo_cal"].astype(float).to_numpy(), 0.0) + eps
    yhi = np.maximum(p["y_pred_hi_cal"].astype(float).to_numpy(), 0.0) + eps

    width_factor = np.maximum(yhi / yhat, yhat / ylo)
    width_log10 = np.log10(width_factor)

    tmp = pd.DataFrame({
        "time_years": p["time_years"].astype(float).to_numpy(),
        "uq_width_log10": width_log10,
        "uq_half_width_factor": width_factor,
    })
    return (
        tmp.groupby("time_years", as_index=False)
        .agg(
            uq_width_log10=("uq_width_log10", "median"),
            uq_half_width_factor=("uq_half_width_factor", "median"),
        )
        .sort_values("time_years")
    )


def load_per_time_df(run_dir: Path, target: str, method: str, batch_key: str) -> pd.DataFrame:
    """Load per-time metrics and merge with width data."""
    df = _load_per_time(run_dir, target).copy()
    shutdown_years = _read_shutdown_years(run_dir)
    df["time_bin"] = df["time_years"].map(assign_time_bin)
    df["phase"] = df["time_years"].map(
        lambda t: "irradiation" if float(t) <= shutdown_years else "cooling"
    )
    w = per_time_width_from_predictions(run_dir, target)
    df = df.merge(w, on="time_years", how="left")
    # Drop duplicates if floating point precision caused multiple width matches
    df = df.drop_duplicates(subset=["time_years"], keep="first")
    df["method"] = method
    df["batch"] = batch_key
    df["target"] = target
    return df


def per_bin_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-time data into per-bin summaries."""
    w = df["n"].astype(float) if "n" in df.columns else pd.Series(np.ones(len(df)), index=df.index)
    rows = []
    for b in BIN_LABELS:
        sub = df[df["time_bin"] == b]
        if sub.empty:
            continue
        wsub = w.loc[sub.index]
        cov = float(np.average(sub["coverage_90_cal"].astype(float), weights=wsub))
        wl = sub["uq_width_log10"].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
        wf = sub["uq_half_width_factor"].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({
            "time_bin": b,
            "coverage_90_cal": cov,
            "uq_width_log10_median": float(wl.median()) if not wl.empty else float("nan"),
            "uq_half_width_factor_median": float(wf.median()) if not wf.empty else float("nan"),
        })
    return pd.DataFrame(rows)


# ===================================================================
# FIGURE 1: Method Comparison Bar Chart
# ===================================================================
def make_method_comparison_figure() -> None:
    """Generate D224 method comparison bar chart."""
    print("\n=== Fig_P3_method_comparison.pdf ===")

    baseline_csv = TABLES_DIR / "d224_baseline_comparison_results.csv"
    if not baseline_csv.exists():
        print(f"  [SKIP] {baseline_csv} not found")
        return

    b = pd.read_csv(baseline_csv)
    target_order = ["atoms_H3", "atoms_Li6", "atoms_Li7", "decay_heat_W"]

    # CQR data from D224 CQR runs
    cqr_runs = {
        "atoms_H3": REPORT_ROOT / "d224_cqr_timebin_H3",
        "atoms_Li6": REPORT_ROOT / "d224_cqr_timebin_Li6",
        "atoms_Li7": REPORT_ROOT / "d224_cqr_timebin_Li7",
        "decay_heat_W": REPORT_ROOT / "d224_cqr_timebin_DH",
    }

    # Build method data: method_label -> {target -> global_picp}
    methods_data: Dict[str, Dict[str, float]] = {
        "Proposed CV+\n(time-binned)": {},
        "CQR\n(time-binned)": {},
        "ACI\n(time-binned)": {},
        "EnbPI\n(time-binned)": {},
    }

    # Proposed: from baseline CSV "proposed" rows (time_bin only)
    for t in target_order:
        row = b[(b["method"] == "proposed") & (b["target"] == t) & (b["grouping"] == "time_bin")]
        if not row.empty:
            methods_data["Proposed CV+\n(time-binned)"][t] = float(row.iloc[0]["global_picp"])

    # CQR: from metrics.json
    for t, rd in cqr_runs.items():
        mf = rd / "metrics.json"
        if mf.exists():
            with open(mf) as f:
                m = json.load(f)
            tcol = list(m["targets"].keys())[0]
            methods_data["CQR\n(time-binned)"][t] = float(m["targets"][tcol]["uq"]["coverage_90_cal"])

    # ACI / EnbPI: from baseline CSV (time_bin only)
    for method_key, label in [("aci", "ACI\n(time-binned)"), ("enbpi", "EnbPI\n(time-binned)")]:
        for t in target_order:
            row = b[(b["method"] == method_key) & (b["target"] == t) & (b["grouping"] == "time_bin")]
            if not row.empty:
                methods_data[label][t] = float(row.iloc[0]["global_picp"])

    # Plot
    x = np.arange(len(target_order), dtype=float)
    width = 0.18
    method_order = list(methods_data.keys())
    colors = ["#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    fig, ax = plt.subplots(figsize=(10.8, 5.2), constrained_layout=True)

    for i, m_label in enumerate(method_order):
        vals = [methods_data[m_label].get(t, np.nan) for t in target_order]
        offset = (i - (len(method_order) - 1) / 2.0) * width
        bars = ax.bar(x + offset, vals, width=width, label=m_label.replace("\n", " "),
                       color=colors[i], alpha=0.9)
        # Add value labels on top
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + 0.005,
                        f"{val:.1%}", ha="center", va="bottom", fontsize=7.5)

    ax.axhline(0.90, color="k", linestyle="--", linewidth=1.2, label="Nominal 90\\%")
    ax.set_ylim(0.75, 1.08)
    ax.set_ylabel("Global PICP@90\\%")
    ax.set_xticks(x)
    ax.set_xticklabels([TARGET_PRETTY.get(t, t) for t in target_order])
    ax.set_title(r"Method comparison across targets ($\mathcal{D}_{224}$, time-binned calibration)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=9)

    out_pdf = FIGURES_DIR / "Fig_P3_method_comparison.pdf"
    out_png = FIGURES_DIR / "Fig_P3_method_comparison.png"
    fig.savefig(out_pdf, dpi=300)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  WROTE {out_pdf}")
    print(f"  WROTE {out_png}")

    # Also save merged table
    rows = []
    for m_label in method_order:
        for t in target_order:
            v = methods_data[m_label].get(t, np.nan)
            rows.append({"method": m_label.replace("\n", " "), "target": t, "global_picp": v})
    out_csv = TABLES_DIR / "P3_method_comparison_d224_fig.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"  WROTE {out_csv}")


# ===================================================================
# FIGURE 2–3: Per-time and per-bin (D224)
# ===================================================================
def make_per_time_figure() -> None:
    """Generate per-time coverage + width figure for D224."""
    print("\n=== Fig_P3_timeaware_per_time_d224.pdf ===")

    methods = ["global", "time_bin", "phase_time_bin"]
    colors = {"global": "C0", "time_bin": "C1", "phase_time_bin": "C2"}

    frames = []
    for target in TARGETS:
        for method in methods:
            rd = D224_GROUPING_RUNS[target][method]
            if not rd.exists():
                print(f"  [WARN] missing: {rd}")
                continue
            frames.append(load_per_time_df(rd, target, method, "d224"))

    if not frames:
        print("  [SKIP] no data")
        return

    all_df = pd.concat(frames, ignore_index=True)

    fig, axes = plt.subplots(nrows=len(TARGETS), ncols=2, figsize=(12, 9), sharex=True)

    for r, target in enumerate(TARGETS):
        for method in methods:
            sub = all_df[(all_df["target"] == target) & (all_df["method"] == method)].sort_values("time_years")
            if sub.empty:
                continue

            t_years = sub["time_years"].astype(float).to_numpy()
            y_cov = sub["coverage_90_cal"].astype(float).to_numpy()
            y_w = sub["uq_half_width_factor"].astype(float).to_numpy()

            axes[r, 0].plot(t_years, y_cov, label=method if r == 0 else None,
                            color=colors[method], marker=".", markersize=5, linewidth=1.5, alpha=0.8)
            axes[r, 1].plot(t_years, y_w, label=method if r == 0 else None,
                            color=colors[method], marker=".", markersize=5, linewidth=1.5, alpha=0.8)

        axes[r, 0].set_ylabel(f"{TARGET_PRETTY[target]}\nPICP@90")
        axes[r, 0].set_ylim(-0.02, 1.02)
        axes[r, 0].axhline(0.9, color="grey", ls=":", lw=0.6, alpha=0.5)
        axes[r, 0].grid(True, alpha=0.3)

        axes[r, 1].set_ylabel(f"{TARGET_PRETTY[target]}\nUQ half-width")

        # Determine y-scale for width panel
        lines = axes[r, 1].get_lines()
        all_y = [y for line in lines for y in line.get_ydata() if np.isfinite(y)]
        y_min = min(all_y) if all_y else 1.0
        y_max = max(all_y) if all_y else 1.0
        y_range = y_max - y_min

        if y_range < 0.2:
            axes[r, 1].set_yscale("linear")
            def _fmt(y, pos):
                if abs(y - 1.0) < 1e-7:
                    return "1"
                if 0.999 < y < 1.001:
                    diff = y - 1.0
                    return f"1{diff:+.1e}".replace("e-0", "e-").replace("e+0", "e+")
                if 0.9 < y < 1.1:
                    diff = y - 1.0
                    return f"1{diff:+.4f}".rstrip("0")
                return f"{y:g}"

            axes[r, 1].yaxis.set_major_locator(MaxNLocator(integer=False, nbins=5))
            axes[r, 1].yaxis.set_major_formatter(FuncFormatter(_fmt))
            if y_range >= 1e-5:
                axes[r, 1].set_ylim(y_min - y_range * 0.2, y_max + y_range * 0.2)
        else:
            axes[r, 1].set_yscale("log")
            axes[r, 1].yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))

        axes[r, 1].grid(True, alpha=0.3)

    # Panel labels
    letters = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for i, ax in enumerate(axes.flatten()):
        ax.text(0.03, 0.92, letters[i], transform=ax.transAxes, fontweight="bold", fontsize=12, va="top")

    for c, title in enumerate(["Coverage (per time)", "Interval half-width (factor, per time)"]):
        axes[0, c].set_title(title)

    for ax in axes[-1, :]:
        ax.set_xlabel("time (years)")
        ax.set_xscale("log")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        leg = fig.legend(handles, labels, loc="upper center",
                         bbox_to_anchor=(0.5, 0.995), ncol=3, frameon=True)
        leg.get_frame().set_edgecolor("lightgray")
        leg.get_frame().set_linewidth(0.8)
        leg.get_frame().set_facecolor("white")

    # Save — using the OLD filename (n80) so main.tex doesn't need changes
    # We also save a _d224 variant for clarity
    for suffix in ["n80", "d224"]:
        out_pdf = FIGURES_DIR / f"Fig_P3_timeaware_per_time_{suffix}.pdf"
        out_png = FIGURES_DIR / f"Fig_P3_timeaware_per_time_{suffix}.png"
        fig.savefig(out_pdf, dpi=250)
        fig.savefig(out_png, dpi=150)
        print(f"  WROTE {out_pdf}")

    plt.close(fig)

    # Per-bin summary
    bin_rows = []
    for method in methods:
        for target in TARGETS:
            sub = all_df[(all_df["method"] == method) & (all_df["target"] == target)]
            s = per_bin_summary(sub)
            if s.empty:
                continue
            s["method"] = method
            s["target"] = target
            bin_rows.append(s)

    bin_df = pd.concat(bin_rows, ignore_index=True)
    out_csv = TABLES_DIR / "P3_timeaware_per_bin_summary_d224.csv"
    bin_df.to_csv(out_csv, index=False)
    print(f"  WROTE {out_csv}")

    # Per-bin figure
    agg = (
        bin_df.groupby(["method", "time_bin"], as_index=False)
        .agg(
            coverage_90_cal=("coverage_90_cal", "mean"),
            uq_half_width_factor_median=("uq_half_width_factor_median", "mean"),
        )
    )
    _plot_per_bin("d224", agg)


def _plot_per_bin(batch_key: str, bin_df: pd.DataFrame) -> None:
    """Per-bin aggregated figure."""
    print(f"\n=== Fig_P3_timeaware_per_bin_{batch_key}.pdf ===")

    methods = ["global", "time_bin", "phase_time_bin"]
    colors = {"global": "C0", "time_bin": "C1", "phase_time_bin": "C2"}

    bins = [b for b in BIN_LABELS if b in set(bin_df["time_bin"].tolist())]
    x = np.arange(len(bins))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4.2), sharex=True)

    for method in methods:
        sub = bin_df[bin_df["method"] == method].set_index("time_bin")
        if sub.empty:
            continue
        y_cov = [float(sub.loc[b, "coverage_90_cal"]) if b in sub.index else np.nan for b in bins]
        y_w = [float(sub.loc[b, "uq_half_width_factor_median"]) if b in sub.index else np.nan
               for b in bins]

        axes[0].plot(x, y_cov, marker="o", linewidth=2, color=colors[method], label=method)
        axes[1].plot(x, y_w, marker="o", linewidth=2, color=colors[method], label=method)

    axes[0].set_ylabel("PICP@90 (bin-avg)")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].axhline(0.9, color="grey", ls=":", lw=0.6, alpha=0.5)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_ylabel("UQ half-width factor (bin-median)")
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)

    axes[0].set_title("Coverage by time bin")
    axes[1].set_title("Width by time bin")

    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    axes[0].text(0.04, 0.95, "(a)", transform=axes[0].transAxes, fontweight="bold", fontsize=12, va="top")
    axes[1].text(0.04, 0.95, "(b)", transform=axes[1].transAxes, fontweight="bold", fontsize=12, va="top")

    axes[0].set_xticks(x)
    axes[0].set_xticklabels([BIN_DISPLAY.get(b, b) for b in bins], rotation=30, ha="right")

    # Save — overwrite n80 filename and also save _d224 variant
    for suffix in ["n80", "d224"]:
        out_pdf = FIGURES_DIR / f"Fig_P3_timeaware_per_bin_{suffix}.pdf"
        out_png = FIGURES_DIR / f"Fig_P3_timeaware_per_bin_{suffix}.png"
        fig.savefig(out_pdf, dpi=250)
        fig.savefig(out_png, dpi=150)
        print(f"  WROTE {out_pdf}")

    plt.close(fig)


# ===================================================================
# FIGURE 4–6: Trajectory UQ bands (D224)
# ===================================================================
TIME_BIN_EDGES = [1.0 / 365.0, 1.0, 10.0, 100.0, 1000.0, 10000.0]
TIME_BIN_LINELABELS = ["1 d", "1 y", "10 y", "100 y", "1 ky", "10 ky"]

BAND_COLORS = {"global": "#4c72b0", "time_bin": "#dd8452", "phase_time_bin": "#55a868"}
BAND_ALPHAS = {"global": 0.18, "time_bin": 0.25, "phase_time_bin": 0.15}
LINE_STYLES = {"global": "--", "time_bin": "-", "phase_time_bin": ":"}


def _abs_log10_ratio(y_pred: np.ndarray, y_true: np.ndarray) -> np.ndarray:
    ratio = (np.maximum(y_pred, 0.0) + EPS) / (np.maximum(y_true, 0.0) + EPS)
    ratio = np.clip(ratio, 1e-300, 1e300)
    return np.abs(np.log10(ratio))


def _pick_samples(df: pd.DataFrame) -> List[str]:
    """Pick median-error and max-error sample_ids."""
    per_sample = (
        df.groupby("sample_id")
        .apply(
            lambda g: _abs_log10_ratio(g["y_pred"].values, g["y_true"].values).mean(),
            include_groups=False,
        )
        .rename("mean_alr")
        .sort_values()
    )
    ids = list(per_sample.index)
    picked = []
    if ids:
        picked.append(ids[len(ids) // 2])  # median
    if len(ids) > 1:
        picked.append(ids[-1])  # worst
    elif ids:
        picked.append(ids[0])
    return picked


def _plot_single_trajectory(
    ax: plt.Axes,
    sample_id: str,
    target: str,
    runs: Dict[str, pd.DataFrame],
    run_labels: List[str],
    label_suffix: str = "",
) -> None:
    """Plot one sample's trajectory + UQ bands from all grouping runs."""
    first_label = run_labels[0]
    df0 = runs[first_label]
    g0 = df0.loc[df0["sample_id"] == sample_id].sort_values("time_years")
    t = g0["time_years"].values
    y_true = g0["y_true"].values

    ax.scatter(t, y_true, marker="o", s=14, color="black", zorder=10,
               label="Reference (FISPACT-II)", alpha=0.8, linewidths=0)

    for label in run_labels:
        df = runs[label]
        g = df.loc[df["sample_id"] == sample_id].sort_values("time_years")
        if g.empty:
            continue
        tg = g["time_years"].values
        y_pred = g["y_pred"].values
        y_lo = g["y_pred_lo_cal"].values if "y_pred_lo_cal" in g.columns else None
        y_hi = g["y_pred_hi_cal"].values if "y_pred_hi_cal" in g.columns else None

        color = BAND_COLORS.get(label, "#999999")
        ls = LINE_STYLES.get(label, "-")

        ax.plot(tg, y_pred, ls=ls, color=color, lw=1.2, label=f"Pred ({label})")

        if y_lo is not None and y_hi is not None:
            valid = np.isfinite(y_lo) & np.isfinite(y_hi)
            if valid.any():
                alpha = BAND_ALPHAS.get(label, 0.2)
                ax.fill_between(tg[valid], np.maximum(y_lo[valid], 0), y_hi[valid],
                                alpha=alpha, color=color, label=f"90\\% PI ({label})")

    for edge, lbl in zip(TIME_BIN_EDGES, TIME_BIN_LINELABELS):
        if t.min() <= edge <= t.max():
            ax.axvline(edge, color="grey", ls=":", lw=0.5, alpha=0.4)
            ax.text(edge, 0.97, lbl, transform=ax.get_xaxis_transform(),
                    fontsize=5.5, ha="center", va="bottom", color="grey", clip_on=True)

    ax.set_xscale("log")

    # Log-y only if dynamic range is wide enough
    y_pos = y_true[np.isfinite(y_true) & (y_true > 0)]
    use_log_y = y_pos.size >= 2 and (float(np.max(y_pos) / np.min(y_pos)) >= 2.0)

    if use_log_y:
        ax.set_yscale("log")
    else:
        ax.set_yscale("linear")
        yfmt = ScalarFormatter(useMathText=True)
        yfmt.set_scientific(True)
        yfmt.set_powerlimits((-2, 2))
        ax.yaxis.set_major_formatter(yfmt)

    ax.set_xlabel("Time (years)")
    ax.set_ylabel(TARGET_DISPLAY.get(target, target))
    short_id = sample_id.split("::")[-1] if "::" in sample_id else sample_id
    ax.set_title(f"{short_id}{label_suffix}", fontsize=9)


def make_trajectory_figures() -> None:
    """Generate trajectory UQ band figures for H3, Li6, Li7 on D224."""
    print("\n=== Trajectory UQ band figures (D224) ===")

    run_labels = ["global", "time_bin", "phase_time_bin"]

    for target in TARGETS:
        print(f"\n  Target: {target}")

        runs: Dict[str, pd.DataFrame] = {}
        for label in run_labels:
            rd = D224_GROUPING_RUNS[target][label]
            df = _load_predictions(rd, target)
            if df is None:
                print(f"    [WARN] no predictions for {label}: {rd}")
                continue
            runs[label] = df

        if not runs:
            print(f"    [SKIP] no prediction data")
            continue

        # Pick samples from first available run
        first_label = next(iter(runs))
        picked = _pick_samples(runs[first_label])
        if not picked:
            print(f"    [SKIP] no test samples")
            continue

        available_labels = [l for l in run_labels if l in runs]
        n_panels = len(picked)
        fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 4.2), squeeze=False)

        suffixes = [" (median error)", " (worst error)"]
        for i, sid in enumerate(picked):
            suffix = suffixes[i] if i < len(suffixes) else ""
            _plot_single_trajectory(axes[0, i], sid, target, runs, available_labels,
                                    label_suffix=suffix)

        # Shared legend
        handles, labels_ = axes[0, 0].get_legend_handles_labels()
        seen = set()
        unique_h, unique_l = [], []
        for h, lb in zip(handles, labels_):
            if lb not in seen:
                seen.add(lb)
                unique_h.append(h)
                unique_l.append(lb)
        fig.legend(unique_h, unique_l, loc="lower center",
                   ncol=min(len(unique_l), 4), fontsize=7.5, frameon=True,
                   bbox_to_anchor=(0.5, -0.02))

        display = TARGET_DISPLAY.get(target, target)
        fig.suptitle(
            rf"{display}: predicted trajectory with 90\% prediction intervals ($\mathcal{{D}}_{{224}}$)",
            fontsize=11,
        )
        fig.tight_layout(rect=[0, 0.06, 1, 0.95])

        fname = f"Fig_P3_trajectory_uq_bands_{target}"
        for ext, dpi in [(".pdf", 300), (".png", 150)]:
            out = FIGURES_DIR / (fname + ext)
            fig.savefig(out, dpi=dpi, bbox_inches="tight")
            print(f"    WROTE {out}")

        plt.close(fig)


# ===================================================================
# Main
# ===================================================================
def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # Verify all required run dirs exist
    missing = []
    for target, groupings in D224_GROUPING_RUNS.items():
        for g, rd in groupings.items():
            if not rd.exists():
                missing.append(f"  {target}/{g}: {rd}")
    if missing:
        print("ERROR: Missing D224 run directories:")
        print("\n".join(missing))
        return 1

    make_method_comparison_figure()
    make_per_time_figure()
    make_trajectory_figures()

    print("\n=== ALL DONE ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
