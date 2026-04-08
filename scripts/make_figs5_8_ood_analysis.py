#!/usr/bin/env python3
"""Regenerate Figures 5-8 with improved presentation.

Fig 5: Projection ablation (^3H) — compact 1×3 layout
Fig 6: OOD failure gate — cleaner bars + outside legend
Fig 7: OOD per-time UQ — step-lines + irradiation shading + clean axes
Fig 8: OOD example curves — shared y-axis + clearer layout
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 150,
    "lines.linewidth": 1.5,
    "lines.markersize": 4,
})
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FixedLocator, NullFormatter, LogLocator

# ── Paths ──
REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "figures"
TABLES_DIR = REPO_ROOT / "tables"

# Projection ablation
ABLATION_RUN_DIR = REPO_ROOT / "data" / "ablation"

# Failure mode / OOD runs
N70_VALID_RUNS: Dict[str, Path] = {
    "baseline": REPO_ROOT / "data" / "n70_valid" / "baseline",
    "candidate": REPO_ROOT / "data" / "n70_valid" / "candidate",
}

TARGETS = ["atoms_H3", "atoms_Li6", "atoms_Li7"]
TARGET_LATEX = {"atoms_H3": "\u00b3H", "atoms_Li6": "\u2076Li", "atoms_Li7": "\u2077Li"}

# Consistent method colors
COLORS = {"Global": "#4C72B0", "Time-binned": "#DD8452", "Phase×time": "#55A868"}
VARIANT_COLORS = {"baseline": "#4C72B0", "candidate": "#DD8452"}
VARIANT_NAMES = {"baseline": "Global", "candidate": "Time-binned"}

SHUTDOWN_YEAR = 10.0

BIN_EDGES = [0.0, 1 / 365, 1, 10, 100, 1000, 10000, float("inf")]
BIN_LABELS = ["≤1d", "1d–1y", "1y–10y", "10y–100y", "100y–1ky", "1ky–10ky", ">10ky"]


def assign_time_bin(t: float) -> str:
    for (lo, hi), label in zip(zip(BIN_EDGES[:-1], BIN_EDGES[1:]), BIN_LABELS):
        if t > lo and t <= hi:
            return label
        if t == 0.0 and lo == 0.0:
            return label
    return BIN_LABELS[-1]


# ═══════════════════════════════════════════════════════════════════════
# Figure 5: Projection Ablation
# ═══════════════════════════════════════════════════════════════════════

def _max_abs_second_difference(z: np.ndarray) -> float:
    z = np.asarray(z, dtype=float)
    if z.size < 3:
        return float("nan")
    d2 = z[:-2] - 2.0 * z[1:-1] + z[2:]
    return float(np.max(np.abs(d2)))


def _load_ablation_mode(run_dir: Path, mode: str):
    mdir = run_dir / f"run__{mode}"
    metrics = json.loads((mdir / "metrics.json").read_text())
    df = pd.read_csv(mdir / "predictions_test_atoms_H3.csv")
    return metrics, df


def _common_d2_max(run_dir: Path) -> float:
    for mode in ["aware", "after"]:
        metrics, _ = _load_ablation_mode(run_dir, mode)
        h3 = metrics["targets"]["atoms_H3"]
        d2_max = h3.get("projection", {}).get("d2_max")
        if d2_max is not None and np.isfinite(float(d2_max)):
            return float(d2_max)
    return float("nan")


def make_fig5() -> None:
    """Projection ablation: compact 1×3 layout with marker+line."""
    print("\n=== Fig 5: Projection Ablation ===")

    d2_max_common = _common_d2_max(ABLATION_RUN_DIR)
    modes = ["none", "after", "aware"]
    mode_labels = ["No projection", "Post-hoc\nprojection", "Projection-\naware"]
    mode_colors = ["#DD8452", "#4C72B0", "#55A868"]

    results = {}
    for mode in modes:
        metrics, df = _load_ablation_mode(ABLATION_RUN_DIR, mode)
        h3 = metrics["targets"]["atoms_H3"]
        eps = float(h3.get("eps", 1.0))

        y = df["y_true"].to_numpy(float)
        yhat = df["y_pred"].to_numpy(float)
        lo = df["y_pred_lo_cal"].to_numpy(float)
        hi = df["y_pred_hi_cal"].to_numpy(float)

        inside = (y >= lo) & (y <= hi)
        coverage = float(np.mean(inside))

        # Width factor
        yhat_e = np.maximum(yhat, 0) + eps
        ylo_e = np.maximum(lo, 0) + eps
        yhi_e = np.maximum(hi, 0) + eps
        w = np.maximum(yhi_e / yhat_e, yhat_e / ylo_e)
        width_log10 = float(np.median(np.log10(np.maximum(w, 1.0))))

        # Smoothness
        z = np.log10(np.maximum(yhat, 0.0) + eps)
        sid = df["sample_id"].astype(str).to_numpy()
        t = df["time_years"].to_numpy(float)
        d2_vals = []
        for s in np.unique(sid):
            m = sid == s
            idx = np.argsort(t[m])
            d2_vals.append(_max_abs_second_difference(z[m][idx]))
        d2a = np.asarray(d2_vals)
        viol_rate = float(np.mean(d2a > d2_max_common + 1e-12)) if np.isfinite(d2_max_common) else float("nan")

        results[mode] = {"coverage": coverage, "width_log10": width_log10, "violation_rate": viol_rate}

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    x = np.arange(len(modes))

    # (a) Coverage
    ax = axes[0]
    vals = [results[m]["coverage"] for m in modes]
    bars = ax.bar(x, vals, color=mode_colors, width=0.6, edgecolor="white", linewidth=0.8)
    ax.axhline(0.90, ls="--", lw=1.0, color="gray", alpha=0.6, label="90% nominal")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0.72, 0.95)
    ax.set_ylabel("PICP (empirical)")
    ax.set_xticks(x)
    ax.set_xticklabels(mode_labels, fontsize=8)
    ax.text(0.02, 0.95, "(a) Coverage", transform=ax.transAxes, fontweight="bold", fontsize=10, va="top")
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)
    ax.legend(fontsize=8, loc="upper right")

    # (b) Interval width
    ax = axes[1]
    vals = [results[m]["width_log10"] for m in modes]
    bars = ax.bar(x, vals, color=mode_colors, width=0.6, edgecolor="white", linewidth=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.0005, f"{v:.4f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("Median log$_{10}$(w)")
    ax.set_xticks(x)
    ax.set_xticklabels(mode_labels, fontsize=8)
    ax.text(0.02, 0.95, "(b) Interval width", transform=ax.transAxes, fontweight="bold", fontsize=10, va="top")
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)
    # Zoom in to show differences
    ymin = min(vals) - 0.002
    ymax = max(vals) + 0.002
    ax.set_ylim(ymin, ymax)

    # (c) Smoothness violations
    ax = axes[2]
    vals = [results[m]["violation_rate"] for m in modes]
    bars = ax.bar(x, vals, color=mode_colors, width=0.6, edgecolor="white", linewidth=0.8)
    for i, v in enumerate(vals):
        label = f"{v*100:.0f}%" if v > 0 else "0%"
        y_pos = v + 0.02 if v > 0 else 0.02
        ax.text(i, y_pos, label, ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("Violation rate")
    ax.set_xticks(x)
    ax.set_xticklabels(mode_labels, fontsize=8)
    ax.text(0.02, 0.98, "(c) Smoothness", transform=ax.transAxes, fontweight="bold", fontsize=10, va="top")
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)

    fig.suptitle("Projection-strategy ablation ($^{3}$H, $D_{224}$, CV+)", fontsize=11)
    fig.subplots_adjust(top=0.88, bottom=0.15, wspace=0.35)

    for ext in ["pdf", "png"]:
        out = FIGURES_DIR / f"Fig_P3_projection_ablation_H3.{ext}"
        fig.savefig(out, dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
        print(f"  WROTE {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Figure 6: OOD Failure Gate
# ═══════════════════════════════════════════════════════════════════════

def _read_gate_json(target: str) -> dict:
    p = N70_VALID_RUNS["candidate"] / f"acceptance_gate_{target}_baseline_vs_minimal_preset.json"
    with p.open("r") as f:
        return json.load(f)


def make_fig6() -> None:
    """OOD failure gate: cleaner bars, outside legend, zoomed y-axis."""
    print("\n=== Fig 6: OOD Failure Gate ===")

    rows = []
    for tgt in TARGETS:
        d = _read_gate_json(tgt)
        tail_name = None
        cand = d.get("slices", {}).get("candidate", {})
        for k in ("tail_window_100_to_min120", "tail_bin_100y_1ky"):
            if k in cand:
                tail_name = k
                break
        if tail_name is None:
            raise KeyError(f"No tail slice found for {tgt}")

        base_tail = d["slices"]["baseline"][tail_name]
        cand_tail = d["slices"]["candidate"][tail_name]

        rows.append({
            "target": tgt,
            "baseline_midterm_cov": float(d["baseline"]["coverage_90_cal_midterm"]),
            "candidate_midterm_cov": float(d["candidate"]["coverage_90_cal_midterm"]),
            "baseline_tail_q99": float(base_tail["q99_factor"]),
            "candidate_tail_q99": float(cand_tail["q99_factor"]),
            "baseline_tail_max": float(base_tail["max_factor"]),
            "candidate_tail_max": float(cand_tail["max_factor"]),
        })

    df = pd.DataFrame(rows)

    fig, (ax_cov, ax_tail) = plt.subplots(1, 2, figsize=(11, 3.8))
    x = np.arange(len(TARGETS))
    w = 0.3
    x_labels = [TARGET_LATEX[t] for t in TARGETS]

    # (a) Midterm coverage
    ax_cov.bar(x - w/2, df["baseline_midterm_cov"], w, color=VARIANT_COLORS["baseline"],
               label="Global", edgecolor="white", linewidth=0.8, alpha=0.85)
    ax_cov.bar(x + w/2, df["candidate_midterm_cov"], w, color=VARIANT_COLORS["candidate"],
               label="Time-binned", edgecolor="white", linewidth=0.8, alpha=0.85)
    ax_cov.axhline(0.90, ls="--", lw=1.0, color="gray", alpha=0.6)
    # Annotate values
    for i, row in df.iterrows():
        ax_cov.text(i - w/2, row["baseline_midterm_cov"] + 0.008, f"{row['baseline_midterm_cov']:.2f}",
                    ha="center", va="bottom", fontsize=7.5, color=VARIANT_COLORS["baseline"])
        ax_cov.text(i + w/2, row["candidate_midterm_cov"] + 0.008, f"{row['candidate_midterm_cov']:.2f}",
                    ha="center", va="bottom", fontsize=7.5, color=VARIANT_COLORS["candidate"])

    ax_cov.set_ylim(0.55, 1.05)
    ax_cov.set_ylabel("Midterm PICP")
    ax_cov.set_xticks(x)
    ax_cov.set_xticklabels(x_labels)
    ax_cov.text(0.02, 0.95, "(a) Midterm coverage", transform=ax_cov.transAxes,
                fontweight="bold", fontsize=10, va="top")
    ax_cov.grid(axis="y", alpha=0.2, linewidth=0.5)

    # (b) Tail-window error factors (log scale)
    bar_w = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5]
    labels_done = set()

    def _bar(ax, x_pos, vals, color, label, hatch=None):
        lbl = label if label not in labels_done else None
        if lbl:
            labels_done.add(label)
        ax.bar(x_pos, vals, bar_w, color=color, label=lbl, edgecolor="white",
               linewidth=0.5, alpha=0.85, hatch=hatch)

    _bar(ax_tail, x + offsets[0] * bar_w, df["baseline_tail_q99"],
         VARIANT_COLORS["baseline"], "Global $q_{99}$")
    _bar(ax_tail, x + offsets[1] * bar_w, df["candidate_tail_q99"],
         VARIANT_COLORS["candidate"], "Time-binned $q_{99}$")
    _bar(ax_tail, x + offsets[2] * bar_w, df["baseline_tail_max"],
         VARIANT_COLORS["baseline"], "Global max", hatch="//")
    _bar(ax_tail, x + offsets[3] * bar_w, df["candidate_tail_max"],
         VARIANT_COLORS["candidate"], "Time-binned max", hatch="//")

    ax_tail.axhline(2.0, ls=":", lw=1.0, color="gray", alpha=0.6, label="Gate threshold (2×)")
    ax_tail.set_yscale("log")
    ax_tail.set_ylabel("Tail error factor")
    ax_tail.set_xticks(x)
    ax_tail.set_xticklabels(x_labels)
    ax_tail.text(0.02, 0.95, "(b) Tail-window error factors", transform=ax_tail.transAxes,
                fontweight="bold", fontsize=10, va="top")
    ax_tail.grid(axis="y", alpha=0.2, linewidth=0.5)

    # Legend at bottom, outside
    handles_a, labels_a = ax_cov.get_legend_handles_labels()
    handles_b, labels_b = ax_tail.get_legend_handles_labels()
    fig.legend(handles_a + handles_b, labels_a + labels_b,
               loc="lower center", bbox_to_anchor=(0.5, -0.12),
               ncol=4, frameon=True, edgecolor="lightgray", fontsize=8,
               fancybox=False, columnspacing=1.2)

    fig.suptitle("Out-of-family failure mode ($\\mathcal{D}_{70,\\mathrm{valid}}$)", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.13)

    for ext in ["pdf", "png"]:
        out = FIGURES_DIR / f"Fig_P3_failure_mode_n70_valid_gate.{ext}"
        fig.savefig(out, dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
        print(f"  WROTE {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Figure 7: OOD Per-Time UQ Diagnostics
# ═══════════════════════════════════════════════════════════════════════

def _read_eps(run_dir: Path) -> float:
    cfg_path = run_dir / "timeaware_uq_config.json"
    if cfg_path.exists():
        cfg = json.loads(cfg_path.read_text())
        return float(cfg.get("eps", 1.0))
    return 1.0


def _per_time_width(run_dir: Path, target: str) -> pd.DataFrame:
    eps = _read_eps(run_dir)
    p = pd.read_csv(run_dir / f"predictions_test_{target}.csv")

    yhat = np.maximum(p["y_pred"].astype(float).to_numpy(), 0.0) + eps
    ylo = np.maximum(p["y_pred_lo_cal"].astype(float).to_numpy(), 0.0) + eps
    yhi = np.maximum(p["y_pred_hi_cal"].astype(float).to_numpy(), 0.0) + eps

    wf = np.maximum(yhi / yhat, yhat / ylo)
    tmp = pd.DataFrame({"time_years": p["time_years"].astype(float), "uq_half_width_factor": wf})
    return tmp.groupby("time_years", as_index=False).agg(
        uq_half_width_factor=("uq_half_width_factor", "median")
    ).sort_values("time_years")


def _load_per_time(variant: str, target: str) -> pd.DataFrame:
    run_dir = N70_VALID_RUNS[variant]
    df = pd.read_csv(run_dir / f"metrics_per_time_{target}.csv")
    w = _per_time_width(run_dir, target)
    df = df.merge(w, on="time_years", how="left")
    df["variant"] = variant
    df["target"] = target
    return df


def make_fig7() -> None:
    """OOD per-time UQ: step-lines, irradiation shading, clean axes."""
    print("\n=== Fig 7: OOD Per-Time UQ ===")

    frames = []
    for v in ["baseline", "candidate"]:
        for tgt in TARGETS:
            frames.append(_load_per_time(v, tgt))
    all_df = pd.concat(frames, ignore_index=True)

    fig, axes = plt.subplots(len(TARGETS), 2, figsize=(12, 8), sharex=True)

    panel_idx = 0
    for r, tgt in enumerate(TARGETS):
        t_label = TARGET_LATEX[tgt]

        for v in ["baseline", "candidate"]:
            sub = all_df[(all_df["target"] == tgt) & (all_df["variant"] == v)].sort_values("time_years")
            if sub.empty:
                continue
            x = sub["time_years"].astype(float).to_numpy()
            y_cov = sub["coverage_90_cal"].astype(float).to_numpy()
            y_w = sub["uq_half_width_factor"].astype(float).to_numpy()
            dname = VARIANT_NAMES[v]
            color = VARIANT_COLORS[v]

            # Coverage (left column)
            ax = axes[r, 0]
            ax.step(x, y_cov, where="mid", color=color, linewidth=1.5,
                    label=dname if r == 0 else None)
            ax.plot(x[::3], y_cov[::3], "o", color=color, markersize=3, alpha=0.6)

            # Width (right column)
            ax = axes[r, 1]
            ax.step(x, y_w, where="mid", color=color, linewidth=1.5,
                    label=dname if r == 0 else None)
            ax.plot(x[::3], y_w[::3], "o", color=color, markersize=3, alpha=0.6)

        # Coverage panel decoration
        ax_c = axes[r, 0]
        ax_c.axhline(0.90, ls="--", lw=0.8, color="gray", alpha=0.6)
        ax_c.set_ylim(0.35, 1.05)
        ax_c.set_ylabel(f"{t_label} PICP")
        ax_c.grid(axis="y", alpha=0.2, linewidth=0.5)
        label_c = chr(97 + panel_idx)
        ax_c.text(0.02, 0.95, f"({label_c})", transform=ax_c.transAxes,
                  fontweight="bold", fontsize=10, va="top")
        panel_idx += 1

        # Irradiation shading
        xlim = ax_c.get_xlim()
        ax_c.axvspan(xlim[0], SHUTDOWN_YEAR, color="lightgray", alpha=0.15, zorder=0)
        ax_c.axvline(SHUTDOWN_YEAR, ls=":", lw=0.8, color="gray", alpha=0.5)

        # Width panel decoration
        ax_w = axes[r, 1]
        # Collect plotted width data for this target to decide scale
        sub_all = all_df[all_df["target"] == tgt]
        w_vals = sub_all["uq_half_width_factor"].dropna().astype(float)
        w_min, w_max = float(w_vals.min()), float(w_vals.max())
        w_ratio = w_max / w_min if w_min > 0 else 1.0

        if w_ratio < 2.0:
            # Near-constant data: use linear scale with padded limits
            ax_w.set_yscale("linear")
            pad = max(abs(w_max) * 0.1, 1.0)
            ax_w.set_ylim(w_min - pad, w_max + pad)
        else:
            ax_w.set_yscale("log")
            ax_w.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0), numticks=12))
            sf = ScalarFormatter()
            sf.set_scientific(False)
            sf.set_useOffset(False)
            ax_w.yaxis.set_major_formatter(sf)
            ax_w.yaxis.set_minor_formatter(NullFormatter())
        ax_w.set_ylabel(f"{t_label} width factor")
        ax_w.grid(axis="y", alpha=0.2, linewidth=0.5)
        label_w = chr(97 + panel_idx)
        ax_w.text(0.02, 0.95, f"({label_w})", transform=ax_w.transAxes,
                  fontweight="bold", fontsize=10, va="top")
        panel_idx += 1

        # Irradiation shading on width panel
        xlim_w = ax_w.get_xlim()
        ax_w.axvspan(xlim_w[0], SHUTDOWN_YEAR, color="lightgray", alpha=0.15, zorder=0)
        ax_w.axvline(SHUTDOWN_YEAR, ls=":", lw=0.8, color="gray", alpha=0.5)

    # X-axis labels and scales
    for c in range(2):
        axes[-1, c].set_xlabel("Time (years)")
        for r in range(len(TARGETS)):
            axes[r, c].set_xscale("log")

    # Column titles
    axes[0, 0].set_title("Coverage (per time point)", fontsize=10)
    axes[0, 1].set_title("Interval width factor (per time point)", fontsize=10)

    # Legend at bottom
    handles, labels = axes[0, 0].get_legend_handles_labels()
    # Add reference lines
    from matplotlib.lines import Line2D
    handles.append(Line2D([0], [0], ls="--", lw=0.8, color="gray", alpha=0.6))
    labels.append("90% nominal")
    handles.append(Line2D([0], [0], ls=":", lw=0.8, color="gray", alpha=0.5))
    labels.append("Shutdown")

    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=4, frameon=True, edgecolor="lightgray", fontsize=9, fancybox=False)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.07)

    for ext in ["pdf", "png"]:
        out = FIGURES_DIR / f"Fig_P3_failure_mode_n70_valid_uq.{ext}"
        fig.savefig(out, dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
        print(f"  WROTE {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Figure 8: OOD Example Curves
# ═══════════════════════════════════════════════════════════════════════

def _find_worst_samples(candidate_dir: Path, target: str, top_k: int = 2) -> List[str]:
    """Find worst-tail OOD samples by max prediction error factor."""
    eps = _read_eps(candidate_dir)
    p = pd.read_csv(candidate_dir / f"predictions_test_{target}.csv")
    p["time_years"] = p["time_years"].astype(float)

    # Focus on tail window
    tail = p[(p["time_years"] >= 100.0) & (p["time_years"] <= 120.0)].copy()
    if tail.empty:
        tail = p[p["time_years"] >= 50.0].copy()
    if tail.empty:
        return p["sample_id"].unique()[:top_k].tolist()

    y = np.maximum(tail["y_true"].astype(float).to_numpy(), 0) + eps
    yhat = np.maximum(tail["y_pred"].astype(float).to_numpy(), 0) + eps
    tail["factor"] = np.power(10, np.abs(np.log10(yhat / y)))

    worst = (tail.groupby("sample_id")["factor"].max()
             .sort_values(ascending=False).head(top_k))
    return worst.index.tolist()


def make_fig8() -> None:
    """OOD example curves: 2×2 with shared structure, clearer colors."""
    print("\n=== Fig 8: OOD Example Curves ===")

    candidate_dir = N70_VALID_RUNS["candidate"]
    baseline_dir = N70_VALID_RUNS["baseline"]
    target = "atoms_H3"

    sample_ids = _find_worst_samples(candidate_dir, target, top_k=2)
    if not sample_ids:
        print("  No samples found, skipping Fig 8")
        return

    def _load(run_dir: Path) -> pd.DataFrame:
        return pd.read_csv(run_dir / f"predictions_test_{target}.csv")

    p_base = _load(baseline_dir)
    p_cand = _load(candidate_dir)

    n = len(sample_ids)
    fig, axes = plt.subplots(n, 2, figsize=(11, 4.0 * n), sharey=False)
    if n == 1:
        axes = np.array([axes])

    variant_data = [
        ("Global", p_base, VARIANT_COLORS["baseline"]),
        ("Time-binned", p_cand, VARIANT_COLORS["candidate"]),
    ]

    for i, sid in enumerate(sample_ids):
        for j, (vname, p, color) in enumerate(variant_data):
            sub = p[p["sample_id"] == sid].sort_values("time_years")
            if sub.empty:
                continue

            x = sub["time_years"].astype(float).to_numpy()
            y = sub["y_true"].astype(float).to_numpy()
            yhat = sub["y_pred"].astype(float).to_numpy()
            lo = sub["y_pred_lo_cal"].astype(float).to_numpy()
            hi = sub["y_pred_hi_cal"].astype(float).to_numpy()

            ax = axes[i, j]

            # Reference (truth)
            ax.plot(x, y, color="#222222", linewidth=2, label="Reference", zorder=3)

            # Prediction
            ax.plot(x, yhat, color=color, linewidth=1.8, label="Prediction", zorder=3)

            # Prediction interval
            ax.fill_between(x, lo, hi, color=color, alpha=0.18,
                            label="90% prediction interval", zorder=1)

            # Irradiation shading
            ax.axvspan(x.min(), SHUTDOWN_YEAR, color="lightgray", alpha=0.15, zorder=0)
            ax.axvline(SHUTDOWN_YEAR, ls=":", lw=0.8, color="gray", alpha=0.5)

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, alpha=0.15, linewidth=0.5)
            ax.set_xlabel("Time (years)")

            if j == 0:
                ax.set_ylabel("$^{3}$H inventory (atoms)")

            label_char = chr(97 + i * 2 + j)
            ax.text(0.02, 0.95, f"({label_char})", transform=ax.transAxes,
                    fontweight="bold", fontsize=10, va="top",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.8, edgecolor="none"))

            ax.set_title(f"{vname} — sample {sid}", fontsize=9)

    # Shared legend at bottom
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=True, edgecolor="lightgray", fontsize=9, fancybox=False)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.06)

    for ext in ["pdf", "png"]:
        out = FIGURES_DIR / f"Fig_P3_failure_mode_n70_valid_examples_H3.{ext}"
        fig.savefig(out, dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
        print(f"  WROTE {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    make_fig5()
    make_fig6()
    make_fig7()
    make_fig8()
    print("\nAll done.")
