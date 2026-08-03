#!/usr/bin/env python3
r"""Paper 3 — projection-aware conformal ablation summary (table + figure).

This script consumes a *single* ablation output directory produced by
`ml/run_p3_projection_aware_ablation.py` (or the SugonHB wrapper) and writes:

- a paper-ready CSV table summarizing calibration + physics-consistency metrics
- a compact figure comparing none/after/aware

By default it is pinned to the current ablation evidence directory under the
SSOT (`/home/gw/ComputeData/CFETR/COOL-PbLi-Burnup/reports/wp2_timeaware_uq`).

Notes
-----
- We focus on `atoms_H3` for the physics projection (smoothness constraint in
  model space) because H3 is log-transformed and has the most physics/shape
  sensitivity.
- Smoothness is measured by max absolute second difference in model space:

  $d2(z) = \max_i |z_{i} - 2 z_{i+1} + z_{i+2}|$.

  We report the per-sample violation rate under a *common* threshold $d2_{max}$
  taken from the projected runs (after/aware).

Outputs
-------
- figures/Fig_P3_projection_ablation_H3.(svg|png)
- tables/P3_projection_ablation_H3_summary.csv

"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
# NOTE (R1 revision): This script reads the original ablation experiment outputs
# from a local working directory (COOL-PbLi-Burnup reports tree). The final
# ablation results used in the manuscript (Table 7) are frozen in
# tables/P3_projection_ablation_H3_summary.csv; the processed per-run outputs
# are archived under data/ablation/ (run__after / run__aware / run__none).

from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


DEFAULT_ABLATION_RUN_DIR = Path(
    "/mnt/nas_smb/CrossDevice/ASIPPWORK/项目课题/2023-基于机器学习的聚变堆增殖包层燃耗演变快速预测研究-院长基金/github/coolburnup-bench-fed/reports/wp2_timeaware_uq/"
    "n80_deepsmoke_addon_v1_addon2_mat_v3_absct_recalc_20260118/"
    "p3_proj_ablation_20260213T032109Z_j23876165"
)


@dataclass(frozen=True)
class ModeSummary:
    mode: str
    n_rows: int
    n_samples: int
    coverage_empirical: float
    coverage_midterm: float
    width_log10_median: float
    width_log10_mean: float
    qhat_global: float
    pred_neg_rate: float
    d2_max_common: float
    d2_violation_rate: float
    d2_median: float
    d2_p95: float
    d2_max_observed: float


def _max_abs_second_difference(z: np.ndarray) -> float:
    z = np.asarray(z, dtype=float)
    if z.size < 3:
        return float("nan")
    d2 = z[:-2] - 2.0 * z[1:-1] + z[2:]
    return float(np.max(np.abs(d2)))


def _width_log10_half_factor(
    *, yhat: np.ndarray, lo: np.ndarray, hi: np.ndarray, eps: float
) -> np.ndarray:
    """Paper-consistent half-width factor in log10 units.

    w = max( (U+eps)/(yhat+eps), (yhat+eps)/(L+eps) )
    width_log10 = log10(w)
    """

    yhat = np.asarray(yhat, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)

    num1 = (np.maximum(hi, 0.0) + eps) / (np.maximum(yhat, 0.0) + eps)
    num2 = (np.maximum(yhat, 0.0) + eps) / (np.maximum(lo, 0.0) + eps)
    w = np.maximum(num1, num2)
    return np.log10(np.maximum(w, 1.0))


def _load_h3_mode_dir(run_dir: Path, mode: str) -> Tuple[dict, pd.DataFrame]:
    mdir = run_dir / f"run__{mode}"
    metrics = json.loads((mdir / "metrics.json").read_text())
    df = pd.read_csv(mdir / "predictions_test_atoms_H3.csv")
    return metrics, df


def _common_d2_max(run_dir: Path) -> float:
    """Choose a common d2_max threshold for fair violation comparison."""

    for mode in ["aware", "after"]:
        metrics, _df = _load_h3_mode_dir(run_dir, mode)
        h3 = metrics["targets"]["atoms_H3"]
        d2_max = h3.get("projection", {}).get("d2_max")
        if d2_max is not None and np.isfinite(float(d2_max)):
            return float(d2_max)
    return float("nan")


def _summarize_mode(run_dir: Path, mode: str, *, d2_max_common: float) -> ModeSummary:
    metrics, df = _load_h3_mode_dir(run_dir, mode)

    h3 = metrics["targets"]["atoms_H3"]
    eps = float(h3.get("eps", 1.0))

    y = df["y_true"].to_numpy(float)
    yhat = df["y_pred"].to_numpy(float)
    lo = df["y_pred_lo_cal"].to_numpy(float)
    hi = df["y_pred_hi_cal"].to_numpy(float)

    inside = (y >= lo) & (y <= hi)
    coverage_emp = float(np.mean(inside))

    cov_mid = h3.get("uq", {}).get("coverage_90_cal_midterm")
    coverage_midterm = float(cov_mid) if cov_mid is not None else float("nan")

    width_log10 = _width_log10_half_factor(yhat=yhat, lo=lo, hi=hi, eps=eps)

    qhat = h3.get("uq", {}).get("qhat_conformal_global")
    qhat_global = float(qhat) if qhat is not None else float("nan")

    pred_neg_rate = float(np.mean(yhat < 0.0))

    # Smoothness diagnostics in model space z=log10(y_pred+eps)
    z = np.log10(np.maximum(yhat, 0.0) + eps)
    sid = df["sample_id"].astype(str).to_numpy()
    t = df["time_years"].to_numpy(float)

    d2_vals: List[float] = []
    for s in np.unique(sid):
        m = sid == s
        idx = np.argsort(t[m])
        d2_vals.append(_max_abs_second_difference(z[m][idx]))
    d2a = np.asarray(d2_vals, dtype=float)

    if np.isfinite(d2_max_common):
        viol = np.isfinite(d2a) & (d2a > float(d2_max_common) + 1e-12)
        d2_viol_rate = float(np.mean(viol))
    else:
        d2_viol_rate = float("nan")

    return ModeSummary(
        mode=mode,
        n_rows=int(len(df)),
        n_samples=int(len(d2a)),
        coverage_empirical=coverage_emp,
        coverage_midterm=coverage_midterm,
        width_log10_median=float(np.nanmedian(width_log10)),
        width_log10_mean=float(np.nanmean(width_log10)),
        qhat_global=qhat_global,
        pred_neg_rate=pred_neg_rate,
        d2_max_common=float(d2_max_common),
        d2_violation_rate=d2_viol_rate,
        d2_median=float(np.nanmedian(d2a)),
        d2_p95=float(np.nanquantile(d2a, 0.95)),
        d2_max_observed=float(np.nanmax(d2a)),
    )


def _write_table(out_csv: Path, summaries: Iterable[ModeSummary]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for s in summaries:
        rows.append(
            {
                "mode": s.mode,
                "n_rows": s.n_rows,
                "n_samples": s.n_samples,
                "coverage_empirical": s.coverage_empirical,
                "coverage_90_cal_midterm": s.coverage_midterm,
                "width_log10_median": s.width_log10_median,
                "width_log10_mean": s.width_log10_mean,
                "qhat_conformal_global": s.qhat_global,
                "pred_negative_rate": s.pred_neg_rate,
                "d2_max_common": s.d2_max_common,
                "d2_violation_rate": s.d2_violation_rate,
                "d2_median": s.d2_median,
                "d2_p95": s.d2_p95,
                "d2_max_observed": s.d2_max_observed,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


def _make_figure(out_svg: Path, out_pdf: Path, table: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    modes = ["none", "after", "aware"]
    table = table.set_index("mode").loc[modes].reset_index()

    fig, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(7.2, 7.4),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0]},
    )

    x = np.arange(len(modes))

    # Coverage
    ax = axes[0]
    ax.plot(x, table["coverage_empirical"], marker="o", lw=2)
    ax.axhline(0.9, ls="--", lw=1.3, color="0.4", label="target 90%")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("PICP@90 (empirical)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    ax.text(
        0.02,
        0.9,
        "(a)",
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=12,
        va="top",
    )

    # Width (paper-consistent half-width factor)
    ax = axes[1]
    ax.plot(x, table["width_log10_median"], marker="o", lw=2)
    ax.set_ylabel(r"median $\log_{10}(w)$")
    ax.grid(True, alpha=0.3)
    ax.text(
        0.02,
        0.9,
        "(b)",
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=12,
        va="top",
    )

    # Smoothness violations (bounded second differences)
    ax = axes[2]
    ax.bar(x, table["d2_violation_rate"], color=["#d55e00", "#0072b2", "#009e73"])
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("smoothness violation rate")
    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(
        0.02,
        0.9,
        "(c)",
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=12,
        va="top",
    )

    fig.suptitle("Projection-aware conformal ablation ($^3$H)")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg)

    fig.savefig(out_pdf, dpi=200)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Generate Paper 3 projection-aware ablation summary (H3)."
    )
    ap.add_argument(
        "--ablation-run-dir",
        type=Path,
        default=DEFAULT_ABLATION_RUN_DIR,
        help="Directory containing run__none/run__after/run__aware",
    )
    ap.add_argument(
        "--out-fig-svg",
        type=Path,
        default=Path("../figures/Fig_P3_projection_ablation_H3.svg"),
    )
    ap.add_argument(
        "--out-fig-png",
        type=Path,
        default=Path("../figures/Fig_P3_projection_ablation_H3.pdf"),
    )
    ap.add_argument(
        "--out-table-csv",
        type=Path,
        default=Path("../tables/P3_projection_ablation_H3_summary.csv"),
    )
    args = ap.parse_args()

    run_dir = args.ablation_run_dir
    if not run_dir.is_dir():
        raise SystemExit(f"ablation-run-dir not found: {run_dir}")

    d2_max_common = _common_d2_max(run_dir)

    summaries = [
        _summarize_mode(run_dir, "none", d2_max_common=d2_max_common),
        _summarize_mode(run_dir, "after", d2_max_common=d2_max_common),
        _summarize_mode(run_dir, "aware", d2_max_common=d2_max_common),
    ]

    out_csv = (Path(__file__).resolve().parent / args.out_table_csv).resolve()
    out_svg = (Path(__file__).resolve().parent / args.out_fig_svg).resolve()
    out_png = (Path(__file__).resolve().parent / args.out_fig_png).resolve()

    table = _write_table(out_csv, summaries)
    _make_figure(out_svg, out_png.with_suffix(".pdf"), table)

    print("WROTE")
    print(" table:", out_csv)
    print(" fig:", out_svg)
    print(" fig:", out_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
