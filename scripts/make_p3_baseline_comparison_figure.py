#!/usr/bin/env python3
# ruff: noqa: E402,E501
# flake8: noqa: E402,E501
"""Generate P3 method-comparison figure and table from baseline CSV outputs.

Inputs
------
- tables/baseline_comparison_results.csv   (proposed, aci, enbpi)
- tables/cqr_expanded_targets_results.csv  (cqr)

Outputs
-------
- figures/Fig_P3_method_comparison.pdf
- tables/P3_method_comparison.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _canonical_target_order() -> List[str]:
    return ["atoms_H3", "atoms_Li6", "atoms_Li7", "decay_heat_W"]


def _target_pretty(t: str) -> str:
    m = {
        "atoms_H3": "H-3",
        "atoms_Li6": "Li-6",
        "atoms_Li7": "Li-7",
        "decay_heat_W": "decay_heat_W",
    }
    return m.get(t, t)


def _load_merged_table(
    baseline_csv: Path,
    cqr_csv: Path,
) -> pd.DataFrame:
    b = pd.read_csv(baseline_csv)
    c = pd.read_csv(cqr_csv)

    need = [
        "method",
        "target",
        "grouping",
        "global_picp",
        "tail_picp_1ky_10ky",
        "median_width",
        "mean_width",
    ]
    for col in need:
        if col not in b.columns:
            raise ValueError(f"{baseline_csv}: missing column {col}")
        if col not in c.columns:
            raise ValueError(f"{cqr_csv}: missing column {col}")

    # Build one compact comparison slice used in manuscript figure:
    # - proposed_global
    # - proposed_time_bin
    # - cqr_time_bin
    # - aci_time_bin
    # - enbpi_time_bin
    rows: List[Dict[str, object]] = []

    for t in _canonical_target_order():
        for _, r in b[(b["method"] == "proposed") & (b["target"] == t)].iterrows():
            g = str(r["grouping"])
            if g == "global":
                label = "proposed_global"
            elif g == "time_bin":
                label = "proposed_time_bin"
            else:
                continue
            rows.append(
                {
                    "method_label": label,
                    "method": r["method"],
                    "target": t,
                    "grouping": g,
                    "global_picp": float(r["global_picp"]),
                    "tail_picp_1ky_10ky": float(r["tail_picp_1ky_10ky"]),
                    "median_width": float(r["median_width"]),
                    "mean_width": float(r["mean_width"]),
                    "source": "baseline_comparison_results.csv",
                }
            )

        for m in ["aci", "enbpi"]:
            d = b[(b["method"] == m) & (b["target"] == t) & (b["grouping"] == "time_bin")]
            if d.empty:
                continue
            r = d.iloc[0]
            rows.append(
                {
                    "method_label": m,
                    "method": m,
                    "target": t,
                    "grouping": "time_bin",
                    "global_picp": float(r["global_picp"]),
                    "tail_picp_1ky_10ky": float(r["tail_picp_1ky_10ky"]),
                    "median_width": float(r["median_width"]),
                    "mean_width": float(r["mean_width"]),
                    "source": "baseline_comparison_results.csv",
                }
            )

        d = c[(c["method"] == "cqr") & (c["target"] == t) & (c["grouping"] == "time_bin")]
        if not d.empty:
            r = d.iloc[0]
            rows.append(
                {
                    "method_label": "cqr",
                    "method": "cqr",
                    "target": t,
                    "grouping": "time_bin",
                    "global_picp": float(r["global_picp"]),
                    "tail_picp_1ky_10ky": float(r["tail_picp_1ky_10ky"]),
                    "median_width": float(r["median_width"]),
                    "mean_width": float(r["mean_width"]),
                    "source": "cqr_expanded_targets_results.csv",
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No rows merged for method comparison")

    method_order = [
        "proposed_global",
        "proposed_time_bin",
        "cqr",
        "aci",
        "enbpi",
    ]
    out["_m"] = out["method_label"].map({m: i for i, m in enumerate(method_order)})
    out["_t"] = out["target"].map({t: i for i, t in enumerate(_canonical_target_order())})
    out = out.sort_values(["_m", "_t"]).drop(columns=["_m", "_t"]).reset_index(drop=True)
    return out


def _plot_global_picp_bar(df: pd.DataFrame, out_pdf: Path) -> None:
    target_order = _canonical_target_order()
    method_order = [
        "proposed_global",
        "proposed_time_bin",
        "cqr",
        "aci",
        "enbpi",
    ]

    x = np.arange(len(target_order), dtype=float)
    width = 0.16

    colors = {
        "proposed_global": "#1f77b4",
        "proposed_time_bin": "#2ca02c",
        "cqr": "#ff7f0e",
        "aci": "#9467bd",
        "enbpi": "#8c564b",
    }

    fig, ax = plt.subplots(figsize=(10.8, 5.2), constrained_layout=True)

    for i, m in enumerate(method_order):
        vals: List[float] = []
        for t in target_order:
            d = df[(df["method_label"] == m) & (df["target"] == t)]
            vals.append(float(d["global_picp"].iloc[0]) if not d.empty else np.nan)

        offset = (i - (len(method_order) - 1) / 2.0) * width
        ax.bar(x + offset, vals, width=width, label=m, color=colors[m], alpha=0.9)

    ax.axhline(0.90, color="k", linestyle="--", linewidth=1.2, label="nominal 90%")
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Global PICP@90")
    ax.set_xticks(x)
    ax.set_xticklabels([_target_pretty(t) for t in target_order])
    ax.set_title("Method comparison across targets (PICP@90)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="lower left", ncol=3, frameon=False)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate P3 method comparison figure/table")
    ap.add_argument(
        "--baseline-csv",
        type=Path,
        default=Path("tables/baseline_comparison_results.csv"),
    )
    ap.add_argument(
        "--cqr-csv",
        type=Path,
        default=Path("tables/cqr_expanded_targets_results.csv"),
    )
    ap.add_argument(
        "--out-pdf",
        type=Path,
        default=Path("figures/Fig_P3_method_comparison.pdf"),
    )
    ap.add_argument(
        "--out-table-csv",
        type=Path,
        default=Path("tables/P3_method_comparison.csv"),
    )
    args = ap.parse_args()

    df = _load_merged_table(
        baseline_csv=args.baseline_csv,
        cqr_csv=args.cqr_csv,
    )

    args.out_table_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_table_csv, index=False)

    _plot_global_picp_bar(df, args.out_pdf)

    print(f"[OK] wrote table: {args.out_table_csv} ({len(df)} rows)")
    print(f"[OK] wrote figure: {args.out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
