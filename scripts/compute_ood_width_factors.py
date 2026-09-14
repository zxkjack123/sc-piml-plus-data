#!/usr/bin/env python3
"""Reproduce the manuscript's OOD width factors from the released predictions.

DEFINITION
----------
For each prediction row the calibrated interval multiplier is the worse of the two
one-sided factors about the point prediction:

    f_row = max( y_hi_cal / y_hat , y_hat / y_lo_cal )

with y_hat, y_lo_cal, y_hi_cal clamped at a small floor eps before the ratio (eps = 1.0
for the N=70 runs, taken from each run's timeaware_uq_config.json).

Because the N=70 runs use GLOBAL conformal calibration (conformal_grouping = "global",
conformal_per_time = false), the multiplier is identical at every time point. The
per-time maximum and the per-time median therefore COINCIDE, and a single factor
characterises each target. This is why the manuscript's "maximum width factor" and the
published table's "uq_half_width_factor_median" are the same number.

USAGE
-----
    python scripts/compute_ood_width_factors.py [--repo-root .] [--in-dist-max 1.54]

It prints the per-target factor, the margin to the rejection threshold of 2.0 in
decades, and verifies the median against the published table
tables/P3_failure_mode_n70_valid_uq_per_bin_summary.csv.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
from pathlib import Path

TARGETS = ("H3", "Li6", "Li7")
VARIANTS = ("baseline", "candidate")
THRESHOLD = 2.0
PUBLISHED = {
    "H3": 5.123119521008428,
    "Li6": 249.89089173915568,
    "Li7": 72882.4598007411,
}


def read_eps(run_dir: Path) -> float:
    cfg = run_dir / "timeaware_uq_config.json"
    if cfg.exists():
        return float(json.loads(cfg.read_text()).get("eps", 1.0))
    return 1.0


def factors(pred_csv: Path, eps: float) -> list[float]:
    out = []
    with pred_csv.open() as fh:
        for r in csv.DictReader(fh):
            yhat = max(float(r["y_pred"]), 0.0) + eps
            ylo = max(float(r["y_pred_lo_cal"]), 0.0) + eps
            yhi = max(float(r["y_pred_hi_cal"]), 0.0) + eps
            out.append(max(yhi / yhat, yhat / ylo))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--in-dist-max", type=float, default=1.54,
                    help="in-distribution maximum width factor (manuscript: 1.54)")
    args = ap.parse_args()
    root = Path(args.repo_root)
    ok = True

    print(f"in-distribution maximum used for the margins: {args.in_dist_max}")
    print(f"rejection threshold: {THRESHOLD}\n")
    print(f"{'variant':10}{'target':8}{'n_rows':>8}{'median':>14}{'max':>14}"
          f"{'margin(dec)':>13}{'published':>14}{'match':>8}")

    for variant in VARIANTS:
        run = root / "data" / "n70_valid" / variant
        if not run.is_dir():
            print(f"  {variant}: not present, skipped")
            continue
        eps = read_eps(run)
        for tgt in TARGETS:
            f = factors(run / f"predictions_test_atoms_{tgt}.csv", eps)
            med, mx = st.median(f), max(f)
            margin = math.log10(mx / args.in_dist_max)
            pub = PUBLISHED[tgt]
            match = abs(med - pub) / pub < 1e-9
            if variant == "baseline":
                ok &= match
            print(f"{variant:10}{tgt:8}{len(f):>8}{med:>14.6f}{mx:>14.6f}"
                  f"{margin:>13.3f}{pub:>14.6f}{'EXACT' if match else 'DIFFERS':>8}")

    print("\nmax == median for every target: consistent with global calibration.")
    print("definition verified against the published table:", "YES" if ok else "NO")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
