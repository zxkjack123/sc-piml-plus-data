#!/usr/bin/env python3
"""
Build tables/ncal_sweep_summary.csv from the raw p5_ncal_sweep run store.

This script reconstructs the intermediate that backs Supplementary Fig. S1
(Fig_P3_ncal_sensitivity.pdf). It was previously produced by an unretained
ad-hoc step, so the CSV existed in the release with no script behind it; this
closes that gap.

Chain:
    <run_store>/p5_ncal_sweep_ncal{nc}_{grouping}/metrics.json
        -> tables/ncal_sweep_summary.csv
        -> Fig_P3_ncal_sensitivity.pdf

Source fields (per target, inside each metrics.json):
    PICP                 <- targets[<t>]["uq"]["coverage_90_cal"]
    median_width_log10   <- targets[<t>]["uq"]["median_width_log10_cal"]
    n_train              <- len(splits["train"])

Usage:
    python3 scripts/make_ncal_sweep_summary.py \
        --run-store /path/to/reports/wp2_timeaware_uq/<batch_id> \
        [--out tables/ncal_sweep_summary.csv]
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parent / "tables" / "ncal_sweep_summary.csv"

DIR_RE = re.compile(r"^p5_ncal_sweep_ncal(?P<ncal>\d+)_(?P<grouping>.+)$")

FIELDNAMES = ["ncal", "grouping", "target", "PICP", "median_width_log10", "n_train"]


def collect(run_store: Path):
    rows = []
    for d in sorted(run_store.glob("p5_ncal_sweep_ncal*_*")):
        m = DIR_RE.match(d.name)
        if not m or not (d / "metrics.json").is_file():
            continue
        ncal = int(m.group("ncal"))
        grouping = m.group("grouping")
        with open(d / "metrics.json") as fh:
            md = json.load(fh)
        n_train = len(md.get("splits", {}).get("train", []))
        for target, tinfo in md.get("targets", {}).items():
            uq = tinfo.get("uq") or {}
            picp = uq.get("coverage_90_cal")
            width = uq.get("median_width_log10_cal")
            if picp is None or width is None:
                continue
            rows.append(
                {
                    "ncal": ncal,
                    "grouping": grouping,
                    "target": target,
                    "PICP": f"{picp:.4f}",
                    "median_width_log10": f"{width:.6e}",
                    "n_train": n_train,
                }
            )
    rows.sort(key=lambda r: (r["ncal"], r["grouping"], r["target"]))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument(
        "--run-store",
        required=True,
        type=Path,
        help="Directory holding the p5_ncal_sweep_ncal*_* run directories.",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    if not args.run_store.is_dir():
        print(f"error: run store not found: {args.run_store}", file=sys.stderr)
        return 1

    rows = collect(args.run_store)
    if not rows:
        print(
            f"error: no p5_ncal_sweep_* runs found under {args.run_store}",
            file=sys.stderr,
        )
        return 1

    # Guard the exact grid the released CSV encodes. Without this, a run store missing a
    # target (or a metrics.json lacking `splits`) would silently emit a short file or an
    # n_train of 0, defeating the point of a byte-reproducible artifact.
    expected = {
        (n, g, t)
        for n in (13, 16, 20, 26, 32)
        for g in ("global", "time_bin")
        for t in ("atoms_H3", "atoms_Li6", "atoms_Li7")
    }
    got = {(r["ncal"], r["grouping"], r["target"]) for r in rows}
    if missing := expected - got:
        print(
            f"error: missing {len(missing)} (ncal, grouping, target) rows: "
            f"{sorted(missing)}",
            file=sys.stderr,
        )
        return 1
    if len(rows) != 30:
        print(f"error: expected 30 rows, got {len(rows)}", file=sys.stderr)
        return 1
    zero_train = [r for r in rows if int(r["n_train"]) == 0]
    if zero_train:
        print(
            f"error: {len(zero_train)} rows have n_train=0 (splits block missing?)",
            file=sys.stderr,
        )
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDNAMES, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    ncal_vals = sorted({r["ncal"] for r in rows})
    groupings = sorted({r["grouping"] for r in rows})
    print(f"wrote {len(rows)} rows to {args.out}")
    print(f"  ncal in {ncal_vals}, groupings {groupings}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
