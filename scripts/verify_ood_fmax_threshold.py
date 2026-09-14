#!/usr/bin/env python3
"""
Verify the OOD gate width-factor statistics quoted in the manuscript.

The manuscript (section "Out-of-distribution threshold rejection") states three numbers:

    (i)   "The 99th percentile of in-distribution f_max across all targets is 8.2,
           conservatively rounded to 10."
    (ii)  "All 44 in-distribution test samples exhibit f_max < 3.2 (observed in the
           sparse >10 ky tail bin for ^55Fe)."
    (iii) "On the OOD batch, all breeding targets exceed f_max = 10 in at least one
           tail bin."

This script recomputes (i) and (ii) from the retained prediction files, and reports the
in-distribution and OOD width-factor distributions side by side so the separation claimed
in the paper can be checked directly.

DEFINITION (recovered from scripts/make_figs5_8_ood_analysis.py, function _per_time_width):

    eps   = the run's `eps` from timeaware_uq_config.json (1.0 for these runs)
    yhat  = max(y_pred,        0) + eps
    ylo   = max(y_pred_lo_cal, 0) + eps
    yhi   = max(y_pred_hi_cal, 0) + eps
    width_factor = max(yhi / yhat, yhat / ylo)

`f_max` for a sample is the maximum of that factor over the sample's time points. The
factor is reference-free and scale-free: it involves only the predicted interval and its
midpoint, so it is available at deployment without the FISPACT-II reference solution.

Usage:
    python3 scripts/verify_ood_fmax_threshold.py [--repo-root .] [--json out.json]
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent

# In-distribution runs: the D224+act9 canonical partition. Targets are spread across two
# locations in this release.
IN_DIST_SOURCES = {
    # 9 activation products (includes 55Fe, named in the manuscript)
    "alltargets_dir": REPO / "data" / "d224_alltargets",
    # breeding isotopes + decay heat
    "phase_timebin_dir": REPO / "data" / "d224_runs",
    "phase_timebin_suffixes": [
        "d224_cvplus_phase_timebin_H3",
        "d224_cvplus_phase_timebin_Li6",
        "d224_cvplus_phase_timebin_Li7",
        "d224_cvplus_decay_heat",
    ],
}

# Out-of-distribution batch (composition/geometry/spectrum shift), for the separation check.
OOD_DIRS = [
    REPO / "data" / "n70_valid" / "baseline",
    REPO / "data" / "n70_valid" / "candidate",
]

BREEDING = ("atoms_H3", "atoms_Li6", "atoms_Li7")


def read_eps(run_dir: Path) -> float:
    cfg = run_dir / "timeaware_uq_config.json"
    if cfg.exists():
        try:
            return float(json.loads(cfg.read_text()).get("eps", 1.0))
        except Exception:
            pass
    return 1.0


def fmax_for_run(run_dir: Path, target: str):
    """Return (f_max per sample as a dict, n_timepoints) or None if unavailable."""
    p = run_dir / f"predictions_test_{target}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    need = {"y_pred", "y_pred_lo_cal", "y_pred_hi_cal"}
    if not need.issubset(df.columns):
        return None
    eps = read_eps(run_dir)
    yhat = np.maximum(df["y_pred"].to_numpy(float), 0.0) + eps
    ylo = np.maximum(df["y_pred_lo_cal"].to_numpy(float), 0.0) + eps
    yhi = np.maximum(df["y_pred_hi_cal"].to_numpy(float), 0.0) + eps
    wf = np.maximum(yhi / yhat, yhat / ylo)
    df = df.assign(_wf=wf)
    # group by sample when a sample id column exists, else treat the whole file as one unit
    key = next((c for c in ("sample_id", "geometry_id", "id") if c in df.columns), None)
    if key is None:
        return {"<all>": float(np.nanmax(wf))}, len(df)
    per = df.groupby(key)["_wf"].max()
    return {str(k): float(v) for k, v in per.items()}, len(df)


def collect_in_distribution():
    out = {}  # target -> {sample: fmax}
    d = IN_DIST_SOURCES["alltargets_dir"]
    if d.is_dir():
        for p in sorted(d.glob("predictions_test_*.csv")):
            tgt = p.stem.replace("predictions_test_", "")
            r = fmax_for_run(d, tgt)
            if r:
                out[tgt] = r[0]
    d2 = IN_DIST_SOURCES["phase_timebin_dir"]
    for sub in IN_DIST_SOURCES["phase_timebin_suffixes"]:
        run = d2 / sub
        if not run.is_dir():
            continue
        for p in sorted(run.glob("predictions_test_*.csv")):
            tgt = p.stem.replace("predictions_test_", "")
            r = fmax_for_run(run, tgt)
            if r:
                out[tgt] = r[0]
    return out


def collect_ood():
    out = {}
    for d in OOD_DIRS:
        if not d.is_dir():
            continue
        for p in sorted(d.glob("predictions_test_*.csv")):
            tgt = p.stem.replace("predictions_test_", "")
            r = fmax_for_run(d, tgt)
            if r:
                out.setdefault(tgt, {}).update(r[0])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--json", default="", help="optional path to write a JSON report")
    args = ap.parse_args(argv)

    ind = collect_in_distribution()
    ood = collect_ood()

    if not ind:
        print("error: no in-distribution prediction files found", file=sys.stderr)
        return 1

    print("=" * 78)
    print("f_max = max over time of max(y_hi/y_hat, y_hat/y_lo), eps=1.0")
    print("=" * 78)

    print("\n--- per-target in-distribution f_max (over the test samples) ---")
    print(f"{'target':16}{'n_samples':>10}{'max_f_max':>12}{'p99':>10}{'median':>10}")
    all_vals = []
    per_target_max = {}
    for tgt in sorted(ind):
        v = np.array(list(ind[tgt].values()), dtype=float)
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        all_vals.append(v)
        per_target_max[tgt] = float(v.max())
        print(
            f"{tgt:16}{v.size:>10}{v.max():>12.4f}{np.percentile(v, 99):>10.4f}{np.median(v):>10.4f}"
        )

    flat = np.concatenate(all_vals)
    print("\n--- manuscript claims vs recomputation (in-distribution) ---")
    print(
        f"  (ii) max in-distribution f_max over ALL samples   = {flat.max():.4f}"
        f"   claim: < 3.2   -> {'MATCH' if flat.max() < 3.2 else 'MISMATCH'}"
    )
    print(
        f"  (i)  99th percentile of in-distribution f_max     = {np.percentile(flat, 99):.4f}"
        f"   claim: 8.2     -> {'MATCH' if abs(np.percentile(flat, 99) - 8.2) <= 0.5 else 'MISMATCH'}"
    )
    print(
        f"       (99th pct over per-target maxima instead    = "
        f"{np.percentile(list(per_target_max.values()), 99):.4f})"
    )

    # which target carries the largest f_max (manuscript names 55Fe for the max)
    worst = max(per_target_max, key=per_target_max.get)
    print(
        f"  worst target by max f_max: {worst} ({per_target_max[worst]:.4f})"
        f"   manuscript names: atoms_Fe55"
    )

    if ood:
        print("\n--- OOD batch (for the separation check) ---")
        print(f"{'target':16}{'n_samples':>10}{'max_f_max':>12}{'>10?':>8}")
        for tgt in sorted(ood):
            v = np.array(list(ood[tgt].values()), dtype=float)
            v = v[np.isfinite(v)]
            if not v.size:
                continue
            print(
                f"{tgt:16}{v.size:>10}{v.max():>12.4f}{('YES' if v.max() > 10 else 'no'):>8}"
            )
        print("\n  (iii) OOD claim is scoped to BREEDING targets in the manuscript:")
        for tgt in BREEDING:
            if tgt in ood:
                v = np.array(list(ood[tgt].values()), dtype=float)
                print(
                    f"       {tgt:16} max f_max = {v.max():.4f}"
                    f"  -> exceeds 10: {'YES' if v.max() > 10 else 'NO'}"
                )

    n_test = max((len(v) for v in ind.values()), default=0)
    print(
        f"\n  (ii) manuscript says 'all 44 in-distribution test samples': "
        f"largest sample count found = {n_test}"
    )

    if args.json:
        Path(args.json).write_text(
            json.dumps(
                {
                    "definition": "f_max = max over time of max(yhi/yhat, yhat/ylo), eps=1.0",
                    "in_distribution": {
                        t: {"n": len(v), "max": max(v.values()) if v else None}
                        for t, v in ind.items()
                    },
                    "ood": {
                        t: {"n": len(v), "max": max(v.values()) if v else None}
                        for t, v in ood.items()
                    },
                    "claim_3p2": float(flat.max()) < 3.2,
                    "claim_8p2_p99": float(np.percentile(flat, 99)),
                },
                indent=2,
            )
        )
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
