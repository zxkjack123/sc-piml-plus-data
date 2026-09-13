# Component-ablation harness (recovered and re-run, 2026-09-14)

This document records the harness behind the component-ablation table, which was
previously unreproducible: the table's three-split values had no retained artifacts, and
the driver version that produced them was lost. Everything below was recovered from the
runs' own `timeaware_uq_config.json` files and validated by exact reproduction.

---

## 1. Why this was needed

The manuscript's `tab:component-ablation` reports success rates over three splits
(seeds 42, 100, 105) against a "90% simultaneous coverage threshold". Only seed-42
artifacts existed on disk, and no invocation of the ablation was recorded in any script.
An audit found the two knobs that define the rows, but the driver version that implements
one of them (`conformal_geometry_aggregation`) was **not present in any retained code** —
not the repo copy (2026-04-01), not the two staged copies (2026-02-13), not git history
(`reports/` is gitignored and only 6 of `ml/`'s files are tracked). It has therefore been
**reimplemented** and validated by reproduction.

## 2. Harness

| Item | Value |
|---|---|
| Driver | `ml/train_timeaware_uq_curve.py` (133,821 B, 2026-04-01) |
| Reimplementation | `ml/train_timeaware_uq_curve_aggaug.py` (adds `--conformal-geometry-aggregation {none,max,mean}`, default `none`) |
| Data | `datasets/wp2_batches/n224_merged_v1/labels_long_bzsum_elapsed_total_with_yref_logratio.csv` (n_rows = 27328) |
| Splits | `datasets/wp2_batches/n224_merged_v1/splits_d224_seed{42,100,105}.json`, 149/31/44 |
| Targets | `atoms_H3`, `atoms_Li6`, `atoms_Li7`, `decay_heat_W` |
| Fixed config | `--conformal-grouping phase_time_bin --uq-alpha 0.1 --eps 1.0 --uq-method residual_conformal --log-targets atoms_Li6 atoms_Li7 atoms_H3` |

### Split protocol (verified)

`run_timeaware_uq_stability_sweep.py::_make_split`: sort the unique `sample_id` values,
`np.random.RandomState(seed).permutation(ids)`, then `test = perm[:44]`,
`val = perm[44:75]`, `train = perm[75:]`, with the sizes taken from a reference split.

**Verified**: `seed = 42` regenerates the existing `splits_d224_seed42.json` exactly
(149/31/44 over 224 unique sample_ids), which is what licenses the seeds 100/105 splits
generated for this re-run:

| File | sha256[:16] |
|---|---|
| `splits_d224_seed100.json` | `e50b62f9e2f0b10f` |
| `splits_d224_seed105.json` | `e269e9290db69a47` |

## 3. The two knobs

Recovered by diffing each run's `timeaware_uq_config.json`; every other key is identical
across the ablation runs.

| Row | `conformal_cv_folds` | `conformal_geometry_aggregation` |
|---|---|---|
| SC-PIML+ (full) | 5 | max |
| No CV+ (split CP) | 0 | max |
| 3-fold CV | 3 | max |
| No aggregation (point-wise) | 5 | none |

### Semantics of `conformal_geometry_aggregation`

Applied to the per-bin residual-score array, **per geometry, before the quantile**:

- `none` — every (geometry, time) residual is a score; `n_eff = n`
- `max` — one score per geometry, `max` over its residuals in the bin; `n_eff` = #geometries
- `mean` — one score per geometry, `mean` over its residuals in the bin

Evidence (originals, `uq_qhat_by_phase_time_bin_atoms_H3.csv`, group `cooling|10y-100y`,
n = 2700 = 180 geometries × 15 time points):

| mode | n | n_eff | qhat |
|---|---|---|---|
| none | 2700 | 2700 | 0.04405661093370128 |
| max | 2700 | 180 | 0.09483213938538171 |
| mean | 2700 | 180 | 0.03911222500585675 |

Bins with one time point per geometry must coincide across modes — confirmed for
`cooling|1ky-10ky` (n = 34, n_eff = 34, qhat = 4.483864527571902 in all three).

**Consequence**: `max` gives *wider* intervals than `none`, matching the width ratios in
Figure S3.

## 4. Acceptance gate — PASSED

The reimplementation was required to reproduce three existing seed-42 runs before any new
run was made. It reproduced **all 4 targets × 3 runs at 0.00e+00 relative deviation** on
`n`/`n_eff`/`qhat`, `metrics_per_time_*`, `coverage_90_cal` and `median_width_log10_cal`:

- `d224_seed42_aggregation_max` (cv 5, max)
- `d224_seed42_aggregation_none` (cv 5, none)
- `d224_seed42_scpiml_max` (cv 0, max)

An environment note: the runs must be executed with an interpreter carrying the same
sklearn as the originals — `/home/gw/opt/ai_writer/.venv/bin/python` (Python 3.12.3,
sklearn 1.7.2) reproduces them exactly; sklearn 1.7.1 differs at ~1e-13.

## 5. Pass criterion (validated)

A target passes if `min` over the seven temporal bins of the **bin-wise simultaneous**
PICP is ≥ 0.90, where a geometry is "covered" in a bin only if *all* its time points in
that bin are covered. A split succeeds if ≥3 of 4 targets pass.

Validated by reproducing `reports/wp2_timeaware_uq/simultaneous_coverage_seed42_max.csv`
**exactly** — 24/24 rows across all four targets, on `n_geometries`, `n_points`,
`simultaneous_PICP` and `pointwise_PICP`.

Simultaneous PICP must be computed from `predictions_test_<target>.csv`. The driver's
`coverage_by_time_*.csv` is *marginal per time point*, not simultaneous.

## 6. Runs performed

Nine new runs, all local, 11–45 s each, each writing to its own report directory. No
existing directory was overwritten. `cv_folds = 3` passes through unchanged (the
force-to-5 path only fires when `--conformal-method cv_plus`; the default is `split_cp`).

## 7. Results

`tables/component_ablation_rerun.csv` (12 rows = 4 configurations × 3 seeds).

| Configuration | cv folds | aggregation | Re-run success rate | Manuscript |
|---|---|---|---|---|
| SC-PIML+ (full) | 5 | max | **1/3 = 33%** | 70% (14/20)† |
| No CV+ (split CP) | 0 | max | **2/3 = 67%** | **0%** |
| No aggregation (point-wise) | 5 | none | **0/3 = 0%** | **67% (2/3)** |
| 3-fold CV | 3 | max | **0/3 = 0%** | **67% (2/3)** |

† the 20-split robustness reference, a different evaluation.

### The re-run does not reproduce the table, and two rows invert

This is not a reproduction failure — the implementation reproduces the seed-42 runs at
`0.00e+00`, and the pass-criterion code reproduces the official
`simultaneous_coverage_seed42_max.csv` exactly. The pattern is also physically coherent:
`max` aggregation gives wider intervals (qhat 0.0948 vs 0.0441), so removing it *narrows*
intervals and *lowers* simultaneous coverage. That matches the manuscript's own Figure S3
width argument, but contradicts its claim that removing aggregation leaves performance
"comparable to MAX".

A plausible diagnosis for the manuscript's `0%` row: the only seed-42 artifact that fails
outright is `d224_seed42_sc_none`, which turns **both** CV+ **and** aggregation off
(0/4 targets). Attributing that failure to CV+ removal alone confounds the two knobs. The
`67% (2/3)` entries for the other two rows have no producing artifact on disk.

**This needs a ruling before the table is rewritten** — either the table is corrected to
the re-run values (which requires rewriting the surrounding text, since the component
ranking changes), or the metric is revisited.
