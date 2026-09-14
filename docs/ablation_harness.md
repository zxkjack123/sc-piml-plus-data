# Component-ablation harness (recovered and re-run, 2026-09-14)

> **Scope notice.** This document describes work performed in the upstream
> `COOL-PbLi-Burnup` working tree, **not in this data release**. None of `ml/`,
> `datasets/`, the driver, the splits or the run directories referenced below are
> contained in this repository. The only repo artifact from this work is
> `tables/component_ablation_rerun.csv`; everything else needs the upstream tree plus a
> specific interpreter (see §4). Paths written as `ml/...`, `datasets/...` or `reports/...`
> are relative to the upstream repo root.

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
(in the **upstream** repository `reports/` is gitignored and only 6 of `ml/`'s files are
tracked; this data repo is a separate repository). It has therefore been
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

### The re-run does not reproduce the table

This is not a reproduction failure. Three independent checks establish the instrument:

1. The reimplementation reproduces three existing seed-42 runs at **0.00e+00** relative
   deviation (§4).
2. The pass-criterion code reproduces the official `simultaneous_coverage_seed42_max.csv`
   exactly (24/24 rows, §5).
3. **The instrument is calibrated against a number the manuscript itself publishes.** The
   manuscript states "The canonical partition (seed~42, fixed \textit{a priori}) achieves
   4/4 targets passing" (main text, split-robustness section). That is reproduced exactly.

### Both inverted rows hold on verified ground alone

The strongest form of the finding does not depend on the seeds 100/105 splits generated for
this re-run. On the **seed-42 split alone** — the one split whose regeneration provably
matches the retained `splits_d224_seed42.json` — the four configurations give:

| Configuration | cv | aggregation | seed 42 | Manuscript |
|---|---|---|---|---|
| SC-PIML+ (full) | 5 | max | **4/4 PASS** | 4/4 (the §3 calibration anchor) |
| No CV+ (split CP) | 0 | max | **3/4 PASS** | **0/3 splits fail** |
| No aggregation (point-wise) | 5 | none | **0/4 fail** | **2/3 splits pass** |
| 3-fold CV | 3 | max | 1/4 fail | 2/3 splits pass |

`No CV+ (split CP)` **passes** on the verified split while the manuscript reports it as a
catastrophic total failure; `No aggregation` **fails** on the verified split while the
manuscript reports it as performing comparably to MAX. Both inversions are therefore
established without relying on the generated splits.

### The criterion is coverage-only, per the manuscript's own Methods

A challenge worth pre-empting: if the pass rule were joint (coverage *and* width), the wider
max-aggregated intervals would be penalised and the table could be internally consistent.
The manuscript's Methods rules this out. It specifies the aggregation as "**per-geometry
maxima** within each temporal bin", with the quantile "estimated from the set of per-geometry
scores `{S_{i,g}}` **rather than** the full set of per-time-point scores", `n_cal,g ≈ 180`,
and "the finite-sample correction ... ensuring exchangeability at the **geometry level**".
That is exactly what the reimplementation does (§3), and it matches the observed
`n_eff = 180`. The criterion is coverage-only.

### Mechanism

`max` aggregation widens intervals (qhat 0.0948 vs 0.0441 for `cooling|10y-100y`), so
removing it *narrows* intervals and *lowers* simultaneous coverage. That is consistent with
the manuscript's own Figure S3 width argument but contradicts its claim that removing
aggregation leaves performance "comparable to MAX". Note that only `max` widens; `mean`
aggregation narrows relative to `none` (qhat 0.0391 < 0.0441).

A finite-sample explanation was tested and excluded: the quantile-level shift between
`n_eff = 2700` and `180` is `(1+1/180)/(1+1/2700) ≈ 1.005` (about +0.5%), whereas the
observed qhat change is 2.15x (+115%). And `max` and `mean` share `n_eff = 180` while
differing 2.4x in qhat, so the operator rather than the sample count dominates.

### A separate issue found while checking: the 20-split evidence is also thin

The manuscript's `70% (14/20)` reference describes 20 partitions with 149/31/44 sizes on
N=224 (seeds 100-119), naming seeds 100, 101, 103, 111, 112, 114, 116 as the 4/4 cases. The
artifacts matching that description are **not** present: the only split-robustness runs on
disk (`p5_split_robust_r00..r19_{global,time_bin}`) are on **N80** with **51/13/16** splits,
`conformal_cv_folds = 0`, no aggregation knob, and no `phase_time_bin` grouping. The
`d224_stability_max` directory that does correspond to a 20-split stability run retains only
**seed 100's derived tables**. The `14/20` figure therefore cannot be reproduced from
retained artifacts either, and per-seed claims such as "seed 100 achieved 4/4" cannot be
cross-checked. This is noted for the same reason as the ablation discrepancy: it is
checkable in principle, and currently is not checkable in practice.

### Recommendation

Correct the table to the re-run values, and re-frame the surrounding prose so that no stated
contribution depends on the component ranking. The re-run's ranking (aggregation and CV+
configuration matter; fold count does not rescue the row as previously described) differs
from the published one, so the "CV+ is critical" narrative needs revisiting alongside the
numbers. Until that re-framing is agreed, the honest interim position is that the ablation
table's success rates are not supported by the retained artifacts, and that this is
demonstrable on the one split whose identity can be verified rather than inferred.


---

## 8. Follow-up: the 20-split robustness claim was also re-run, and does not reproduce

The manuscript's other simultaneous-coverage result is the 20-split robustness study:
"14/20 splits (70%)" with mean "3.0 +/- 0.8" targets passing, seven named 4/4 seeds
(100, 101, 103, 111, 112, 114, 116), and per-target rates 7Li 18/20, 6Li 17/20, 3H 16/20,
decay heat 10/20. Those artifacts were never retained, so the study was re-run with the
same recovered harness (full configuration: `cv_folds=5`, `aggregation=max`,
`grouping=phase_time_bin`) over seeds 100-119, with the split protocol and pass criterion
both re-validated first (seed-42 split reproduced exactly; criterion reproduced
`simultaneous_coverage_seed42_max.csv` exactly, 24/24 rows).

Result, and every published quantity disagrees:

| Quantity | Published | Re-run |
|---|---|---|
| splits passing >=3/4 | 14/20 | **2/20** (seeds 107, 114) |
| mean targets passing | 3.0 +/- 0.8 | **1.15 +/- 1.01** |
| 4/4 seeds | 100,101,103,111,112,114,116 | **none** |
| 7Li pass rate | 18/20 | **9/20** |
| 6Li pass rate | 17/20 | **10/20** |
| 3H pass rate | 16/20 | **3/20** |
| decay heat pass rate | 10/20 | **1/20** |

The seed-42 anchor does reproduce (4/4), so the instrument is sound; the published claim
looks to be anchored on a favourable split while the other 19 mostly fail.

### The headline is threshold-sensitive, and only matches at a looser threshold

Sweeping the pass threshold over the same 20 runs
(`tables/split_robust_20_threshold_sweep.csv`):

| threshold | splits >=3/4 | mean | 3H / 6Li / 7Li / decay |
|---|---|---|---|
| **0.90 (as documented)** | **2/20** | 1.15 | 3 / 10 / 9 / 1 |
| 0.85 | 7/20 | 1.95 | 8 / 16 / 11 / 4 |
| **0.80** | **14/20** | 2.85 | 11 / 19 / 16 / 11 |
| published | 14/20 | 3.00 | 16 / 17 / 18 / 10 |

So the published 70% is only reached at a threshold about 10 percentage points looser than
the one the manuscript states. At the stated 90% the figure is 10% (2/20). No threshold,
however, reproduces the published per-target rates or the seven named 4/4 seeds, so the
published detail is not reconstructible from this configuration under any threshold tested.

### Caveat on an earlier diagnostic

An earlier automated pass reported that a "geometry-level all-time coverage" variant
reached 14/20 (mean 3.05). I could not reproduce that variant from its description: my
closest reading gives 0/20 (requiring every time point of a geometry to be covered, no
binning) and the marginal per-(geometry,time) reading gives 20/20. That diagnostic should
therefore be treated as **unconfirmed** and is not relied on here. The threshold sweep above
is the reproducible form of the same question.

### Files

- `tables/split_robust_20_seed100_119.csv` - 20 rows, per-seed min bin-wise simultaneous PICP and pass counts.
- `tables/split_robust_20_threshold_sweep.csv` - the threshold sweep above.
- Runs: `reports/wp2_timeaware_uq/d224_seed{100..119}_aggregation_max/`; splits `splits_d224_seed{100..119}.json`.


---

## 9. Independent numerical spot-check of the results tables

A separate independent pass re-derived the numbers from the raw `predictions_test_*.csv`
files. All checks passed; no arithmetic error, transcription error, fabricated value or
cross-file contradiction was found.

| Check | Result |
|---|---|
| All PICP values in [0, 1], no NaN/blank | pass (128 values across both tables) |
| Every value equals an achievable rational k/n_geometries | pass (80/80 once the correct per-bin denominators are used) |
| `targets_passing` equals count of the four columns >= 0.90; `success` equals `>= 3` | pass (20/20 and 12/12 rows) |
| Threshold sweep recomputed from the raw 20-row table | pass (all 7 rows; std agrees to 0.004, consistent with 2-dp storage) |
| Monotonicity: lower threshold cannot reduce passing splits or mean targets | pass (no violation) |
| Cross-file: seeds 100/105 agree between the two tables for the same configuration | pass |
| Seed 42 values and `targets_passing = 4` | pass |
| **Bin-wise criterion reproduces the retained `simultaneous_coverage_seed42_max.csv`** | **pass: 24/24 rows, 0 cell mismatches across 24 x 4 comparisons** |
| End-to-end raw traces, seed 114 (claimed 3/4) and seed 102 (claimed 0/4) | pass, all 8 values identical to 6 dp |
| At threshold 0.90 exactly 2/20 splits pass; at 0.80 exactly 14/20 | pass |

### One assumption that must not be restated in the manuscript

Per-bin geometry counts are **not** uniform. For `atoms_Li6`, `atoms_Li7` and
`decay_heat_W` every bin carries all 44 test geometries, but for **`atoms_H3`** the two tail
bins carry between **3 and 11** geometries depending on the split (measured across seeds
100-119: observed counts 3, 4, 6, 7, 8, 9, 10, 11). This is why some reported values are
not multiples of 1/44 or 1/6 (e.g. seed 103 gives 6/7 = 0.8571, seed 116 gives 7/10 = 0.7000);
both were traced back to the raw predictions and are correct.

Consequence: any statement of the form "the tail bins contain six geometries" is true for
the canonical seed-42 split only, and is **false in general for `atoms_H3`**. Do not restate
it in the manuscript.
