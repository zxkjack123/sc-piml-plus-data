# Provenance manifest

Purpose: make every numeric artifact in the manuscript traceable to the file that
produced it, so a reviewer (or a future co-author) can answer "which file backs this
number?" without manual archaeology. Verified 2026-09-13 against manuscript R2.

This file also carries an **exclusion list** of artifacts that must not be used to
derive published numbers.

---

## 1. Verified provenance chains

Each chain was verified by machine comparison of the manuscript table against the source
file, not by inspection.

| Manuscript artifact | Source | Verification | Status |
|---|---|---|---|
| `tab:benchmark` (6 methods x 15 targets) | `tables/benchmark_6methods_unified_summary.csv` | **90/90** PICP cells within 0.05 pp; derived `Pass (>=95%)` row reproduces exactly as 12/11/4/3/13/15 | traceable, **in-repo** |
| `tab:benchmark`, Mean-WIS row | `tables/benchmark_mean_wis_per_target.csv` | **6/6** grand means reproduce: 0.757 / 2.025 / 0.910 / 0.873 / 0.665 / 1.076 vs published 0.76 / 2.02 / 0.91 / 0.87 / 0.67 / 1.08 | traceable, **in-repo** |
| `tab:base-model-summary` (7 learners x 5 metrics) | `tables/base_model_leaderboard.csv` | **35/35** values match | traceable, **in-repo** |
| `tab:error-by-bin` | `data/d224_stability_max/stability_20260721T183813Z/seed_100/table9_error_by_bin.csv` | **24/24** cells (4 targets x 6 bins) | traceable, in-repo |
| `tab:error-stats` | `data/d224_stability_max/stability_20260721T183813Z/seed_100/table10_error_statistics.csv` | **24/24** cells (4 targets x 6 statistics) | traceable, in-repo |
| Supplementary Fig. S1 (`ncal_sensitivity`) | raw: `reports/wp2_timeaware_uq/p5_ncal_sweep_ncal{nc}_{grouping}/metrics.json` -> `tables/ncal_sweep_summary.csv` (30 rows) -> figure via `ml/plot_p3_ncal_sensitivity.py` | chain identified from the generator's `read_csv`/`savefig` calls | traceable |
| Supplementary Fig. S3 / `Fig_S13_binwise_width_comparison` | `tables/fig_s13_binwise_width_data.csv` (28 rows = 4 targets x 7 bins), consumed by `scripts/make_fig_s5_binwise_width.py` | CSV reproduces the previously inlined arrays **exactly**; script now runs end-to-end | traceable, **in-repo** |

### Notes on individual artifacts

**`tab:error-by-bin`** previously collapsed the six source bins into four by merging the
three tail bins. That masked the fact that each tail bin holds a single time point
(n = 1) and that the tritium tail error is 1.7870, not 0.0704. The R2 revision restores
the six-bin structure and adds an `n` row.

**Mean WIS definition** (needed to recompute the row): mean over the 15 targets of the
mean Winkler interval score in `log10(y + 1)` space with `alpha = 0.05`, i.e.
`WIS = (z_hi - z_lo) + (2/alpha) * (z_lo - z_true)_+ + (2/alpha) * (z_true - z_hi)_+`
where `z = log10(max(y,0) + 1.0)`. The per-target values are in
`tables/benchmark_mean_wis_per_target.csv`; grand means match the published row.

**Fig. S3 width ratios**: recomputed from `tables/fig_s13_binwise_width_data.csv`, the
max/no-aggregation width ratio spans **0.52x to 3.00x** across the 17 non-zero
(target, bin) cells, and **7 of those 17 are below unity** (max-aggregation is narrower
there). This confirms the manuscript's stated range and justifies its wording.

---

## 2. Known gaps (reproducibility not yet closed)

| Artifact | Status | Detail |
|---|---|---|
| `tab:component-ablation` | **partial** | The four-row table mixes denominators: the `70% (14/20)` entry is a cross-reference to the 20-split robustness result (fine). The other entries (`0%`, `67% (2/3)`, `67% (2/3)`) describe a **three-split** ablation (seeds 42/100/105 per the section prose), but **only seed-42 artifacts exist on disk** (`reports/wp2_timeaware_uq/d224_seed42_{aggregation_max,aggregation_mean,aggregation_none,sc_none,scpiml_max}`). No producing artifact was found for the `0%` / `2/3` values, and no seeds 100/105 ablation directories exist. The numbers are consistent with the section prose but are **not independently reproducible** from the retained artifacts. |
| `tables/ncal_sweep_summary.csv` | raw source external | The 30-row summary is in-repo, but the raw inputs (`p5_ncal_sweep_*` metrics) live in the external run store. |

**Recommended fix for `tab:component-ablation`:** either record the ablation harness and
re-run the three seeds so the artifacts exist, or restate the table on the seed-42 result
that is retained. Until then the table should be read as reported-but-uncorroborated.

---

## 3. Exclusion list

Do **not** derive any published number from the following.

| Artifact | Reason | Superseded by |
|---|---|---|
| `reports/comprehensive_comparison/run_20260404T134324Z_alpha0.05` (external) | **Degenerate-eps contamination of the SC-PIML column.** Decay heat reports 52.20% where the paper publishes 97.8%; contact dose reports 93.39% where the paper publishes 99.8%. Its SC / CQR / ACI / EnbPI columns match the paper exactly, so only the SC-PIML column is affected. The deviation tracks eps dominance: decay heat -45.6 pp (eps dominates; max 2.5 W, 91% of test points below 1 W), contact dose -6.4 pp (eps comparable at the low end), activity +0.04 pp (eps negligible). Root cause: `log10(y + 1.0)` is effectively linear in y for these low-magnitude aggregates, which collapses the multiplicative conformal score. | `tables/benchmark_6methods_unified_summary.csv` |
| `reports/comprehensive_comparison/run_20260404T153155Z_alpha0.05` (external) | Superseded partial run: matches 86/90 cells. Its **SC-PIML+** column disagrees on four targets (Po-210 96.0 vs 96.4, Hg-203 95.8 vs 97.7, Pb-203 97.3 vs 97.6, Mn-54 95.0 vs 95.9). Its other five columns are correct (75/75). | `tables/benchmark_6methods_unified_summary.csv` |
| `tables/baseline_comparison_results_OLD_bugged_aci.csv` | Already marked as bugged (ACI). | `tables/baseline_comparison_results_alpha005.csv` |

---

## 4. Method-name reconciliation

The manuscript names are SC, CQR, ACI, EnbPI, SC-PIML, SC-PIML+. Intermediate files use
different labels, and one of them is ambiguous.

| Label found in CSVs | Manuscript column | Note |
|---|---|---|
| `proposed`, `Proposed CV+` | **ambiguous** | maps to SC-PIML or SC-PIML+ depending on the file; resolve via section 1, never by label |
| `proposed_cvplus_unprojected` | SC-PIML | |
| `proposed_cvplus_projaware` | SC-PIML+ | |
| `cqr`, `CQR` | CQR | |
| `aci`, `ACI` | ACI | |
| `enbpi`, `EnbPI` | EnbPI | |
| (no label present in older tables) | SC | `tables/benchmark_6methods_unified_summary.csv` is the only in-repo file carrying the explicit `SC` label |

Because `tables/P3_method_comparison_d224.csv` (93.4% for decay heat) and
`tables/alpha0.05_comparison.csv` (97.6%) both use a `Proposed CV+` label for decay heat
but disagree, **the label alone does not identify the manuscript column**.

---

## 5. How to re-verify

```bash
# 1. Error-decomposition tables (fully in-repo)
python3 -c "
import csv
for f in ('table9_error_by_bin','table10_error_statistics'):
    p=f'data/d224_stability_max/stability_20260721T183813Z/seed_100/{f}.csv'
    print(p, len(list(csv.DictReader(open(p)))), 'rows')"

# 2. Mean-WIS row (fully in-repo)
python3 -c "
import csv, collections
rows=list(csv.DictReader(open('tables/benchmark_mean_wis_per_target.csv')))
g=collections.defaultdict(list)
for r in rows: g[r['method']].append(float(r['mean_wis_log10']))
for m in ('SC','CQR','ACI','EnbPI','SC-PIML','SC-PIML+'):
    print(m, round(sum(g[m])/len(g[m]),3))"

# 3. Figure S3 (fully in-repo)
python3 scripts/make_fig_s5_binwise_width.py
```

The benchmark and base-model chains are now fully in-repo (section 1). Only the
`tab:component-ablation` gap in section 2 remains open.
