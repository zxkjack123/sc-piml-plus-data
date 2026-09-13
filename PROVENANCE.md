# Provenance manifest

Purpose: make every numeric artifact in the manuscript traceable to the file that
produced it, so a reviewer (or a future co-author) can answer "which file backs this
number?" without manual archaeology. Verified 2026-09-13 against manuscript R2.

This file also carries an **exclusion list** of artifacts that must not be used to
derive published numbers.

---

## 1. Verified provenance chains

Each chain below was verified by machine comparison of the manuscript table against the
source file, not by inspection.

| Manuscript artifact | Source | Verification | Verdict |
|---|---|---|---|
| `tab:benchmark` (6 methods x 15 targets, PICP + derived pass counts) | `reports/comprehensive_comparison/run_20260404T171802Z_alpha0.05/alpha0.05_unified_summary.csv` **(external run store, not in this repo)** | **90/90** PICP cells within 0.05 pp; the derived "Pass (>=95%)" row reproduces exactly as 12/11/4/3/13/15 | traceable, but the source is NOT in this repository (see section 2) |
| `tab:base-model-summary` (7 base learners x pass/PICP/WIS/MAE/time) | `reports/comprehensive_comparison/extended_basemodel_benchmark/extended_basemodel_leaderboard.csv` **(external run store)** | **35/35** values match (5 metrics x 7 models) | traceable, source external |
| `tab:error-by-bin` (MAE_log10 by target and temporal bin) | `data/d224_stability_max/stability_20260721T183813Z/seed_100/table9_error_by_bin.csv` | **24/24** cells (4 targets x 6 bins) | traceable, in-repo |
| `tab:error-stats` (error distribution statistics) | `data/d224_stability_max/stability_20260721T183813Z/seed_100/table10_error_statistics.csv` | **24/24** cells (4 targets x 6 statistics) | traceable, in-repo |

Note on `tab:error-by-bin`: the manuscript previously collapsed the six source bins into
four by merging the three tail bins. That masked the fact that the tail bins each contain
a single time point (n = 1) and that the tritium tail error is 1.7870, not 0.0704. The
R2 revision restores the six-bin structure and adds an `n` row.

---

## 2. Known gaps (reproducibility not yet closed)

| Artifact | Status | Detail |
|---|---|---|
| `tab:benchmark` | **source missing from this repo** | No file under `tables/` contains the labels `SC-PIML` or `SC-PIML+`; the labels used here are `proposed`, `Proposed CV+`, etc. The README's Data Availability claim that this release enables reproduction of all tables is therefore **not currently true** for the paper's central benchmark table. The source exists only in the external run store referenced above. |
| `tab:component-ablation` | untraceable | No producing source identified for the four-row ablation table. The `14/20` figure in the surrounding prose is a cross-reference to the 20-split robustness result, not a 3-split measurement, and is fine. |
| Supplementary Fig. S1 (`ncal_sensitivity`) | partial | `tables/ncal_sweep_summary.csv` is a plausible candidate but the generation chain is not established. |
| Supplementary Fig. S3 / `Fig_S13_binwise_width_comparison` | hardcoded | `scripts/make_fig_s5_binwise_width.py` plots from hardcoded arrays; no run directory or raw data is consumed. |
| mean-WIS row of `tab:benchmark` | partial | Not present in the summary CSV (which carries width, not WIS). WIS must be recomputed from the per-sample predictions. |

**Recommended fix (not yet done):** copy the two external summary files into `tables/`
under explicit names (for example `tables/benchmark_6methods_unified_summary.csv` and
`tables/base_model_leaderboard.csv`) and record them here, so the chain stays inside the
public release.

---

## 3. Exclusion list

Do **not** derive any published number from the following.

| Artifact | Reason | Superseded by |
|---|---|---|
| `reports/comprehensive_comparison/run_20260404T134324Z_alpha0.05` (external) | **Degenerate-eps contamination of the SC-PIML column.** Decay heat reports 52.20% where the paper publishes 97.8%; contact dose reports 93.39% where the paper publishes 99.8%. Its SC / CQR / ACI / EnbPI columns match the paper exactly (97.56 / 96.42 / 94.71 / 94.02 vs 97.6 / 96.4 / 94.7 / 94.0), so only the SC-PIML column is affected. The deviation tracks eps dominance: decay heat -45.6 pp (eps dominates; max 2.5 W, 91% of test points below 1 W), contact dose -6.4 pp (eps comparable at the low end), activity +0.04 pp (eps negligible). Root cause: `log10(y + 1.0)` is effectively linear-in-y for these low-magnitude aggregates, which collapses the multiplicative conformal score. | `run_20260404T171802Z_alpha0.05` |
| `reports/comprehensive_comparison/run_20260404T153155Z_alpha0.05` (external) | Superseded partial run: matches 86/90 cells. Its **SC-PIML+** column disagrees with the manuscript on four targets (Po-210, Hg-203, Pb-203, Mn-54), for example Hg-203 95.8 vs 97.7. Its other five columns are correct. | `run_20260404T171802Z_alpha0.05` |
| `tables/baseline_comparison_results_OLD_bugged_aci.csv` (in-repo) | Already marked as bugged (ACI). | `tables/baseline_comparison_results_alpha005.csv` |

---

## 4. Method-name reconciliation

The manuscript names are SC, CQR, ACI, EnbPI, SC-PIML, SC-PIML+. Intermediate files use
different labels, and one of them is ambiguous.

| Label found in CSVs | Manuscript column | Note |
|---|---|---|
| `proposed`, `Proposed CV+` | **ambiguous** | maps to SC-PIML or SC-PIML+ depending on the file; never assume |
| `proposed_cvplus_unprojected` | SC-PIML | |
| `proposed_cvplus_projaware` | SC-PIML+ | |
| `cqr`, `CQR` | CQR | |
| `aci`, `ACI` | ACI | |
| `enbpi`, `EnbPI` | EnbPI | |
| (no label present) | SC | the split-conformal baseline has no counterpart label in this repository |

Because `tables/P3_method_comparison_d224.csv` (93.4% for decay heat) and
`tables/alpha0.05_comparison.csv` (97.6%) both use a `Proposed CV+` label for decay heat
but disagree, **the label alone does not identify the manuscript column**. Always resolve
via the provenance chain in section 1.

---

## 5. How to re-verify

```bash
# error-decomposition tables (in-repo chain)
python3 - <<'PY'
import csv
for f in ("table9_error_by_bin","table10_error_statistics"):
    p=f"data/d224_stability_max/stability_20260721T183813Z/seed_100/{f}.csv"
    print(p, len(list(csv.DictReader(open(p)))), "source rows")
PY
```

The benchmark and base-model chains require the external run store; the comparison
procedure is recorded in the project audit report `E_provenance_audit.md`.
