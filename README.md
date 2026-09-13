# Open Data: Time-aware Conformal Uncertainty Quantification for Fusion Blanket Activation Time-Series Prediction

This repository contains the data and figure-generation scripts accompanying the paper:

> **Time-aware conformal uncertainty quantification with physics-consistent projection for fusion blanket activation time-series prediction**
>
> Xiaokang Zhang, Xilong Tong, Yanshi Wei and Shanliang Zheng
>
> *Machine Learning: Science and Technology* (MLST), manuscript MLST-105507 (revised version R1)

The data and scripts in this release correspond to the **R1 revision** (August 2026), which
addresses all reviewer comments from the first round of review. Key additions over the
initial submission include: bin-wise simultaneous coverage guarantees via geometry-level
score aggregation, a 20-split robustness study (`tables/split_robustness.csv`), a component
ablation study (`data/ablation/`, `tables/P3_projection_ablation_H3_summary.csv`), a point
prediction error decomposition (`scripts/analyze_error_decomposition.py`,
`data/d224_stability_max/`), and per-bin interval-width comparisons between aggregation modes
(`scripts/make_fig_s5_binwise_width.py`, Figure S13).

## Repository structure

```
├── data/                         # Processed data used by figure scripts
│   ├── dataset/                  # D224 dataset (labels + train/cal/test split)
│   ├── reference/                # FISPACT-II reference inventory (cell 3)
│   ├── d224_runs/                # SC-PIML+ calibration variant runs on D224
│   │   ├── d224_cvplus_smoke/    # Time-binned CV+ (H3)
│   │   ├── d224_cvplus_global_*/  # Global calibration baselines
│   │   ├── d224_cvplus_phase_timebin_*/  # Phase × time-bin calibration
│   │   ├── d224_cvplus_decay_heat/  # Decay heat target
│   │   └── d224_cqr_timebin_*/   # CQR baselines
│   ├── ablation/                 # Projection ablation experiment (3 modes)
│   │   ├── run__none/            # No projection
│   │   ├── run__after/           # Post-hoc projection
│   │   └── run__aware/           # Projection-aware (proposed)
│   ├── n70_valid/                # Out-of-distribution validation (N70)
│   │   ├── baseline/             # Standard pipeline
│   │   └── candidate/            # SC-PIML+ with acceptance gate
│   ├── d224_alltargets/          # Full 15-target run (9 nuclides + 6 derived)
│   └── d224_stability_max/       # 20-split stability run (seed 100 shown; error-decomposition inputs)
├── scripts/                      # Python scripts to regenerate paper figures
│   ├── make_fig2_actual_evolution.py      # Fig 2: temporal evolution curves
│   ├── make_figs3_4_per_time_bin.py       # Figs 3–4: per-time-step & per-bin calibration
│   ├── make_figs3_4_9_regenerate.py       # Figs 3, 4, 9: alternative generation + method comparison
│   ├── make_figs5_8_ood_analysis.py       # Figs 5–8: projection ablation & OOD analysis
│   ├── analyze_error_decomposition.py     # R1: Tables 8–9, error decomposition (seed 100)
│   ├── make_fig_s5_binwise_width.py       # R1: Figure S13, binwise width comparison
│   ├── make_p3_baseline_comparison_figure.py  # R1: baseline comparison figure
│   └── generate_base_model_comparison_fig.py  # R1: base-learner comparison figure
├── figures/                      # Generated figures (PDF)
│   ├── method_overview_v2.tex    # TikZ source for Fig 1 (method overview)
│   └── Fig_P3_*.pdf              # Figures 1–9 as used in the paper
├── tables/                       # CSV tables referenced in the paper
├── LICENSE
├── requirements.txt
├── PROVENANCE.md                 # Artifact-to-source manifest + exclusion list
└── README.md
```

## Data provenance

`PROVENANCE.md` records, for every numeric artifact in the manuscript, the file that
produced it, together with an **exclusion list** of artifacts that must not be used to
derive published numbers. It currently flags one run whose SC-PIML column is contaminated
by a degenerate log-transform on low-magnitude aggregates, and one superseded partial run.
It also records the known reproducibility gaps in this release. Consult it before reusing
any CSV in `tables/`.

## Reproducing the figures

```bash
pip install -r requirements.txt

# Fig 2: Temporal evolution curves (5×3 grid)
python scripts/make_fig2_actual_evolution.py

# Figs 3–4: Per-time-step and per-bin calibration diagnostics
python scripts/make_figs3_4_per_time_bin.py

# Figs 5–8: Projection ablation & OOD failure-mode analysis
python scripts/make_figs5_8_ood_analysis.py

# Fig 9: Method comparison bar chart
python scripts/make_figs3_4_9_regenerate.py

# R1 Tables 8–9: Point prediction error decomposition (requires data/d224_stability_max/)
python scripts/analyze_error_decomposition.py

# R1 Figure S13: Per-bin width comparison (no-aggregation vs max-aggregation)
python scripts/make_fig_s5_binwise_width.py
```

Fig 1 (method overview) is a standalone TikZ diagram; compile with:
```bash
cd figures && pdflatex method_overview_v2.tex
```

## Data description

### D224 dataset (`data/dataset/`)

- **labels_long_bzsum_elapsed_total_act9_with_yref_logratio.csv** — Merged inventory records for 224 CFETR COOL blanket geometries, covering 15 targets: ³H, ⁶Li (log-ratio), ⁷Li (log-ratio), decay heat, activity, contact dose, and 9 activation products (²¹⁰Po, ²¹⁰Bi, ²⁰³Hg, ²⁰⁴Tl, ²⁰⁷Bi, ²⁰³Pb, ⁵⁴Mn, ⁶⁰Co, ⁵⁵Fe).
- **splits_d224_seed42.json** — Canonical train/calibration/test split (149/31/44 geometries, seed 42).

### SC-PIML+ run outputs (`data/d224_runs/`, `data/ablation/`, `data/n70_valid/`)

Each run directory contains:
- `metrics.json` — Global metrics (PICP, WIS, RMSE, etc.)
- `metrics_per_time_*.csv` — Per-time-step coverage and interval width
- `predictions_test_*.csv` — Test set predictions with uncertainty bounds
- `timeaware_uq_config.json` — Configuration used for the run (where available)
- `acceptance_gate_*.json` — OOD acceptance gate results (N70 candidate only)

### 20-split stability inputs (`data/d224_stability_max/stability_20260721T183813Z/seed_100/`)

- `metrics_per_time_*.csv` — Per-time-step metrics for the four representative targets (³H, ⁶Li, ⁷Li, decay heat) used by `analyze_error_decomposition.py` to produce Tables 8–9.
- `table9_error_by_bin.csv`, `table10_error_statistics.csv` — Frozen outputs reported in the manuscript.

### Reference data (`data/reference/`)

- **cell3_true_full.csv** — FISPACT-II reference inventory for cell 3 (used for true evolution curves in Fig 2).

## Data generation workflow

The dataset is generated with a coupled neutronics–activation workflow:
1. **OpenMC** (with the ENDF/B-VIII.0 nuclear data library) performs neutron transport on a standardized CFETR COOL PbLi blanket sector, producing per-geometry neutron flux spectra.
2. **FISPACT-II v5.0** (with the EAF-2010 cross-section library) solves the Bateman equations for transmutation and decay chains, yielding time-dependent inventories over a deterministic 122-point time grid (irradiation 0–10 y; cooling 10–10 010 y).

## Citation

If you use this data or code, please cite the paper:

```bibtex
@article{zhang2026timeaware,
  title={Time-aware conformal uncertainty quantification with physics-consistent
         projection for fusion blanket activation time-series prediction},
  author={Zhang, Xiaokang and Tong, Xilong and Wei, Yanshi and Zheng, Shanliang},
  journal={Machine Learning: Science and Technology},
  year={2026},
  note={MLST-105507, revised version}
}
```

## License

Code: MIT License. Data: CC-BY-4.0.
