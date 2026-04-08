# Open Data: Time-aware Conformal Uncertainty Quantification for Radioactivity Time-Series Prediction

This repository contains the data and figure-generation scripts accompanying the paper:

> **Time-aware conformal uncertainty quantification with physics-consistent projection for radioactivity-related time-series prediction**
>
> Xiaokang Zhang et al.
>
> Submitted to *Machine Learning: Science and Technology* (MLST)

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
│   └── d224_alltargets/          # Full 15-target run (9 nuclides + 6 derived)
├── scripts/                      # Python scripts to regenerate paper figures
│   ├── make_fig2_actual_evolution.py      # Fig 2: temporal evolution curves
│   ├── make_figs3_4_per_time_bin.py       # Figs 3–4: per-time-step & per-bin calibration
│   ├── make_figs3_4_9_regenerate.py       # Figs 3, 4, 9: alternative generation + method comparison
│   └── make_figs5_8_ood_analysis.py       # Figs 5–8: projection ablation & OOD analysis
├── figures/                      # Generated figures (PDF)
│   ├── method_overview_v2.tex    # TikZ source for Fig 1 (method overview)
│   └── Fig_P3_*.pdf              # Figures 1–9 as used in the paper
├── tables/                       # CSV tables referenced in the paper
├── LICENSE
├── requirements.txt
└── README.md
```

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

### Reference data (`data/reference/`)

- **cell3_true_full.csv** — FISPACT-II reference inventory for cell 3 (used for true evolution curves in Fig 2).

## Citation

If you use this data or code, please cite the paper:

```bibtex
@article{zhang2025timeaware,
  title={Time-aware conformal uncertainty quantification with physics-consistent
         projection for radioactivity-related time-series prediction},
  author={Zhang, Xiaokang and others},
  journal={Machine Learning: Science and Technology},
  year={2025},
  note={Submitted}
}
```

## License

Code: MIT License. Data: CC-BY-4.0.
