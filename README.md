# ML Pipeline for Topographical Feature Extraction from Turbulent Shear Flow LES Data

**Author:** Tennisraja Vijayaram  
**Programme:** MSc Artificial Intelligence — Aston University  
**Supervisor:** Dr Andrew McMullan  
**Submission:** September 2026

---

## Overview

This repository contains the full MATLAB codebase, raw dataset, and output figures for the MSc dissertation:

> *"Machine Learning Pipeline for Topographical Feature Extraction from Turbulent Shear Flow Large Eddy Simulation Data"*

The pipeline automatically identifies and classifies coherent vortex structures (Kelvin–Helmholtz roll-up vortices) from LES flow field data across isothermal and reactive combustion cases (equivalence ratios φ = 0.125–8.0).

The key output is the **structure spacing parameter k = mean(l/xₘ) / λ**, validated against the Browand & Troutt (1985) theoretical range of **k ∈ [0.46, 0.60]**.

---

## Repository Structure

```
EP4DIS-ML-Pipeline-Dissertation/
│
├── code/
│   ├── main_ml_pipeline.m          Main pipeline: load → clean → compare → validate (isothermal)
│   ├── MLstructureFinder.m         Core ML structure detection engine (all 8 cases)
│   ├── statistical_baseline.m      Statistical threshold cleaning for all 8 cases
│   ├── reactive_analysis.m         Statistical analysis for 4 reactive cases
│   ├── ml_comparison.m             4-algorithm ML comparison (RF, SVM, k-NN, DT)
│   ├── ml_comparison1.m            Extended ML comparison variant
│   ├── figures.m                   Additional figure generation utilities
│   └── tripletStatsLoaderGUI.m     GUI tool for loading and inspecting triplet data
│
├── data/
│   ├── isothermal_all_triplets.csv     Isothermal case — raw vortex triplet data
│   ├── hs1f-hs1h-all_triplets.csv      φ = 1.0  — Stoichiometric reactive case
│   ├── hs1f-ls8h_all_triplets.csv      φ = 8.0  — Very fuel-rich reactive case
│   ├── hs8h-ls1f_all_triplets.csv      φ = 0.125 — Very lean reactive case
│   ├── labelled.csv                    Training data: 12 features + manual label (1/0)
│   ├── manual_labels.csv               Manual labels: 3 geometric features + label
│   └── manual_label_sample.csv         Sample of manual labelling output
│
├── figures/                            All output figures (numbered, descriptive names)
│   └── (see Figures section below)
│
├── .gitignore
└── README.md
```

> **Note:** Raw data for φ = 0.25, 0.5, 2.0, and 4.0 are not included (original LES files exceeded storage limits). Their derived `_clean.csv` outputs exist but are excluded from this repo. The four cases above are sufficient to reproduce all dissertation results.

---

## Dataset Description

Each `*_all_triplets.csv` contains vortex triplet detections from the LES flow field. Each row is one detected triplet:

| Column | Description |
|---|---|
| `l` | Physical structure spacing (m) |
| `l_over_xminusx0` | Normalised spacing l / (xₘ − x₀) |
| `b_over_a` | Asymmetry ratio — b/a = 1 for a symmetric vortex |
| `x_upstream`, `x_downstream` | Upstream / downstream vortex x-positions (m) |
| `detJ_*`, `traceJ_*`, `discJ_*` | Velocity gradient Jacobian invariants at upstream, core, downstream |
| `label` | Manual classification: 1 = genuine structure, 0 = spurious |

---

## How to Run

Set the MATLAB working directory to the **root of this repository**, then run scripts in this order:

### Step 1 — Statistical baseline (all 8 cases)
```matlab
run('code/statistical_baseline.m')
```
Outputs: `*_full_clean.csv`, `*_filt_clean.csv`, Figures 11–16.

### Step 2 — Main ML pipeline (isothermal case)
```matlab
run('code/main_ml_pipeline.m')
```
Outputs: `stat_clean_triplets.csv`, Figures 1–6. Requires `labelled.csv` and `manual_labels.csv` in `data/`.

### Step 3 — 4-algorithm ML comparison
```matlab
run('code/ml_comparison.m')
```
Outputs: Figures 17–18. Requires `labelled.csv` in `data/`.

### Step 4 — Reactive cases (4 cases)
```matlab
run('code/reactive_analysis.m')
```
Outputs: Figures 7–10.

### Step 5 — Full per-case ML analysis (all available cases)
```matlab
run('code/MLstructureFinder.m')
```
Outputs: Per-case histograms and contour PDFs (Figures 19–48).

---

## ML Methods

| Model | MATLAB Function | Validation |
|---|---|---|
| Random Forest | `TreeBagger(100, X, Y)` | OOB error |
| Support Vector Machine | `fitcsvm(..., 'KernelFunction', 'rbf')` | 5-fold CV |
| k-Nearest Neighbours | `fitcknn(..., 'NumNeighbors', 5)` | 5-fold CV |
| Decision Tree | `fitctree(X, Y)` | 5-fold CV |

**Features (12 total):** `l`, `l/(xₘ−x₀)`, `b/a`, and 9 Jacobian invariants (det J, trace J, discriminant of J at upstream, core, and downstream vortex positions).

---

## Figures

| File | Description |
|---|---|
| `01_isothermal_lx_pdf_comparison.png` | l/(x−x₀) PDF: Original vs Statistical vs ML cleaned |
| `02_isothermal_ba_pdf_comparison.png` | b/a PDF: Original vs Statistical vs ML cleaned |
| `03_isothermal_joint_pdf_comparison.png` | Joint PDF (l/x vs b/a): three cleaning methods |
| `04_isothermal_downstream_trends.png` | Downstream trends: l/x, b/a, structure size |
| `05_rf_feature_importance.png` | Random Forest OOB feature importance |
| `06_isothermal_fourway_cleaning_comparison.png` | 4-way: Raw / Stat / ML-rule / ML-manual |
| `07_reactive_4cases_k_values.png` | k values — 4 reactive cases |
| `08_reactive_4cases_ba_values.png` | Mean b/a — 4 reactive cases |
| `09_reactive_4cases_lx_pdfs.png` | l/(x−x₀) PDFs — 4 cases side by side |
| `10_reactive_4cases_ba_pdfs.png` | b/a PDFs — 4 cases side by side |
| `11_all8cases_k_full_domain.png` | k — all 8 cases, full domain |
| `12_all8cases_ba_full_domain.png` | Mean b/a — all 8 cases, full domain |
| `13_all8cases_k_filtered_region.png` | k — all 8 cases, 0 < x < 0.16 m |
| `14_all8cases_ba_filtered_region.png` | Mean b/a — all 8 cases, 0 < x < 0.16 m |
| `15_all8cases_lx_pdf_filtered.png` | l/(x−x₀) PDFs — all 8 cases, filtered region |
| `16_all8cases_ba_pdf_filtered.png` | b/a PDFs — all 8 cases, filtered region |
| `17_ml_algorithm_accuracy_comparison.png` | Accuracy: RF vs SVM vs k-NN vs Decision Tree |
| `18_ml_k_value_by_model.png` | k value produced by each ML model |
| `19–26_hist_lx_*.png` | l/(x−x₀) histograms per case (isothermal → φ=8.0) |
| `27–34_hist_ba_*.png` | b/a histograms per case (isothermal → φ=8.0) |
| `35–42_contour_pdf_*.png` | Joint contour PDFs per case |
| `43–45_k_vs_phi_region*.png` | k vs φ by spatial region (regions 1–3) |
| `46_k_all_cases_summary.png` | k summary — all cases |
| `47_accuracy_all_cases.png` | Classification accuracy — all cases |
| `48_k_vs_phi_all_cases.png` | k vs φ — all cases combined |

---

## Key Parameters

| Parameter | Value | Source |
|---|---|---|
| λ (mixing layer thickness) | 0.428 m | McMullan et al. (2015) |
| x₀ (virtual origin) | 0.0 m | Set at splitter plate |
| n_std (cleaning threshold) | ±1σ | Statistical baseline |
| Spatial filter | 0 < x < 0.16 m | Turbulent development region |
| Theoretical k range | [0.46, 0.60] | Browand & Troutt (1985) |

---

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox

---

## Citation

> Vijayaram, T. (2026). *Machine Learning Pipeline for Topographical Feature Extraction from Turbulent Shear Flow Large Eddy Simulation Data.* MSc Dissertation, Aston University.

---

## Acknowledgements

Supervisory support from Dr Andrew McMullan (Aston University).  
LES simulation data provided by the Aston University Fluid Mechanics research group.
