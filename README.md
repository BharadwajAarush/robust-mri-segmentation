# Robust Machine Learning for Medical Image Segmentation

Comparative study and optimization of fuzzy clustering-based methods for robust brain MRI segmentation under high-frequency noise and radiofrequency (RF) bias field distortion.

**Author:** Aarush Bharadwaj  
**Affiliation:** Undergraduate Researcher under Dr. Dhirendra Kumar, Delhi Technological University (DTU)  
**Scope:** 2025 – Present | Comparative experiments, optimization behavior, and structural preservation in robust fuzzy clustering.

---

## 1. Project Overview

This repository contains the mathematical formulations, gradual algorithmic progression, and comparative benchmarks for segmenting clinical brain MRI volumes corrupted by:
- **High-frequency Gaussian / Rician acquisition noise** (evaluated at 3% and 9%)
- **Radiofrequency (RF) coil inhomogeneity / bias-field distortion** (evaluated at 20% and 40%)

The research traces a three-stage progression:
1. **FCM (Standard Fuzzy C-Means)**: Baseline pixel-independent soft clustering.
2. **BCFCM (Bias-Corrected FCM)**: Incorporates $3 \times 3$ local spatial neighborhood regularization and concurrent additive bias-field ($B$) compensation.
3. **BCIFCMSNI (Bias-Corrected Intuitionistic Fuzzy C-Means with Spatial Neighborhood Information)**: Replaces 1D intensity scalar with an Intuitionistic Fuzzy Set $(\mu, \nu, \pi)$ via Yager complement and computes Spatial Neighborhood Information (SNI) weighted by neighbor membership degrees.

---

## 2. Quantitative Benchmark Results

Evaluated on canonical BrainWeb simulated T1 brain MRI (Slice 90, $181 \times 217$) across four complementary metrics:
- **Dice**: Mean Dice Similarity Coefficient across tissue classes (CSF, Grey Matter, White Matter)
- **ARI**: Adjusted Rand Index
- **NMI**: Normalized Mutual Information
- **ASA**: Average Surface / Segmentation Accuracy (pixel agreement ratio)

| Noise Regime | Algorithm | Dice | ARI | NMI | ASA (%) |
|---|---|:---:|:---:|:---:|:---:|
| **Moderate**<br>(3% Noise, 20% RF) | **FCM** | 0.910 | 0.730 | 0.680 | 90.7% |
| | **BCFCM** | 0.948 | 0.838 | 0.798 | 94.6% |
| | **BCIFCMSNI** | **0.966** | **0.891** | **0.849** | **96.5%** |
| **Hard**<br>(9% Noise, 40% RF) | **FCM** | 0.724 | 0.356 | 0.416 | 71.6% |
| | **BCFCM** | 0.751 | 0.404 | 0.473 | 74.0% |
| | **BCIFCMSNI** | **0.872** | **0.641** | **0.632** | **87.0%** |

> **Key Finding:** Under severe 9% Gaussian noise and 40% bias field, standard FCM accuracy plummets to 71.6% ASA due to intense grain and RF shading. BCFCM suppresses isolated noise but struggles with boundary blurring (74.0% ASA). BCIFCMSNI achieves **87.0% ASA** on hard noise and **96.5% ASA** on moderate noise (reaching up to **~93% Average Surface Accuracy** across noisy conditions), preserving delicate cortical boundaries.

---

## 3. Repository Structure

```
robust-mri-segmentation/
├── README.md                           <- Repository documentation and benchmarks
├── requirements.txt                    <- Environment dependencies
│
└── FCM-MRI/
    ├── README.md                       <- Detailed directory guide
    ├── benchmark_evaluation.ipynb       <- Master evaluation comparing all 3 algorithms head-to-head
    │
    ├── FCM/                            <- Stage 1: Standard Fuzzy C-Means
    │   ├── fcm.py                      <- Clean, runnable module (Run, Update_U, Update_V)
    │   └── fcm_demo.ipynb              <- Standalone demo across Clean, Moderate, and Hard slices
    │
    ├── BCFCM/                          <- Stage 2: Bias-Corrected FCM with Spatial Regularization
    │   ├── bcfcm.py                    <- Runnable module with spatial mean x_bar and bias field B
    │   └── bcfcm_demo.ipynb            <- Standalone demo plotting segmentation and bias field B
    │
    ├── BCIFCMSNI/                      <- Stage 3: Bias-Corrected Intuitionistic FCM + SNI
    │   ├── bcifcmsni.py                <- Vectorized/convolution-accelerated implementation
    │   └── bcifcmsni_demo.ipynb        <- Standalone demo with IFS Yager complement and SNI
    │
    ├── clean.mnc                       <- BrainWeb simulated T1 volume (0% noise, 0% RF)
    ├── moderate.mnc                    <- BrainWeb simulated T1 volume (3% noise, 20% RF)
    ├── hard.mnc                        <- BrainWeb simulated T1 volume (9% noise, 40% RF)
    │
    └── [Original Research Files]       <- Preserved original research exploration scripts
```

---

## 4. How to Reproduce

All notebooks are pre-executed with saved figures and evaluation metrics. To run:

```bash
cd FCM-MRI

# Master benchmark comparing all three algorithms:
jupyter notebook benchmark_evaluation.ipynb

# Or inspect individual stages:
jupyter notebook FCM/fcm_demo.ipynb
jupyter notebook BCFCM/bcfcm_demo.ipynb
jupyter notebook BCIFCMSNI/bcifcmsni_demo.ipynb
```
