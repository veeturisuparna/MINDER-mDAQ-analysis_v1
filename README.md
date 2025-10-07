# MINDER system initial signal analysis pipeline — README

* **`ECG_Analysis.m`** — compares dry vs. wet (gel) ECG leads, extracts beat‑level Signal Quality Indices (SQIs), and renders per‑subject and roll‑up figures.
* **`main_EDA_Analysis.m`** — loads DC‑EDA (custom device) and BIOPAC EDA, aligns and cleans both, decomposes **SCL/SCR**, and computes validation metrics (corr/DTW/coherence, ΔZC).

> The EDA workflow mirrors the study protocol and metrics from the mEDA validation paper; the ECG workflow summarizes the arm‑ECG dry‑electrode evaluation (silver‑knit vs. gold‑plated) and bioinstrumentation comparison.

---

## Why this code exists

* **ECG**: to rapidly quantify and visualize how **dry electrodes** (esp. silver‑knit) behave vs. gel across postures/movements and bioinstrumentation choices.
* **EDA**: to align a **custom DC‑EDA** device with a research‑grade BIOPAC system using simple, reproducible metrics over SCL/SCR.


## Requirements

* **MATLAB** (R2020b+ recommended)
* **Signal Processing Toolbox** (`findpeaks`, `pwelch`, `filtfilt`, `butter`, `mscohere`, `xcorr` , `dtw` 
---

## Data layout & assumptions
### ECG (`ECG_Analysis.m`) for ETextiles dataset paper with mDAQ & BIOPAC
* Input: `.mat` files under a root `path`, each containing `data` with:
  * **col 1 = wet (gel) lead**, **col 2 = dry lead**
  * Sampling rate `srate = 125` Hz (set in script, needs to be adjusted for UMASS dataset to 100 Hz)
* Script iterates all `.mat` files recursively for BIOPAC & 'tsv' or 'csv' file for mDAQ.

### EDA (`main_EDA_Analysis.m`)
* **DC‑EDA CSV**: 
  * `VOLTAGE_MV` + `ISI_US` **or** `eda` + (`isi` **or** `ISI_US`) **or** `EDA_uS` + `ISI_US`
* **BIOPAC MAT**: EDA signal in **column 5**, sampled at **1000 Hz**
* Filenames must match `config.participant`, e.g. `P2.csv` and `P2 (BIOPAC_only_All_Sensors_Ch1Sync).mat`
---

## Quick start
### ECG
```matlab
% 1) Set the input folder
Example: 
envpath = 'G:\My Drive\BSN internal testing\ECG+EDA testing\Biopac\4-11-24_electrode_study_ECG';
path = envpath;  % in script
% 2) Run the script
ECG_Analysis
```

**What you’ll get**: per‑subject full‑screen figure (lead traces + detected beats, beat ensembles, SQI boxplots, histograms) and a final **Sorted SQI** dashboard across subjects.

### EDA

```matlab
% 1) Configure participant and data folder in main_EDA_Analysis.m
config.participant = 'P2';
config.dataPath    = '.';  % folder containing P2.csv + BIOPAC .mat
% 2) Run
main_EDA_Analysis
```

**What you’ll get**: three figures — **Signal Overview** (raw/cleaned, filtered, SCL/SCR + metrics text), **Event‑Based Analysis** (means ± SD, % change from baseline), **Decomposition** (SCL/SCR per device).

---

## What the code actually does
### ECG pipeline (dry vs. wet)
1. **Normalize** each lead (125 Hz) and remove DC offset.
2. **HilbertTransform** → analytic signal; peak‑based **IdentifyQRS** locates beats.
3. **Epoch** each beat around the R‑peak; analyze the QRS window (samples 13–39).
4. **SQIs** per subject & lead:
   * **Skewness** (morphology asymmetry)
   * **Kurtosis** (peakedness)
   * **Relative PSD 5–15 Hz / 0–45 Hz** via `pwelch` (QRS band energy)
5. **Plots**: traces with peaks, beat overlays, boxplots (Dry vs. Wet) with reference lines, histogram vs. normal; final cross‑subject “Sorted SQI”.

### EDA pipeline (DC‑EDA vs. BIOPAC)
1. **Load** DC‑EDA & BIOPAC; build time (minutes). BIOPAC assumed 1 kHz; DC device ~100 Hz.
2. **Resample** both to a common grid (default **100 Hz**), pad/truncate to equal length.
3. **Artifact handling** (z‑score domain): sliding‑window std + slope checks → set artifacts to NaN → interpolate → **de‑normalize** to original scale.
4. **Filter**: low‑pass 1.5 Hz (smoothing). (A 10 Hz low‑pass is applied pre‑normalization for visualization.)
5. **Decompose SCL/SCR** via FFT mask at `config.sclCutoffHz` (default **0.05 Hz**; adjust if replicating alternate splits).
6. **Normalize for events**: min–max on the 1.5 Hz filtered signal.
7. **Metrics**
   * **SCL**: max cross‑correlation & lag (minutes), **mean coherence** (≤ `coherence_freq_max`, default 0.05 Hz), **windowed DTW** (downsampled)
   * **SCR**: baseline‑centered **zero‑crossing** counts with refractory period; report **ΔZC = |ZC_BIOPAC − ZC_DC|**
8. **Event analysis**: mean ± SD per labeled interval; % change from the baseline event.

---

## Outputs & interpretation
* **ECG**: For each subject, compare Dry vs. Wet on **Skewness**, **Kurtosis**, **Rel. PSD (5–15 Hz)**.
  * Higher QRS band energy and stable skew/kurtosis → cleaner morphology (often seen with silver‑knit in static tasks).
  * Use "Sorted SQI" to rank subjects/sessions and spot outliers (e.g., movement‑noisy trials).
* **EDA**: Report string in the overview figure includes **Corr**, **Lag (min)**, **Mean Coherence**, **Norm. DTW**; SCR panel shows ZC markers.
* **ΔZC ** : Difference between SCR's from mDAQ & BIOPAC.
---

## Configuration you’ll need to change

```matlab
% EDA sampling / filters
config.resampleRateHz = 100;
config.filterCutoffHz.low  = 1.5;    % smoothing
config.filterCutoffHz.high = 10;     % pre‑norm viz
config.sclCutoffHz         = 0.05;   % SCL/SCR split (set higher if needed)
% Artifact handling (z‑score domain)
config.artifact.sd_window_sec = 60;
config.artifact.sd_threshold  = 0.01;
config.artifact.slope_window_sec = 5;
config.artifact.absolute_slope_tolerance = 1e-4;
% SCL metrics
config.analysis.coherence_freq_max = 0.05;
config.analysis.dtw_window_min     = 0.1;
config.analysis.dtw_overlap_min    = 0;
config.analysis.dtw_downsample_factor = 10;
% SCR metrics
config.analysis.scr_cross_window   = 0.5;  % sec
config.analysis.scr_min_zc_dist_sec= 0.5;  % sec
```
---
