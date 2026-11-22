#  ECG Denoising with Least Squares + Extended Kalman Filter

**Problem.** Ambulatory and wearable ECGs are corrupted by muscle activity (EMG) and baseline
wander from movement/breathing. Muscle noise overlaps the ECG spectrum, so a fixed band-pass
filter cannot remove it without also distorting the P-QRS-T waves.

**Idea.** Use what we *know* about the ECG: every beat has the same P-Q-R-S-T shape.
1. **Least Squares (LS):** fit a 5-Gaussian beat model (McSharry et al., 2003) to the
   phase-averaged noisy beat with a hand-written Levenberg–Marquardt algorithm.
2. **Maximum Likelihood (MLE):** estimate the measurement-noise variance as the residual
   mean square (the ML estimate for Gaussian noise) → Kalman `R`.
3. **Extended Kalman Filter (EKF):** state = [cardiac phase θ, ECG amplitude z], nonlinear
   dynamical ECG model as the process model (Sameni et al., 2007).
4. **EKS (Rauch–Tung–Striebel smoother):** backward pass for offline recordings.

Compared against a **0.5–40 Hz Butterworth band-pass** and **wavelet denoising** (sym4, Bayes).

---

## How to run (in this order)

Open MATLAB in this folder, then:

| Step | Script | What it does 
|---|---|---|---|
| 1 | `prepare_dataset` | Reads raw PhysioNet files in `data/raw/` → `data/ecg_dataset.mat` 
| 2 | `tune_parameters` | Grid-search of EKF process noise on **TRAIN** data → `results/tuned_params.mat` 
| 3 | `run_experiments` | **TEST**: SNR sweep, demo case, real record 105 → `results/experiment_results.mat`, `results/tables/snr_sweep_results.csv` 
| 4 | `make_figures` | 10 figures, 1920×1080 PNG → `results/figures/` | ~30 s |
| 5 | `make_videos` | 2 MP4 clips, 1920×1080 @ 30 fps → `results/videos/` (resumable: frames cached in `frames_*`) 
| – | `demo_quick` | One noisy ECG, prints metrics and plots on screen (live demo; change noise type/SNR at the top)

`main.m` runs steps 1–5 in sequence. `data/ecg_dataset.mat` is already included, so step 1
is only needed if you want to rebuild it from the raw files. Step 2 is optional: its result is
already the default in `src/ekf_ecg.m`.

**Requirements:** MATLAB R2023a+, Signal Processing Toolbox, Wavelet Toolbox (wavelet baseline only).

**If MATLAB crashes** ("heap corruption") during long loops: start MATLAB with software rendering
(`matlab -softwareopengl`) and re-run the script. `make_videos` resumes from the frames it already saved.
If a PNG is open in a viewer and locked, `grab_frame` saves `<name>_new.png` instead.

**Tuned EKF hyper-parameters** (train data only): qa = 0.05, qθ = 0.02, qη = 0.001 (+ fixed qb = 0.05, qω = 0.05).

---

## Data

All from PhysioNet (360 Hz, mV): MIT-BIH Arrhythmia DB and MIT-BIH Noise Stress Test DB (NSTDB).

| File | Content | Used as |
|---|---|---|
| `100`, `103` | Clean ECGs, two different patients (lead MLII) | Ground truth |
| `ma` | Real muscle-artifact noise (electrodes placed so no ECG is visible) | Noise |
| `bw` | Real baseline-wander noise | Noise |
| `105` | ECG with strong *natural* noise, no clean version | Qualitative real-world check |

The clean ECGs contain **no** added noise. Noise is added by us, so the true clean signal is
known and SNR/RMSE can be computed exactly. Noise is scaled with **PhysioNet's `nst`
definition** of SNR:
`S = (trimmed-mean QRS peak-to-peak)^2 / 8`, `N = (trimmed-mean 1-s RMS)^2`,
`a = sqrt(S / (N·10^(SNR/10)))`. For **MA+BW** the noise power budget is split equally between
the two sources.

### Train / test split (`src/data_split.m`)

Only the EKF process-noise hyper-parameters are tuned. Everything else is computed from each
noisy recording itself.

| | Clean ECG | Noise segment |
|---|---|---|
| **Train** (tuning only) | Record 100, 1:00–2:00 | NSTDB starting at sample 400000 |
| **Test** | Record 100, 5:00–6:00 **and** record 103, 5:00–6:00 | 5 random segments (seed 2026), none overlapping the training segment |

Test conditions: noise ∈ {MA, BW, MA+BW} × SNR ∈ {−5, −2.5, …, 20} dB × 2 records × 5 noise
segments → **10 runs per point**.

The LS wave parameters and the ML noise variance are **not** trained. They are re-estimated from
every noisy test signal (no clean signal or labels used), so nothing leaks from the test set.

---

## Folder layout

```
main.m  prepare_dataset.m  tune_parameters.m  run_experiments.m
make_figures.m  make_videos.m  demo_quick.m
src/        functions (WFDB readers, nst scaling, R-peak detector, LS/LM fit, EKF+EKS, metrics)
data/raw/   original PhysioNet files      data/ecg_dataset.mat   converted dataset
results/    figures/  videos/  tables/  experiment_results.mat  tuned_params.mat
slides/     presentation slides and speaker script
```

## Key functions

| File | Purpose |
|---|---|
| `src/fit_gaussians_lm.m` | **LS**: Levenberg–Marquardt fit of the 5-Gaussian PQRST model (box-constrained) |
| `src/ekf_pipeline.m` | Full chain + **MLE** of noise variance |
| `src/ekf_ecg.m` | **EKF** forward pass + **RTS smoother**, analytic Jacobians |
| `src/nst_power.m`, `src/make_noise.m` | PhysioNet `nst` SNR calibration |
| `src/detect_rpeaks.m`, `src/ecg_phase.m` | R peaks → cardiac phase observation |
| `src/eval_metrics.m` | SNR in/out/improvement, RMSE, PRD, correlation |

## References
- P. E. McSharry et al., "A dynamical model for generating synthetic ECG signals," *IEEE TBME*, 2003.
- R. Sameni et al., "A nonlinear Bayesian filtering framework for ECG denoising," *IEEE TBME*, 2007.
- G. B. Moody, W. Muldrow, R. G. Mark, "A noise stress test for arrhythmia detectors," *Computers in Cardiology*, 1984.
- Goldberger et al., "PhysioBank, PhysioToolkit, and PhysioNet," *Circulation*, 2000.
